# Render only accepted quantities; never fit a model or select favorable rows.
args<-commandArgs(TRUE)
if(length(args)!=5L || !args[1]%in%as.character(2:8) || !args[2]%in%c('JTIM','WJG') || !args[3]%in%c('complete','eligible-panels'))
 stop('Usage: render_dual_journal_figures.R FIGURE_NUMBER JTIM|WJG complete|eligible-panels INPUT_MANIFEST.tsv OUTPUT_DIRECTORY')
source('workflow/R/common/load_analysis.R');load_public_analysis()
cx<-public_analysis_context('current_figures',args[4],args[5],create=FALSE)
cx$resource_check();runtime<-cx$runtime('rendering')
n<-as.integer(args[1]);journal<-args[2];render_mode<-args[3]
dj_spec<-jsonlite::read_json(cx$input('publication/figure_spec.json'))
mapping<-cx$read_tsv('publication/PANEL_MAPPING.tsv',colClasses='character');mapping<-mapping[mapping$figure==as.character(n),]
catalog<-cx$read_tsv('publication/figure_inputs.tsv',colClasses='character');catalog<-catalog[catalog$figure==as.character(n),]
if(render_mode=='complete')invisible(cx$input(catalog$logical_input))
included<-if(render_mode=='complete')mapping$panel else mapping$panel[mapping$bundled_mode=='COMPLETE_PANEL']
stopifnot(length(included)>0L)
if(!nzchar(Sys.which('rsvg-convert')))stop('Required existing tool rsvg-convert is missing; see rendering runtime specification.')
style<-dj_spec$journals[[journal]]
fonts<-systemfonts::match_fonts(rep(style$font,2),weight=c('normal','bold'))$path
if(!all(grepl(style$font,fonts,ignore.case=TRUE)))stop('Required journal font unavailable: ',style$font)
source('workflow/R/methods/dual_journal_figures.R')
assign('.F',new.env(parent=emptyenv()),envir=.GlobalEnv)
for(k in c('panels','regions','sources','marks','axes','titles','legends','colorbars','zero_lines','text','cache','operations')).F[[k]]<-list()
.F$n<-n;.F$height<-dj_spec$height_mm[[as.character(n)]];.F$style<-style;.F$included<-included;.F$sources<-dj_sources(n)
blue<-'#0072B2';orange<-'#D55E00';negative<-blue;positive<-orange;individual<-blue;adjacent<-blue;connector<-'#C6C6C6';unavailable<-'#D6D6D6';light<-'#D6E8F3'
stem<-paste0('Figure',n,'_',journal,if(length(included)<nrow(mapping))'_eligible-panels' else '')
scene<-groups<-NULL
svg<-cx$output('',paste0(stem,'.svg'))
cx$atomic(svg,function(p){
 svglite::svglite(p,width=168/25.4,height=.F$height/25.4,bg='white');grid.newpage();grid.rect(gp=gpar(fill='white',col=NA))
 groups<<-get(paste0('dj_fig',n))()
 # Explicit visible qualification in each withheld panel allocation.
 boxes<-list('2_B'=c(104,2,62,90),'4_B'=c(2,66,164,56),'5_A'=c(2,2,80,64),'5_B'=c(86,2,80,64),'8_A'=c(2,2,164,96))
 for(letter in setdiff(mapping$panel,included)){
  panel(letter,mapping$title[match(letter,mapping$panel)],boxes[[paste(n,letter,sep='_')]])
  dj_label('External input required',8,20);dj_label('Eligible-panel export',8,26)
 }
 upViewport(0);scene<<-grid.grab(warn=0);dev.off()
},function(p)stopifnot(file.info(p)$size>1000))
pdf<-cx$output('',paste0(stem,'.pdf'))
cx$atomic(pdf,function(p){stopifnot(system2(Sys.which('rsvg-convert'),c('--format=pdf','--output',shQuote(p),shQuote(svg)))==0)},function(p)stopifnot(file.info(p)$size>1000))
for(dpi in c(600,96,192)){
 dest<-cx$output('',paste0(stem,if(dpi==600)'' else paste0('_inspection_',dpi),'.png'))
 cx$atomic(dest,function(p){grDevices::png(p,width=168,height=.F$height,units='mm',res=dpi,type='cairo-png',bg='white',antialias='subpixel');grid.draw(scene);dev.off();stopifnot(system2('python3',c(shQuote(cx$input('workflow/scripts/set_png_resolution.py')),shQuote(p),as.character(dpi)))==0)},function(p)stopifnot(file.info(p)$size>1000))
}
if(journal=='JTIM')cx$atomic(cx$output('',paste0(stem,'.tiff')),function(p){grDevices::tiff(p,width=168,height=.F$height,units='mm',res=600,type='cairo',bg='white',compression='lzw',antialias='subpixel');grid.draw(scene);dev.off()},function(p)stopifnot(file.info(p)$size>1000))
filter_groups<-function(g)Filter(function(x)length(x)>1L,lapply(g,function(x)intersect(x,names(.F$regions))))
layout<-list(schema_version=1L,backend='R grid deviceLoc',figure=list(width_pt=168/25.4*72,height_pt=.F$height/25.4*72),panels=unname(.F$regions),row_groups=filter_groups(groups$rows),column_groups=filter_groups(groups$cols),exemptions=groups$ex)
for(x in list(c('data-layout.json','layout'),c('rendered-marks.json','.F$marks'),c('panel-layout.json','unname(.F$panels)')))
 cx$write_json(eval(parse(text=x[2])),cx$output('qa',x[1]))
cx$write_json(list(text=.F$text,axes=.F$axes,operations=.F$operations,legends=.F$legends,colorbars=.F$colorbars),cx$output('qa','presentation.json'))
cx$write_json(list(figure=n,journal=journal,requested_mode=render_mode,included_panels=included,external_panels=setdiff(mapping$panel,included),complete_figure=length(included)==nrow(mapping),new_inference=FALSE,runtime=runtime,font_paths=fonts,source_keys=unlist(lapply(.F$sources,function(d)d$binding_key))),cx$output('qa','render_metadata.json'))
cx$write_tsv(cx$verify_inputs(),cx$output('qa','verified_inputs.tsv'));cx$guard();cx$resource_check(6)
cat('Rendered',stem,'; complete figure:',length(included)==nrow(mapping),'\n')
