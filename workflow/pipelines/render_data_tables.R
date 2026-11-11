# Current eligible table-section export; consumes saved display rows, never fits.
args<-commandArgs(TRUE)
if(length(args)!=4L || !args[2]%in%c('JTIM','WJG'))stop('Usage: render_data_tables.R all|SECTION_ID JTIM|WJG INPUT_MANIFEST.tsv OUTPUT_DIRECTORY')
source('workflow/R/common/load_analysis.R');load_public_analysis()
cx<-public_analysis_context('current_table_sections',args[3],args[4],create=FALSE)
cx$resource_check();runtime<-cx$runtime('rendering')
suppressPackageStartupMessages(library(openxlsx2))
sections<-jsonlite::read_json(cx$input('publication/table_layout.json'))
directory<-cx$read_tsv('publication/TABLE_DIRECTORY.tsv')
if(args[1]!='all'){
 sections<-Filter(function(s)s$id==args[1],sections)
 if(length(sections)!=1L)stop('Unknown stable section ID; see publication/TABLE_DIRECTORY.tsv')
 directory<-directory[directory$section_id==args[1],,drop=FALSE]
}
stopifnot(length(sections)>0)
font<-if(args[2]=='JTIM')'Arial' else 'Book Antiqua'
if(!any(grepl(font,systemfonts::system_fonts()$family,fixed=TRUE)))stop('Required table font unavailable: ',font)
wb<-wb_workbook(creator='GC-VOC');wb<-wb_add_worksheet(wb,'Directory',grid_lines=FALSE)
wb<-wb_set_base_font(wb,font_size=if(args[2]=='JTIM')11 else 12,font_name=font,font_color=wb_color(hex='FF000000'))
wb<-wb_add_data(wb,'Directory','Eligible current table sections; held rows are omitted explicitly. This is not the complete manuscript workbook.',start_row=1,col_names=FALSE)
wb<-wb_merge_cells(wb,'Directory','A1:J2');wb<-wb_add_cell_style(wb,'Directory','A1:J2',wrap_text=TRUE,vertical='top')
directory_widths<-c(24,20,28,50,24,14,14,14,50,40)
wb<-wb_add_data(wb,'Directory',directory,start_row=4)
wb<-wb_set_col_widths(wb,'Directory',1:10,directory_widths)
wb<-wb_add_cell_style(wb,'Directory',wb_dims(rows=4:(4+nrow(directory)),cols=1:10),wrap_text=TRUE,vertical='top')
wb<-wb_add_font(wb,'Directory','A4:J4',bold=TRUE)
wb<-wb_set_row_heights(wb,'Directory',4,40)
for(i in seq_len(nrow(directory))) {
 ht<-max(40,18*max(ceiling(nchar(as.character(directory[i,]))/(directory_widths*.7)))+12)
 wb<-wb_set_row_heights(wb,'Directory',4+i,ht)
}
positions<-list();tabs<-unique(vapply(sections,function(s)s$table,''))
for(tab in tabs){
 wb<-wb_add_worksheet(wb,tab,grid_lines=FALSE,orientation='landscape',zoom=80)
 current<-1L;maxcols<-1L
 for(s in Filter(function(s)s$table==tab,sections)){
  visible<-as.integer(unlist(s$visible_columns));if(!length(visible))visible<-seq_along(s$headers)
  ncols<-length(visible);maxcols<-max(maxcols,ncols)
  wb<-wb_add_data(wb,tab,paste0(s$id,' | ',s$label,' | ',s$related_figures),start_row=current,col_names=FALSE)
  wb<-wb_merge_cells(wb,tab,wb_dims(rows=current,cols=seq_len(ncols)))
  wb<-wb_add_font(wb,tab,wb_dims(rows=current,cols=seq_len(ncols)),bold=TRUE)
  wb<-wb_set_row_heights(wb,tab,current,34)
  note<-paste0(if(nzchar(s$section_note))paste0(s$section_note,'\n') else '', 'Eligible rows: ',length(s$rows),' / ',s$source_rows,'; held rows: ',s$held_rows,'. Source/field mapping: publication/table_marks.tsv. Full precision: canonical scientific tables.')
  wb<-wb_add_data(wb,tab,note,start_row=current+1L,col_names=FALSE)
  wb<-wb_merge_cells(wb,tab,wb_dims(rows=current+1L,cols=seq_len(ncols)))
  wb<-wb_add_cell_style(wb,tab,wb_dims(rows=current:(current+1L),cols=seq_len(ncols)),wrap_text=TRUE,vertical='top')
  wb<-wb_set_row_heights(wb,tab,current+1L,max(40,18*ceiling(nchar(note)/(ncols*18))+12))
  if(length(s$rows)){
   dat<-as.data.frame(do.call(rbind,lapply(s$rows,unlist)),stringsAsFactors=FALSE)
   names(dat)<-unlist(s$headers);dat<-dat[,visible,drop=FALSE]
   wb<-wb_add_data(wb,tab,dat,start_row=current+2L,col_names=TRUE,with_filter=FALSE)
   dims<-wb_dims(rows=(current+2L):(current+2L+nrow(dat)),cols=seq_len(ncols))
   wb<-wb_add_cell_style(wb,tab,dims,horizontal='left',vertical='top',wrap_text=TRUE)
   wb<-wb_add_font(wb,tab,wb_dims(rows=current+2L,cols=seq_len(ncols)),bold=TRUE)
   wb<-wb_add_border(wb,tab,wb_dims(rows=current+2L,cols=seq_len(ncols)),top_border='thin',bottom_border='thin',left_border=NULL,right_border=NULL)
   wb<-wb_set_row_heights(wb,tab,current+2L,44)
   for(j in seq_len(nrow(dat))) {
    ht<-max(28,19*max(ceiling(nchar(as.character(dat[j,]))/18))+14)
    wb<-wb_set_row_heights(wb,tab,current+2L+j,ht)
   }
   positions[[s$id]]<-list(table=tab,first_data_row=current+3L,headers=names(dat),rows=unname(split(dat,seq_len(nrow(dat)))),related_figures=s$related_figures,source_rows=s$source_rows,held_rows=s$held_rows)
   current<-current+4L+nrow(dat)
  }else current<-current+4L
 }
 wb<-wb_set_col_widths(wb,tab,seq_len(maxcols),rep(27,maxcols))
 wb<-wb_freeze_pane(wb,tab,first_active_row=4,first_active_col=3)
 wb<-wb_set_page_setup(wb,tab,orientation='landscape',fit_to_width=1,fit_to_height=0,left=.35,right=.35,top=.4,bottom=.4)
 wb<-wb_add_named_region(wb,tab,dims=wb_dims(rows=1:max(1,current-1L),cols=seq_len(maxcols)),name='_xlnm.Print_Area',local_sheet=TRUE)
}
wb<-wb_set_active_sheet(wb,'Directory');wb<-wb_set_selected(wb,'Directory')
path<-cx$allowed(cx$output('',paste0('Eligible_sections_',args[2],'.xlsx')));stage<-paste0(path,'.stage.xlsx')
if(file.exists(stage))stop('Unresolved staged workbook')
wb_save(wb,stage,overwrite=FALSE);stopifnot(length(utils::unzip(stage,list=TRUE)$Name)>5)
stopifnot(file.rename(stage,path))
cx$write_json(list(status='ELIGIBLE_SECTIONS_ONLY',journal=args[2],section_ids=vapply(sections,function(s)s$id,''),held_rows=sum(vapply(sections,function(s)s$held_rows,0)),numeric_cells='Saved display tokens retained as text; full-precision numbers remain in canonical source tables',positions=positions,runtime=runtime),cx$output('qa','table_render.json'))
cx$write_tsv(cx$verify_inputs(),cx$output('qa','verified_inputs.tsv'));cx$guard();cx$resource_check(6)
cat('Exported eligible table sections; full workbook reproduction is not claimed.\n')
