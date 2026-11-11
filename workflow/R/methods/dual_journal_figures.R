# Current Figure 2-8 compositions; definitions only, no private source lookup.
source('workflow/R/methods/figure_rendering.R')
dj_panel_available <- function(letter) letter %in% .F$included
dj_sources <- function(n) {
 catalog<-cx$read_tsv('publication/figure_inputs.tsv',colClasses='character')
 catalog<-catalog[catalog$figure==as.character(n),]
 if(render_mode=='eligible-panels')catalog<-catalog[catalog$availability=='BUNDLED',]
 setNames(lapply(seq_len(nrow(catalog)),function(i){r<-catalog[i,];p<-cx$input(r$logical_input)
   d<-read.delim(p,quote='"',comment.char='',colClasses='character',na.strings=NULL,check.names=FALSE)
   fields<-strsplit(r$schema,';',fixed=TRUE)[[1]]
   if(r$availability=='BUNDLED') {
     if(!identical(names(d),fields))stop('Bundled figure schema mismatch: ',r$logical_input)
   } else {
     contract_path<-cx$input('workflow/config/figure_scientific_contracts.json')
     validator<-cx$input('workflow/pipelines/figure_scientific_inputs.py')
     message<-system2('python3',c(shQuote(validator),shQuote(contract_path),shQuote(r$logical_input),shQuote(p)),stdout=TRUE,stderr=TRUE)
     if(!is.null(attr(message,'status')) && attr(message,'status')!=0L)stop('Scientific figure input rejected: ',r$logical_input,'; ',paste(message,collapse=' '))
     contract<-jsonlite::read_json(contract_path)$inputs[[r$logical_input]]
     stopifnot(identical(fields,unlist(contract$fields,use.names=FALSE)))
     d<-d[,fields,drop=FALSE]
     d<-d[do.call(order,c(d[unlist(contract$keys)],list(method='radix'))),,drop=FALSE]
     if(r$role=='donors')d$display_index<-seq_len(nrow(d))
   }
   d$binding_key<-paste0(if(r$availability=='BUNDLED')r$canonical_path else r$logical_input,'#',seq_len(nrow(d)))
   d}),catalog$role)
}
dj_label <- function(text,x,y,face='plain',hjust=0,rot=0,col=ink) {
  if(!grepl('₂|g_av',text))return(label(text,x,y,.F$style$label,face,hjust,rot,col))
  lines<-strsplit(text,'\n',fixed=TRUE)[[1]];angle<-rot*pi/180
  for(li in seq_along(lines)) {
    chars<-strsplit(gsub('g_av','gₐᵥ',lines[li],fixed=TRUE),'',fixed=TRUE)[[1]]
    subscript<-chars%in%c('₂','ₐ','ᵥ');chars[chars=='₂']<-'2';chars[chars=='ₐ']<-'a';chars[chars=='ᵥ']<-'v'
    sizes<-ifelse(subscript,.F$style$label*.75,.F$style$label)
    widths<-vapply(seq_along(chars),function(i)convertWidth(grobWidth(textGrob(chars[i],gp=gpar(fontfamily=font_family(),fontsize=sizes[i],fontface=face))),'mm',TRUE),0)
    offsets<-c(0,head(cumsum(widths),-1))-hjust*sum(widths)
    yy<-y-(length(lines)-li)*.F$style$label*1.2/72*25.4
    for(i in seq_along(chars)) {
      drop<-if(subscript[i]).45 else 0
      label(chars[i],x+offsets[i]*cos(angle)+drop*sin(angle),yy-offsets[i]*sin(angle)+drop*cos(angle),sizes[i],face,0,rot,col)
    }
  }
}
dj_select <- function(d,...) {z<-d[...,,drop=FALSE];stopifnot(nrow(z)>0);z}
dj_one <- function(d,...) {z<-dj_select(d,...);stopifnot(nrow(z)==1);z}
dj_limits <- function(x,zero=TRUE,pad=.08) {
  x<-num(x);r<-range(c(x[is.finite(x)],if(zero)0));if(diff(r)==0)r<-r+c(-.5,.5)
  r+c(-1,1)*diff(r)*pad
}
dj_ticks <- function(lim,n=4) {a<-pretty(lim,n=n);a[a>=lim[1]&a<=lim[2]]}
dj_na <- function(d,x,y,reason='Unavailable') {
  dj_label('×',x,y+1,hjust=.5,col='#666666')
  mark_record('unavailable',d$binding_key,x,y,list(reason=reason))
}
# Layout uses measured label widths and actual rows. The larger of both
# journal metrics fixes shared wrapping without changing scientific geometry.
dj_text_width <- function(t) max(vapply(c('Arial','Tahoma'),function(f)
  convertWidth(grobWidth(textGrob(t,gp=gpar(fontfamily=f,fontsize=if(f=='Arial')8 else 7))),'mm',TRUE),0))
dj_legend <- function(labels,colors,center,y,types=rep('circle',length(labels)),nrow=1,width=NULL,broad=FALSE) {
  if(is.null(width))width<-.F$panels[[.F$current]]$box_mm[3]-4
  ww<-vapply(labels,dj_text_width,0)+4
  nr<-max(nrow,ceiling((sum(ww)+3*(length(labels)-1))/width))
  # Contiguous, balanced semantic order; 7 items requested in 2 rows gives 4+3.
  repeat {
    groups<-split(seq_along(labels),ceiling(seq_along(labels)/ceiling(length(labels)/nr)))
    if(all(vapply(groups,function(ix)sum(ww[ix])+3*(length(ix)-1)<=width,TRUE)))break
    nr<-nr+1;stopifnot(nr<=length(labels))
  }
  for(q in seq_along(groups)) {
    ix<-groups[[q]];gap<-if(broad&&length(ix)>1)(width-sum(ww[ix]))/(length(ix)-1) else 3
    st<-center-(sum(ww[ix])+gap*(length(ix)-1))/2
    for(k in seq_along(ix)) {
      i<-ix[k];yy<-y+(q-1)*4.8
      if(types[i]=='line')line(st,yy-1,st+2.2,yy-1,colors[i],.2,'dashed')
      else if(types[i]=='ci'){line(st,yy-1,st+2.2,yy-1,colors[i],.2);line(st,yy-1.7,st,yy-.3,colors[i],.2);line(st+2.2,yy-1.7,st+2.2,yy-.3,colors[i],.2)}
      else if(types[i]=='rect')rect(st+1.1,yy-1,2.2,1.7,colors[i])
      else if(types[i]=='hollow')rect(st+1,yy-1,1.1,1.1,'white',colors[i])
      else if(types[i]=='cross')dj_label('×',st+.9,yy,hjust=.5,col=colors[i])
      else dot(st+1,yy-1,.55,colors[i],shape=types[i])
      dj_label(labels[i],st+3,yy)
      .F$legends[[length(.F$legends)+1]]<-list(panel=.F$current,text=labels[i],x_mm=st,y_mm=yy,color=colors[i],type=types[i],band='below title, above data',center_mm=center,available_width_mm=width,row=q,label_width_mm=ww[i]-4)
      st<-st+ww[i]+gap
    }
  }
  invisible(y+(length(groups)-1)*4.8+3)
}
dj_hci <- function(d,r,y,standard=TRUE,col=blue,shape='circle',size=.65) {
  fields<-if(standard)c('g_av','g_av_ci_lower','g_av_ci_upper') else c('effect','ci_lower','ci_upper')
  v<-num(unlist(d[fields]));if(!is.finite(v[1])){xx<-mean(r$z[c(1,3)]);return(dj_na(d,xx,y))}
  stopifnot(v[1]>=r$xlim[1],v[1]<=r$xlim[2])
  if(all(is.finite(v[2:3]))){stopifnot(min(v)>=r$xlim[1],max(v)<=r$xlim[2]);line(r$x(v[2]),y,r$x(v[3]),y,col,.2);for(a in v[2:3])line(r$x(a),y-.65,r$x(a),y+.65,col,.2)}
  dot(r$x(v[1]),y,size,col,shape=shape)
  mark_record('estimate_CI',d$binding_key,r$x(v[1]),y,list(region=r$id,fields=fields,values=v,CI_available=all(is.finite(v[2:3]))))
}
dj_vci <- function(d,r,x,standard=FALSE,col=blue,bar=FALSE,w=2) {
  fields<-if(standard)c('g_av','g_av_ci_lower','g_av_ci_upper') else c('effect','ci_lower','ci_upper')
  v<-num(unlist(d[fields]));if(!is.finite(v[1])){return(dj_na(d,r$x(x),r$y(0)))}
  stopifnot(v[1]>=r$ylim[1],v[1]<=r$ylim[2])
  if(bar)cell_mark(d,r$x(x),mean(r$y(c(0,v[1]))),w,abs(diff(r$y(c(0,v[1])))),v[1],fields[1],r$ylim,'bar',col,list(region=r$id))
  if(all(is.finite(v[2:3]))){stopifnot(min(v)>=r$ylim[1],max(v)<=r$ylim[2]);line(r$x(x),r$y(v[2]),r$x(x),r$y(v[3]),ink,.18);for(a in v[2:3])line(r$x(x)-.8,r$y(a),r$x(x)+.8,r$y(a),ink,.18)}
  if(!bar)dot(r$x(x),r$y(v[1]),.6,col)
  mark_record('estimate_CI_vertical',d$binding_key,r$x(x),r$y(v[1]),list(region=r$id,fields=fields,values=v,CI_available=all(is.finite(v[2:3]))))
}
# Interrupt reference guides around explicit unavailable slots (no data removed).
dj_zero_gaps <- function(r,axis,centers,half=2) {
  lim<-if(axis=='x')r$z[c(2,4)] else r$z[c(1,3)]
  edges<-sort(c(lim,as.vector(rbind(centers-half,centers+half))))
  for(i in seq(1,length(edges)-1,by=2)) {
    if(axis=='x')line(r$x(0),edges[i],r$x(0),edges[i+1],connector,.12,'dashed')
    else line(edges[i],r$y(0),edges[i+1],r$y(0),connector,.12,'dashed')
  }
  .F$operations[[paste0(r$id,'_guide_gaps')]]<-list(operation='Reference guide interrupted around unavailable glyph; no numerical geometry altered',axis=axis,missing_centers_mm=centers,half_gap_mm=half)
}
dj_colorbar <- function(id,x,top,bottom,lim,title,title_x=x+3) {
  grad<-grid::linearGradient(c(blue,'#F1F3F5',orange),stops=c(0,.5,1),x1=0,y1=0,x2=0,y2=1)
  rect(x,mean(c(top,bottom)),2.2,bottom-top,grad)
  vals<-c(lim[1],0,lim[2]);for(v in vals){yy<-bottom-(v-lim[1])/diff(lim)*(bottom-top);line(x+1.1,yy,x+1.8,yy,ink,.1);dj_label(tick_labels(v),x+2.3,yy+.8)}
  dj_label(title,title_x,top-3,hjust=.5)
  .F$colorbars[[id]]<-list(panel=.F$current,limits=lim,title=title,top_mm=top,bottom_mm=bottom,x_mm=x,ticks=vals,direction='low bottom, high top',tick_side='right')
}
dj_heat <- function(d,x,y,w,h,field,lim) {
  v<-num(d[[field]]);cell_mark(d,x,y,w,h,v,field,lim)
  if(!is.finite(v))dj_na(d,x,y)
}
dj_pair <- function(d,r,nat='NAT_score',tum='Tumor_score',id='patient_id',factor=1) {
  d<-d[order(d[[id]],method='radix'),];a<-num(d[[nat]])*factor;b<-num(d[[tum]])*factor
  offsets<-(seq_len(nrow(d))-(nrow(d)+1)/2)/max(1,nrow(d)-1)*4
  for(i in seq_len(nrow(d))){xx<-r$x(c(.2,.8))+offsets[i];line(xx[1],r$y(a[i]),xx[2],r$y(b[i]),connector,.09)
    point_xy(d[i,],r,.2,a[i],yfield=nat,dx=offsets[i],radius=.38,fill=blue);point_xy(d[i,],r,.8,b[i],yfield=tum,dx=offsets[i],radius=.38,fill=orange)
    mark_record('pair_connector',d$binding_key[i],xx[1],r$y(a[i]),list(region=r$id,identity=d[[id]][i],fields=c(nat,tum),values=c(a[i],b[i]),factor=factor))}
  dj_label('Adjacent',r$x(.2),r$z[4]+4,hjust=.5);dj_label('Tumor',r$x(.8),r$z[4]+4,hjust=.5)
}
dj_boxplot <- function(d,r,x,field,col=blue,id='patient_id') {
  v<-num(d[[field]]);stopifnot(all(is.finite(v)),min(v)>=r$ylim[1],max(v)<=r$ylim[2]);x0<-r$x(x)
  stopifnot(!anyDuplicated(d[[id]]))
  q<-as.numeric(quantile(v,c(.25,.5,.75),type=7));iq<-q[3]-q[1]
  wh<-range(v[v>=q[1]-1.5*iq & v<=q[3]+1.5*iq]);w<-3.6
  off<-(((seq_len(nrow(d))-1)%%11)/10-.5)*4.6
  for(i in seq_len(nrow(d)))point_xy(d[i,],r,x,v[i],yfield=field,dx=off[i],radius=.24,fill=adjustcolor(col,alpha.f=.48))
  line(x0,r$y(wh[1]),x0,r$y(q[1]),col,.18);line(x0,r$y(q[3]),x0,r$y(wh[2]),col,.18)
  for(vv in wh)line(x0-w/3,r$y(vv),x0+w/3,r$y(vv),col,.18)
  rect(x0,mean(r$y(q[c(1,3)])),w,abs(diff(r$y(q[c(1,3)]))),adjustcolor('white',.7),col)
  line(x0-w/2,r$y(q[2]),x0+w/2,r$y(q[2]),col,.3)
  mark_record('boxplot',paste(d$binding_key,collapse=';'),x0,r$y(q[2]),list(region=r$id,source_keys=d$binding_key,source_field=field,quartiles=q,whiskers=wh,width_mm=w,n=length(v),outlier_symbols=FALSE))
  .F$operations[[paste0(r$id,'_',x)]]<-list(operation='Boxplot: type-7 quartiles; whiskers are observed extrema within 1.5 IQR; all patient points drawn once with deterministic 11-column offsets; no extra outlier symbols; whiskers are not confidence intervals',source_keys=d$binding_key,field=field,n=length(v),quartiles=q,whiskers=wh)
}

dj_fig2 <- function() {
  if(dj_panel_available('A')) {
  panel('A','RNA program changes',c(2,2,98,90));p<-S('rna')
  breaks<-lapply(c('G','H'),function(g){vv<-c(0,num(p$delta[p$program==g]));raw<-diff(range(vv))/12;base<-10^floor(log10(raw));step<-ceiling(raw/base)*base;seq(floor(min(vv)/step),ceiling(max(vv)/step))*step})
  histograms<-list();for(g in 1:2)for(i in 1:4){z<-p[p$source==c('GSE29272','GSE27342_64','TCGA_STAD','GSE122401')[i]&p$program==c('G','H')[g],];if(nrow(z))histograms[[paste(g,i)]]<-hist(num(z$delta),breaks=breaks[[g]],right=FALSE,include.lowest=TRUE,plot=FALSE)$counts/nrow(z)*100}
  ymax<-ceiling(max(unlist(histograms))/10)*10
  dj_label('Patients (%)',1,11)
  for(g in 1:2){left<-32+(g-1)*35;dj_label(c('G','H')[g],left+14,11,'bold',.5)
    for(i in 1:4){z<-p[p$source==c('GSE29272','GSE27342_64','TCGA_STAD','GSE122401')[i]&p$program==c('G','H')[g],];top<-16+(i-1)*16
      r<-region(paste0('A',g,'_',i),c(left,top,left+28,top+11),range(breaks[[g]]),c(0,ymax))
      if(g==1){dj_label(c('GSE29272','GSE27342','TCGA-STAD','GSE122401')[i],1,top+4);dj_label(c('n = 112','n = 49','H: n = 33','n = 80')[i],1,top+8)}
      if(!nrow(z)){dj_label('Unavailable',left+14,top+7,hjust=.5);ref<-S('primary');ref<-dj_one(ref,ref$display_source==3&ref$display_program=='G');mark_record('unavailable',ref$binding_key,left+14,top+5.5,list(reason='Full TCGA G unavailable'));next}
      vv<-num(z$delta);bins<-findInterval(vv,breaks[[g]],rightmost.closed=TRUE);stopifnot(all(bins>=1&bins<length(breaks[[g]])))
      counts<-tabulate(bins,nbins=length(breaks[[g]])-1);percent<-counts/nrow(z)*100;stopifnot(sum(counts)==nrow(z),abs(sum(percent)-100)<1e-12)
      for(k in seq_along(counts)){xx<-mean(r$x(breaks[[g]][k+c(0,1)]));ww<-diff(r$x(breaks[[g]][k+c(0,1)]));hh<-r$y(0)-r$y(percent[k]);rect(xx,r$y(percent[k]/2),ww,hh,adjustcolor(c(blue,orange)[g],.72),'white')
        mark_record('histogram_bin',z$binding_key[1],xx,r$y(percent[k]/2),list(region=r$id,bin=k,breaks=breaks[[g]][k+c(0,1)],count=counts[k],percentage=percent[k],n=nrow(z),width_mm=ww,height_mm=hh,source_keys=z$binding_key[bins==k]))}
      for(k in seq_len(nrow(z)))mark_record('histogram_observation',z$binding_key[k],r$x(mean(breaks[[g]][bins[k]+c(0,1)])),r$y(percent[bins[k]]),list(region=r$id,bin=bins[k],field='delta',value=vv[k]))
      zero_x(r);yaxis(r,c(0,ymax/2,ymax),show_labels=g==1)
      if(i==4)xaxis(r,dj_ticks(r$xlim,3)) else line(left,top+11,left+28,top+11,ink,.1)
      .F$operations[[r$id]]<-list(operation='Equal-width percentage histogram; target 12 bins. Step = ceiling((range including zero / 12) / 10^floor(log10(range/12))) * 10^floor(log10(range/12)); outer limits rounded outwards to that grid. Left closed/right open, final bin includes upper edge. Same program breaks across cohorts; common percentage axis.',source_keys=z$binding_key,field='delta',breaks=breaks[[g]],counts=counts,percentages=percent,n=nrow(z),common_percentage_limit=ymax)
    }
  };shared_title('Tumor − adjacent rank score',c('A1_4','A2_4'),86,size=.F$style$label)
  }
  if(dj_panel_available('B')) {
  panel('B','Protein program changes',c(104,2,62,90));dj_legend(c('Adjacent','Tumor'),c(blue,orange),36,11,width=52)
  for(g in 1:2){top<-21+(g-1)*36;r<-region(paste0('B',g),c(16,top,59,top+25),c(0,1),c(0,1));d<-S(c('protein_G','protein_H')[g]);dj_label(paste0('CPTAC ',c('G6','H4')[g],' · n = ',nrow(d)),37.5,top-3,'bold',.5);yaxis(r,c(0,.5,1),'Rank score',title_x=3);dj_pair(d,r)}
  }
  if(dj_panel_available('C')) {
  panel('C','Across-cohort effects',c(2,98,98,65));d<-S('primary')
  for(g in 1:2){left<-29+(g-1)*35;r<-region(paste0('C',g),c(left,17,left+30,51),c(-3.6,1.7));if(g==1)dj_zero_gaps(r,'x',32.6) else zero_x(r);dj_label(c('G','H')[g],left+15,12,'bold',.5)
    for(i in 1:5){z<-dj_one(d,d$display_program==c('G','H')[g]&d$display_source==i);y<-20+c(0,7,14,21,31)[i]*.9
      if(g==1)dj_label(c('GSE29272','GSE27342','TCGA-STAD','GSE122401','CPTAC protein')[i],1,y+1);dj_hci(z,r,y,col=c(blue,orange)[g])}
    xaxis(r,c(-3,0,1))
  };shared_title('Hedges’ g_av (95% CI)',c('C1','C2'),62,size=.F$style$label)
  }
  if(dj_panel_available('D')) {
  panel('D','Anatomical strata',c(104,98,62,65));dj_legend(c('G','H'),c(blue,orange),37.5,11,types=c('rect','rect'));d<-S('strata');r<-region('D',c(16,17,59,53),c(.5,2.5),c(-3.6,1.7));zero_y(r);yaxis(r,c(-3,-2,-1,0,1),'Hedges’ g_av',title_x=3)
  for(i in 1:2){for(g in 1:2){z<-dj_one(d,d$construct==c('G_core','H_core')[g]&grepl(c('cardia','non_cardia')[i],d$specification)&!if(i==1)grepl('non',d$specification) else FALSE);dj_vci(z,r,i+c(-.16,.16)[g],TRUE,c(blue,orange)[g],TRUE,5)}
    dj_label(c('Cardia','Non-cardia')[i],r$x(i),58,hjust=.5);dj_label(paste0('n = ',c(48,64)[i]),r$x(i),62,hjust=.5)}
  }
  list(rows=c(lapply(1:4,function(i)paste0('A',1:2,'_',i)),list(c('C1','C2'))),cols=list(paste0('A1_',1:4),paste0('A2_',1:4),c('B1','B2')),ex=list())
}

dj_fig3 <- function() {
  if(dj_panel_available('A')) {
  panel('A','G member patterns',c(2,2,164,53));dj_legend('Unavailable','#777777',80,10,'cross');d<-S('rna_members');p<-S('protein_members')
  rl<-c(-3,3);pl<-c(-1,1);region('A_matrix',c(28,22,128,51));dj_label('RNA',68,15,'bold',.5);dj_label('Protein',118,15,'bold',.5)
  for(j in 1:5){x<-28+(j-.5)*20;dj_label(c('GSE29272','GSE27342','TCGA-STAD','GSE122401','CPTAC')[j],x,20,hjust=.5)
    for(k in 1:7){y<-24+(k-1)*4.2;if(j==1)dj_label(G[k],25,y+1,'italic',1)
      z<-if(j<5)dj_one(d,d$display_column==j&d$display_member==G[k]) else p[p$member==G[k],,drop=FALSE]
      if(!nrow(z)){rect(x,y,19.5,4,unavailable);dj_label('×',x,y+1,hjust=.5);ref<-S('primary');ref<-dj_one(ref,ref$display_source==5&ref$display_program=='G');stopifnot(!grepl('ALOX12',ref$members));mark_record('unavailable',ref$binding_key,x,y,list(reason='ALOX12 absent from accepted CPTAC G6 members',field='members',value=ref$members))} else dj_heat(z,x,y,19.5,4,'effect',if(j<5)rl else pl)
    }}
  dj_label('Assembly',1,31,rot=90);dj_label('Oxidation',1,48,rot=90)
  dj_colorbar('A_RNA',135,29,50,rl,'RNA\nΔlog₂');dj_colorbar('A_protein',152,29,50,pl,'Protein\nΔlog₂')
  }
  if(dj_panel_available('B')) {
  panel('B','Functional submodules',c(2,60,82,65));dj_legend(c('Assembly','Oxidation'),c(blue,orange),46,11,types=c('rect','rect'));d<-S('submodules');r<-region('B',c(14,21,79,53),c(.5,4.5),c(-.8,1.8));dj_zero_gaps(r,'y',r$x(c(3.17,4.17)));yaxis(r,c(-.5,0,1,1.5),'g_av',title_x=2)
  for(i in 1:4)for(j in 1:2){z<-dj_one(d,d$display_source==c(1,2,3,5)[i]&d$display_submodule==c('substrate_assembly','oxidation_support')[j]);dj_vci(z,r,i+c(-.17,.17)[j],TRUE,c(blue,orange)[j],TRUE,4)}
  for(i in 1:4)dj_label(c('GSE29272','GSE27342','TCGA-STAD','CPTAC')[i],r$x(i),59,hjust=.5)
  }
  if(dj_panel_available('C')) {
  panel('C','Scoring sensitivity',c(88,60,78,65));dj_legend(c('GSE29272','GSE27342','GSE122401','CPTAC'),c(blue,blue,blue,orange),44,11,c('circle','square','diamond','circle'),2,width=66)
  r<-region('C',c(25,23,59,57),c(-.2,1.2),c(-.2,1.2));line(r$x(-.2),r$y(-.2),r$x(1.2),r$y(1.2),'#999999',.12,'dashed')
  a<-S('primary');b<-S('alternative');for(i in c(1,2,4,5)){x<-dj_one(a,a$display_source==i&a$display_program=='G');y<-dj_one(b,b$display_source==i&b$display_program=='G');dot(r$x(num(x$g_av)),r$y(num(y$g_av)),.8,if(i==5)orange else blue,shape=c('circle','square','diamond','circle')[match(i,c(1,2,4,5))]);mark_record('concordance',paste(x$binding_key,y$binding_key,sep=';'),r$x(num(x$g_av)),r$y(num(y$g_av)),list(region='C',x_value=num(x$g_av),y_value=num(y$g_av),field='g_av'))}
  xaxis(r,c(0,.5,1),'Rank g_av',68);yaxis(r,c(0,.5,1),'Mean-z g_av',title_x=12)
  }
  if(dj_panel_available('D')) {
  panel('D','Member influence',c(2,134,108,97));dj_legend(c('Omission','Intact','95% CI'),c(blue,ink,ink),58,11,c('circle','line','ci'))
  a<-S('omissions');b<-S('primary');for(j in 1:4){src<-c(1,2,4,5)[j];z<-a[a$display_source==src&a$display_program=='G',];genes<-if(src==5)setdiff(G,'ALOX12') else G;z<-z[match(genes,z$display_omitted),];stopifnot(nrow(z)==length(genes),!anyNA(z$binding_key));ref<-dj_one(b,b$display_source==src&b$display_program=='G')
    left<-19+((j-1)%%2)*53;top<-23+((j-1)%/%2)*39;lim<-dj_limits(c(z$ci_lower,z$ci_upper,ref$effect));r<-region(paste0('D',j),c(left,top,left+31,top+24),lim);zero_x(r);line(r$x(num(ref$effect)),top,r$x(num(ref$effect)),top+24,ink,.13,'dashed')
    dj_label(c('GSE29272','GSE27342','GSE122401','CPTAC (G6)')[j],left+15.5,top-3,'bold',.5)
    for(i in seq_along(genes)){y<-top+1.5+(i-1)*3.5;dj_label(genes[i],left-2,y+1,'italic',1);dj_hci(z[i,],r,y,FALSE)}
    xaxis(r,dj_ticks(lim,2))
    mark_record('intact_reference',ref$binding_key,r$x(num(ref$effect)),top+12,list(region=r$id,field='effect',value=num(ref$effect),orientation='vertical'))
  };shared_title('Δ rank score',c('D3','D4'),97,size=.F$style$label);dj_label('TCGA-STAD: unavailable',54,16,hjust=.5)
  }
  if(dj_panel_available('E')) {
  panel('E','Buffering context',c(114,134,52,97));a<-S('context_rank');b<-S('context_z');dj_legend(c('Rank','Mean-z'),c(blue,orange),28,11,c('circle','square'),width=49)
  r<-region('E',c(24,23,49,86),c(-1.5,1.5));dj_zero_gaps(r,'x',c(27,71),half=4.4)
  for(i in 1:5){y<-27+(i-1)*14;dj_label(c('GSE29272','GSE27342','TCGA-STAD','GSE122401','CPTAC')[i],1,y+1)
    for(j in 1:2){d<-if(j==1)a else b;z<-d[d$display_source==i&d$display_program=='G_buffer',];yy<-y+c(-2,2)[j]
      if(nrow(z))dj_hci(z,r,yy,TRUE,c(blue,orange)[j],c('circle','square')[j]) else {dj_label('×',mean(r$z[c(1,3)]),yy+1,hjust=.5);ref<-dj_one(a,a$display_source==i&a$display_program=='G_buffer');mark_record('unavailable_method',ref$binding_key,mean(r$z[c(1,3)]),yy,list(field='display_program',value='G_buffer',reason='Mean-z not computed in original external context package; complete alternative table retained'))}}}
  xaxis(r,c(-1,0,1),'Hedges’ g_av',97)
  }
  list(rows=list(c('D1','D2'),c('D3','D4')),cols=list(c('D1','D3'),c('D2','D4')),ex=list())
}

dj_fig4 <- function() {
  if(dj_panel_available('A')) {
  panel('A','H member expression',c(2,2,164,59));d<-S('rna_members');d<-d[d$display_member%in%H,]
  stopifnot(nrow(d)==16,all(is.finite(num(d$effect))),all(is.finite(num(d$ci_lower))),all(is.finite(num(d$ci_upper))))
  cols<-c('#6C7680','#009E73','#D55E00','#0072B2')
  dj_legend(c('GSE29272','GSE27342','TCGA-STAD','GSE122401'),cols,88,11,types=rep('rect',4),width=143)
  lim<-dj_limits(c(d$effect,d$ci_lower,d$ci_upper));r<-region('A',c(20,18,160,47),c(.5,4.5),lim);zero_y(r)
  yaxis(r,dj_ticks(lim,4),'Native RNA coefficient (log₂)',title_x=3)
  for(g in 1:4){for(i in 1:4){z<-dj_one(d,d$display_column==i&d$display_member==H[g]);dj_vci(z,r,g+c(-.3,-.1,.1,.3)[i],FALSE,cols[i],TRUE,5)};dj_label(H[g],r$x(g),54,'italic',.5)}
  }
  if(dj_panel_available('B')) {
  panel('B','H protein changes',c(2,66,164,56));for(j in 1:2){d<-S(c('protein_changes','mun_ratios')[j]);if(j==1)d<-d[d$population_reference_id=='protein_H_program'&d$member%in%H,];field<-c('delta','centered_log2_T_over_adjacent_ratio')[j];left<-17+(j-1)*84;lim<-dj_limits(d[[field]]);r<-region(paste0('B',j),c(left,18,left+60,46),c(.5,4.5),lim);zero_y(r)
    dj_label(paste0(c('CPTAC','Mun')[j],' · n = ',c(30,80)[j]),left+30,12,'bold',.5);yaxis(r,dj_ticks(lim,2),if(j==1)'Paired Δlog₂ ratio' else 'Centered log₂ T/N',title_x=left-14)
    for(i in 1:4){z<-d[d$member==H[i],];stopifnot(nrow(z)==c(30,80)[j]);dj_boxplot(z,r,i,field,c(blue,orange)[j]);dj_label(H[i],r$x(i),53,hjust=.5)}
  }
  }
  if(dj_panel_available('C')) {
  panel('C','H robustness',c(2,128,73,64));dj_legend(c('Rank','Mean-z'),c(blue,orange),48,11,c('circle','square'));a<-S('primary');b<-S('alternative');r<-region('C',c(26,21,70,52),c(-3.7,.2));zero_x(r)
  for(i in 1:5){y<-25+(i-1)*6;dj_label(c('GSE29272','GSE27342','TCGA-STAD','GSE122401','CPTAC')[i],1,y+1);for(j in 1:2){d<-if(j==1)a else b;z<-dj_one(d,d$display_source==i&d$display_program=='H');dj_hci(z,r,y+c(-1.2,1.2)[j],TRUE,c(blue,orange)[j],c('circle','square')[j],.55)}}
  xaxis(r,c(-3,-2,-1,0),'Hedges’ g_av (95% CI)',63)
  }
  if(dj_panel_available('D')) {
  panel('D','Related programs',c(79,128,87,64));dj_legend('Unavailable','#777777',40,11,'cross');d<-S('context');region('D_matrix',c(25,22,76,52))
  for(j in 1:3)dj_label(c('Expanded','Reductive','Conjugation')[j],25+(j-.5)*17,18,hjust=.5)
  for(i in 1:5){y<-25+(i-1)*6;dj_label(c('GSE29272','GSE27342','TCGA-STAD','GSE122401','CPTAC')[i],23,y+1,hjust=1);for(j in 1:3){z<-dj_one(d,d$display_source==i&d$display_program==c('H_expanded','H_reductive','H_conjugation')[j]);dj_heat(z,25+(j-.5)*17,y,16.7,5.7,'g_av',c(-3,3))}}
  dj_colorbar('D',81,24,51,c(-3,3),'g_av')
  }
  list(rows=list(c('B1','B2')),cols=list(),ex=list())
}

dj_fig5 <- function() {
  for(j in 1:6){l<-LETTERS[j];g<-if(j<=2)c('G','H')[j] else H[j-2];if(!dj_panel_available(l))next;panel(l,if(j<=2)paste(g,'program association') else g,c(2+((j-1)%%2)*84,c(2,67,129)[1+(j-1)%/%2],80,if(j<=2)64 else 61))
    if(j<=2){d<-S(paste0(g,'_pairs'));e<-S(paste0(g,'_result'));xf<-'RNA_delta';yf<-'protein_delta';rho<-e$rho;src<-'CPTAC';xl<-'RNA Δ rank score';yl<-'Protein Δ rank score'}
    else {d<-S('E_pairs');d<-d[d$member==g,];e<-S('E_results');e<-e[e$member==g,];xf<-'RNA_delta';yf<-'protein_mean_centered_log2_T_N';rho<-e$effect;src<-'Mun';xl<-'RNA Δlog₂(UQ + 1)';yl<-'Protein centered log₂ T/N'}
    stopifnot(nrow(d)==c(19,21,80,80,80,80)[j],!anyDuplicated(d$patient_id))
    if(j<=2){dj_label(paste0(src,' · n = ',nrow(d),' · ρ = ',fmt(rho)),43,11,hjust=.5);dj_label(paste0('95% CI ',fmt(e$ci_lower),' to ',fmt(e$ci_upper)),43,16,hjust=.5)}
    else dj_label(paste0('ρ = ',fmt(rho),' · 95% CI ',fmt(e$ci_lower),' to ',fmt(e$ci_upper)),43,11,hjust=.5)
    mark_record('association_summary',e$binding_key,43,12,list(fields=c(if(j<=2)'rho' else 'effect','ci_lower','ci_upper'),values=c(num(rho),num(e$ci_lower),num(e$ci_upper))))
    shift<-if(j>=3)3 else 0;r<-region(l,c(18,20-shift,77,53-shift),dj_limits(d[[xf]]),dj_limits(d[[yf]]));zero_x(r);zero_y(r)
    for(i in seq_len(nrow(d)))point_xy(d[i,],r,num(d[[xf]][i]),num(d[[yf]][i]),xf,yf,radius=.5,fill=if(j<=2)blue else orange)
    xaxis(r,dj_ticks(r$xlim,4),xl,64-shift);yaxis(r,dj_ticks(r$ylim,3),yl,title_x=3)
  }
  list(rows=list(c('A','B'),c('C','D'),c('E','F')),cols=list(c('A','C','E'),c('B','D','F')),ex=list())
}

dj_fig6 <- function() {
  panel('A','Gastric marker changes',c(2,2,164,48));dj_legend('Unavailable','#777777',84,10,'cross');d<-S('markers');genes<-c('GKN1','GKN2','MUC5AC','TFF1','MUC6','TFF2','PGC');region('A_matrix',c(30,20,138,45))
  for(j in 1:4){dj_label(c('GSE29272','GSE27342','GSE122401','TCGA-STAD')[j],30+(j-.5)*27,18,hjust=.5);for(k in 1:7){y<-22+(k-1)*3.5;z<-dj_one(d,d$source==c('GSE29272','GSE27342_core','GSE122401','TCGA_STAD')[j]&d$member==genes[k]);if(j==1)dj_label(genes[k],27,y+1,'italic',1);dj_heat(z,30+(j-.5)*27,y,26.7,3.3,'standardized_mean_paired_change',c(-1.8,1.8))}}
  dj_label('Pit',3,29,rot=90);dj_label('Gland',3,40,rot=90);dj_label('PGC',3,45);dj_colorbar('A',146,25,44,c(-1.8,1.8),'Mean Δz')
  for(j in 1:2){l<-c('B','C')[j];panel(l,paste0('H–D7: ',c('GSE122401','TCGA-STAD')[j]),c(2+(j-1)*84,54,80,66));d<-S(c('GSE122401_pairs','TCGA_pairs')[j]);f<-S(c('GSE122401_fit','TCGA_fit')[j]);col<-c(blue,orange)[j];dj_label(paste0('n = ',nrow(d)),45,11,hjust=.5)
    r<-region(l,c(18,17,77,50),c(-3.5,2),c(-.45,.22));zero_x(r);zero_y(r)
    grid.polygon(x=mm(r$x(c(num(f$delta_D),rev(num(f$delta_D))))),y=unit(1,'npc')-mm(r$y(c(num(f$ci_lower),rev(num(f$ci_upper))))),gp=gpar(fill=adjustcolor(col,.14),col=NA))
    grid.lines(x=mm(r$x(num(f$delta_D))),y=unit(1,'npc')-mm(r$y(num(f$fit))),gp=gpar(col=col,lwd=.8))
    .F$operations[[paste0(l,'_saved_fit')]]<-list(source_keys=f$binding_key,fields=c('delta_D','fit','ci_lower','ci_upper'),new_fit=FALSE)
    for(i in seq_len(nrow(d)))point_xy(d[i,],r,num(d$delta_D[i]),num(d$delta_H[i]),'delta_D','delta_H',radius=.5,fill=col)
    xaxis(r,-3:2,'Δ D7 score',61);yaxis(r,c(-.4,-.2,0,.2),'Δ H rank score',title_x=3)
  }
  panel('D','Component associations',c(2,126,80,98));dj_legend(c('GSE122401','TCGA-STAD'),c(blue,orange),47.5,11,types=c('rect','rect'),width=64)
  rect(29,16,3,2,adjustcolor(blue,.18),blue);dj_label('95% CI',32,17);line(49,14.8,49,17.2,ink,.3);dj_label('Estimate',52,17)
  d<-S('D_results');r<-region('D',c(18,22,77,87),c(0,.1));zero_x(r)
  for(i in 1:3){y<-29+(i-1)*25;dj_label(c('Pit','Gland','PGC')[i],2,y+1);for(j in 1:2){z<-dj_one(d,d$source==c('GSE122401','TCGA_STAD')[j]&d$component==c('D_pit4','D_gland2','PGC_z')[i]);v<-num(unlist(z[c('effect','ci_lower','ci_upper')]));yy<-y+c(-3,3)[j];thickness<-3
    stopifnot(all(is.finite(v)),min(v)>=0,max(v)<=.1,v[2]<=v[1],v[1]<=v[3])
    rect(mean(r$x(v[2:3])),yy,diff(r$x(v[2:3])),thickness,adjustcolor(c(blue,orange)[j],.18),c(blue,orange)[j])
    line(r$x(v[1]),yy-thickness/2,r$x(v[1]),yy+thickness/2,c(blue,orange)[j],.35)
    mark_record('horizontal_CI_crossbar',z$binding_key,r$x(v[1]),yy,list(region='D',fields=c('effect','ci_lower','ci_upper'),values=v,thickness_mm=thickness,fill=adjustcolor(c(blue,orange)[j],.18),outline=c(blue,orange)[j],estimate_stroke=TRUE))
  }}
  xaxis(r,seq(0,.1,.02),digits=2);shared_title('Rank-H4 change / component unit','D',98,size=.F$style$label)
  panel('E','Component correlations',c(86,126,80,49));d<-S('correlations')
  for(j in 1:2){left<-9+(j-1)*31;top<-17;r<-region(paste0('E',j),c(left,top,left+26,top+26));dj_label(c('GSE122401','TCGA-STAD')[j],left+13,12,'bold',.5)
    for(i in 1:3){if(j==1)dj_label(c('Pit','Gland','PGC')[i],left-2,top+(i-.5)*26/3+1,hjust=1);dj_label(c('Pit','Gland','PGC')[i],left+(i-.5)*26/3,49,hjust=.5)
      for(k in 1:3){x<-left+(k-.5)*26/3;y<-top+(i-.5)*26/3;if(i==k){rect(x,y,8.45,8.45,'#E1E1E1');dj_label('—',x,y+1,hjust=.5)}else {cc<-c('D_pit4','D_gland2','PGC_z');z<-d[d$source==c('GSE122401','TCGA_STAD')[j]&((d$component_1==cc[i]&d$component_2==cc[k])|(d$component_1==cc[k]&d$component_2==cc[i])),];stopifnot(nrow(z)==1);dj_heat(z,x,y,8.45,8.45,'rho',c(-1,1))}}}}
  dj_colorbar('E',72,22,43,c(-1,1),'Spearman\nρ',title_x=71)
  panel('F','Scoring sensitivity',c(86,181,80,43));dj_legend(c('GSE122401','TCGA-STAD'),c(blue,orange),44,11,types=c('rect','rect'),width=73);d<-S('coefficients')
  for(j in 1:2){z<-d[d$model=='M1'&d$coefficient=='delta_D'&d$H_definition==c('accepted_rank_H','standardized_member_mean_H')[j],];lim<-dj_limits(c(z$ci_lower,z$ci_upper));left<-13+(j-1)*39;r<-region(paste0('F',j),c(left,22,left+23,35),c(.5,2.5),lim);dj_label(c('Rank H','Mean-z H')[j],left+11.5,18,'bold',.5);zero_y(r);yaxis(r,dj_ticks(lim,2));for(i in 1:2)dj_vci(dj_one(z,z$source==c('GSE122401','TCGA_STAD')[i]),r,i,FALSE,c(blue,orange)[i],TRUE,6)}
  dj_label('Slope per D7 score unit',44,43,hjust=.5)
  list(rows=list(c('B','C'),c('E1','E2'),c('F1','F2')),cols=list(),ex=list())
}

dj_fig7 <- function() {
  if(dj_panel_available('A')) {
  panel('A','Compartment programs',c(2,2,68,73));d<-S('programs');d<-d[d$display_method=='Rank score',];dj_legend(c('Formal','Exploratory'),c(ink,ink),47.5,11,c('circle','diamond'));r<-region('A',c(30,18,65,61),c(-.12,.12));zero_x(r)
  for(i in 1:5){z<-dj_one(d,d$display_order==i);y<-22+(i-1)*8;dj_label(c('Myeloid · G','Myeloid · H','Fibroblast · H','B/plasma · H','Epithelial · H')[i],1,y-1);dj_label(paste0('n = ',z$biological_n),1,y+3);dj_hci(z,r,y,FALSE,if(i==1)blue else orange,if(i==5)'diamond' else 'circle')}
  xaxis(r,c(-.1,0,.1),'Δ rank score (95% CI)',72)
  }
  if(dj_panel_available('B')) {
  panel('B','Members and coverage',c(74,2,92,69));dj_legend(c('Gate unmet','Covered; not tested'),c('#777777','#777777'),48.5,11,c('cross','hollow'),width=86);d<-S('members')
  region('B_G',c(17,24,28,61));region('B_H',c(48,24,80,61));dj_label('Myeloid',22.5,20,'bold',.5)
  for(j in 1:4)dj_label(c('Myel.','Fibro.','B/pl.','Epith.')[j],48+(j-.5)*8,20,hjust=.5)
  for(i in 1:7){y<-26+(i-1)*5.2;dj_label(G[i],15,y+1,'italic',1);z<-dj_one(d,d$display_program=='G'&d$display_gene==G[i]);dj_heat(z,22.5,y,10.7,5,'effect',c(-.8,.8))}
  for(i in 1:4){y<-27+(i-1)*10;dj_label(H[i],46,y+1,'italic',1);for(j in 1:4){x<-48+(j-.5)*8;if(j==4){rect(x,y,7.7,9.6,'#EEEEEE');dj_label('□',x,y+1,hjust=.5);cov<-S('coverage');cov<-dj_one(cov,cov$compartment=='epithelial'&cov$program=='H');stopifnot(cov[[paste0('member_',H[i],'_available')]]=='TRUE');mark_record('not_tested',cov$binding_key,x,y,list(member=H[i],field=paste0('member_',H[i],'_available'),value='TRUE',state='NA_NOT_APPLICABLE',reason='Covered member; not in original native-member inferential roster'))}else{z<-dj_one(d,d$display_program=='H'&d$display_gene==H[i]&d$display_cell_type==c('myeloid','fibroblast','B_plasma')[j]);dj_heat(z,x,y,7.7,9.6,'effect',c(-.8,.8))}}}
  dj_colorbar('B',85,26,58,c(-.8,.8),'Δlog₂')
  }
  if(dj_panel_available('C')) {
  panel('C','Epithelial H sensitivity',c(2,79,164,85));a<-S('rank');a<-a[a$record_type=='paired_rank_scores',];a<-a[order(a$donor,method='radix'),];b<-S('normalized');b<-b[match(a$donor,b$donor),];stopifnot(identical(a$donor,b$donor),nrow(a)==23)
  dj_label('GSE150290 · 23 paired donors',84,11,hjust=.5)
  for(j in 1:2){d<-if(j==1)a else b;e<-if(j==1){z<-S('rank');z[z$record_type=='accepted_primary_mean_CI',]}else S('normalized_result');top<-23+(j-1)*30;lim<-dj_limits(c(d$delta,e$ci_lower,e$ci_upper));r<-region(paste0('C',j),c(20,top,160,top+20),c(.5,28.5),lim);zero_y(r);dj_label(c('Rank H','Normalized H4')[j],20,top-3,'bold');dj_label('Mean (95% CI)',r$x(26.5),top-3,hjust=.5);yaxis(r,dj_ticks(lim,3),if(j==1)'Δ rank score' else 'Δ mean member-z',title_x=3)
    for(i in 1:23){point_xy(d[i,],r,i,num(d$delta[i]),'', 'delta',radius=.65,fill=if(j==1)blue else orange);line(r$x(i),r$y(0),r$x(i),r$y(num(d$delta[i])),if(j==1)blue else orange,.16)}
    dj_vci(e,r,26.5,FALSE,if(j==1)blue else orange)
    if(j==2)for(i in 1:23)dj_label(sub('Pat','',d$donor[i]),r$x(i),top+25,hjust=.5)
  };dj_label('Donor ID',82,85,hjust=.5)
  }
  if(dj_panel_available('D')) {
  panel('D','Epithelial member changes',c(2,170,164,44));dj_legend(c('Overall epithelium (23)','Enteroendocrine (12)','Unavailable'),c(blue,orange,'#777777'),87.5,11,c('circle','circle','cross'));d<-S('epi_members')
  for(j in 1:4){left<-13+(j-1)*41;r<-region(paste0('D',j),c(left,21,left+26,36),c(.5,2.5),c(-.9,.1));if(j==2)dj_zero_gaps(r,'y',r$x(2)) else zero_y(r);dj_label(H[j],left+13,18,'italic',.5);yaxis(r,c(-.8,-.4,0),show_labels=j==1,title=if(j==1)'Δ log₂ expression' else '',title_x=2)
    for(i in 1:2){z<-dj_one(d,d$endpoint==c('overall_epithelial','enteroendocrine')[i]&d$member==H[j]);dj_vci(z,r,i,FALSE,c(blue,orange)[i]);dj_label(c('Overall','EEC')[i],r$x(i),42,hjust=.5)}}
  }
  list(rows=list(c('B_G','B_H'),paste0('D',1:4)),cols=list(c('C1','C2')),ex=list(list(panels=c('B_G','B_H'),checks=c('panel-width','horizontal-gutter'),reason='Approved one-column G matrix beside four-column H matrix')))
}

dj_fig8 <- function() {
  cc<-c('T_NK','B_plasma','myeloid','fibroblast','endothelial','mast_other_immune','epithelial');labs<-c('T/NK','B/plasma','Myeloid','Fibroblast','Endothelial','Mast/other immune','Epithelial')
  if(dj_panel_available('A')) {
  panel('A','Captured cell composition',c(2,2,164,96))
  dj_legend(labs,unname(cell_colors[cc]),90.5,11,types=rep('rect',7),nrow=2,width=141,broad=TRUE);d<-S('composition');donors<-sort(unique(d$donor),method='radix');stopifnot(length(donors)==24,nrow(d)==336)
  for(j in 1:2){top<-26+(j-1)*33;r<-region(paste0('A',j),c(20,top,161,top+25),c(.5,24.5),c(0,100));yaxis(r,c(0,50,100),if(j==1)'Captured cells (%)' else '',title_x=3);dj_label(c('Adjacent','Tumor')[j],20,top-2,'bold')
    for(i in 1:24){z<-d[d$donor==donors[i]&d$tissue==c('NAT','Tumor')[j],];z<-z[match(cc,z$compartment),];stopifnot(nrow(z)==7,abs(sum(num(z$proportion))-1)<1e-12);base<-0
      for(k in 1:7){v<-num(z$proportion[k])*100;cell_mark(z[k,],r$x(i),r$y(base+v/2),5.3,abs(r$y(v)-r$y(0)),v,'proportion',c(0,100),'stacked_fraction',cell_colors[cc[k]],list(region=r$id,factor=100,donor=donors[i],tissue=c('NAT','Tumor')[j],lower=base,upper=base+v));base<-base+v}
      if(j==2)dj_label(i,r$x(i),top+30,hjust=.5)}
  };dj_label('Donor',92,96,hjust=.5)
  }
  if(dj_panel_available('B')) {
  panel('B','Method sensitivity',c(2,104,100,76));dj_label('Empirical logit:',1,11);dj_legend(c('100 cells','200 cells'),c(blue,orange),67,11,c('circle','square'),width=61);d<-S('methods')
  for(j in 1:3){z<-d[d$composition_method==c('empirical_logit','CLR','quasibinomial')[j],];if(!nrow(z))z<-d[grepl(c('empirical','CLR','quasibinomial')[j],d$primary_method),];stopifnot(nrow(z)==c(14,7,7)[j]);lim<-dj_limits(c(z$ci_lower,z$ci_upper));left<-32+(j-1)*23;r<-region(paste0('B',j),c(left,27,left+19,63),lim);zero_x(r);dj_label(c('Empirical logit\nn = 20 / 18','CLR\nn = 20','Quasibinomial\nn = 24')[j],left+9.5,23,'bold',.5)
    for(i in 1:7){y<-30+(i-1)*5;zz<-z[z$compartment==cc[i],];if(j==1)dj_label(labs[i],30,y+1,hjust=1);for(k in seq_len(nrow(zz))){strict<-grepl('200',zz$role[k]);dj_hci(zz[k,],r,y+if(nrow(zz)==2)if(strict)1 else -1 else 0,FALSE,if(strict)orange else blue,if(strict)'square' else 'circle',.5)}}
    xaxis(r,dj_ticks(lim,2));dj_label(c('Δ logit','Δ CLR','Log odds')[j],left+9.5,74,hjust=.5)
  }
  }
  if(dj_panel_available('C')) {
  panel('C','Enteroendocrine fractions',c(106,104,60,76));dj_legend(c('Adjacent','Tumor'),c(blue,orange),36,11,width=55);d<-S('eec_pairs');dj_label('GSE150290 · n = 20',36,18,hjust=.5)
  r<-region('C',c(15,26,57,70),c(0,1),c(0,max(num(c(d$NAT_proportion,d$Tumor_proportion)))*1.08));yaxis(r,dj_ticks(r$ylim,3),'Fraction of captured epithelium',title_x=3);dj_pair(d,r,'NAT_proportion','Tumor_proportion','donor')
  }
  list(rows=list(paste0('B',1:3)),cols=list(c('A1','A2')),ex=list())
}

