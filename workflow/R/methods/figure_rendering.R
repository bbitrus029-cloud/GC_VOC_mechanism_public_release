# Shared rendering primitives for current dual-journal Figures 2–8.
# Rendering only: no scoring, model fitting, resampling, or P/q calculation.
Sys.setenv(OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1, MKL_NUM_THREADS=1,
           RCPP_PARALLEL_NUM_THREADS=1, TZ='UTC', PYTHONDONTWRITEBYTECODE=1)
suppressPackageStartupMessages({library(grid)})
.nature_alignment_groups <- function(groups) lapply(groups, function(g) as.character(unlist(g,use.names=FALSE)))
font_family <- function() if (exists('.F') && !is.null(.F$style)) .F$style$font else 'Tahoma'
num <- function(x) suppressWarnings(as.numeric(x))
mm <- function(x) unit(x,'mm')
ink <- '#202020'; blue <- '#355C70'; negative <- '#A4BBC8'; positive <- '#D9A79B'
individual <- '#8FA6B3'; adjacent <- '#929DA5'; connector <- '#D0D6DA'
light <- '#DCE5EA'; unavailable <- '#ECEDEF'
G <- c('ACSL4','LPCAT3','ALOX5','ALOX12','ALOX15','ALOX15B','POR')
H <- c('ALDH2','ALDH3A1','ALDH3A2','ALDH7A1')
cells <- c('epithelial','fibroblast','endothelial','myeloid','T_NK','B_plasma','mast_other_immune')
cell_labels <- c('Epithelial','Fibroblast','Endothelial','Myeloid','T/NK','B/plasma','Mast/other immune')
cell_colors <- setNames(c('#DDA599','#C3A775','#B8C8D1','#355C70','#7D909D','#A695B4','#B7B7AF'),cells)
cohorts <- c('GSE29272','GSE27342','TCGA-STAD','GSE122401†')
fmt <- function(x,digits=2) sub('-','−',sprintf(paste0('%.',digits,'f'),num(x)),fixed=TRUE)
tick <- function(x,digits=NULL) sub('-','−',if(is.null(digits)) format(x,trim=TRUE,scientific=FALSE) else sprintf(paste0('% .',digits,'f'),x),fixed=TRUE)
tick_labels <- function(x,digits=NULL) trimws(tick(x,digits))
json <- function(x,p) jsonlite::write_json(x,p,auto_unbox=TRUE,pretty=TRUE,na='null',digits=NA)
S <- function(name) .F$sources[[name]]
label <- function(text,x,y,size=8,face='plain',hjust=0,rot=0,col=ink) {
 if (!is.null(.F$style) && exists('dj_label') && grepl('g_av|₂',text)) return(dj_label(text,x,y,face,hjust,rot,col))
 if (!is.null(.F$style) && face != 'bold' && size >= 7 && size < 10) size <- .F$style$label
 stopifnot(size>=5)
 .F$text[[length(.F$text)+1L]]<-list(panel=.F$current,text=text,x_mm=x,y_mm=y,size_pt=size,hjust=hjust,rotation=rot,face=face)
 grid.text(text,x=mm(x),y=unit(1,'npc')-mm(y),just=c(hjust,0),rot=rot,
           gp=gpar(fontfamily=font_family(),fontsize=size,fontface=face,col=col))
}
line <- function(x1,y1,x2,y2,col=ink,width=.12,lty='solid') grid.segments(
 x0=mm(x1),y0=unit(1,'npc')-mm(y1),x1=mm(x2),y1=unit(1,'npc')-mm(y2),
 gp=gpar(col=col,lwd=width*96/25.4,lty=lty,lineend='butt'))
rect <- function(x,y,w,h,fill,border=NA) grid.rect(x=mm(x),y=unit(1,'npc')-mm(y),
 width=mm(w),height=mm(h),gp=gpar(fill=fill,col=border,lwd=.3))
dot <- function(x,y,r=.6,fill=blue,border=NA,shape='circle') {
 if(shape=='square') rect(x,y,2*r,2*r,fill,border)
 else if(shape=='diamond') grid.polygon(x=mm(x+c(0,r,0,-r)),y=unit(1,'npc')-mm(y+c(-r,0,r,0)),gp=gpar(fill=fill,col=border,lwd=.5))
 else grid.circle(x=mm(x),y=unit(1,'npc')-mm(y),r=mm(r),gp=gpar(fill=fill,col=border,lwd=.5))
}
measure <- function() {
 lo<-deviceLoc(unit(0,'npc'),unit(0,'npc'),valueOnly=TRUE)
 hi<-deviceLoc(unit(1,'npc'),unit(1,'npc'),valueOnly=TRUE)
 c(lo$x,lo$y,hi$x,hi$y)*72
}
panel <- function(letter,title,box) {
 upViewport(0);pushViewport(viewport(x=mm(box[1]),y=unit(1,'npc')-mm(box[2]),width=mm(box[3]),height=mm(box[4]),just=c('left','top'),name=paste0('panel_',letter)))
 .F$current<-letter;.F$origin<-box[1:2];.F$panels[[letter]]<-list(id=letter,bbox_pt=measure(),box_mm=box,title=title)
 label(letter,.5,4,10,'bold');label(title,6,4,if(is.null(.F$style))10 else .F$style$title,'bold')
}
region <- function(id,z,xlim=NULL,ylim=NULL) {
 pushViewport(viewport(x=mm(z[1]),y=unit(1,'npc')-mm(z[2]),width=mm(z[3]-z[1]),height=mm(z[4]-z[2]),just=c('left','top'),name=paste0('data_',id)))
 .F$regions[[id]]<-list(id=id,bbox_pt=measure(),panel=.F$current,local_mm=z,xlim=xlim,ylim=ylim)
 upViewport()
 invisible(list(x=function(v) z[1]+(v-xlim[1])/diff(xlim)*(z[3]-z[1]),
      y=function(v) z[4]-(v-ylim[1])/diff(ylim)*(z[4]-z[2]),z=z,xlim=xlim,ylim=ylim,id=id))
}
# Axis records describe only presentation; accepted scientific fields stay untouched.
xaxis <- function(r,ticks,title='',title_y=NULL,size=7,digits=NULL) {
 z<-r$z;line(z[1],z[4],z[3],z[4],'black',.1);labs<-tick_labels(ticks,digits)
 for(i in seq_along(ticks)) {xx<-r$x(ticks[i]);line(xx,z[4],xx,z[4]+1,'black',.1);label(labs[i],xx,z[4]+4.2,size,hjust=.5)}
 .F$axes[[length(.F$axes)+1L]]<-list(panel=.F$current,region=r$id,orientation='x',ticks=ticks,labels=labs,labelled=TRUE)
 if(nzchar(title)) shared_title(title,r$id,if(is.null(title_y))z[4]+9 else title_y)
}
yaxis <- function(r,ticks,title='',size=7,title_x=3,show_labels=TRUE) {
 z<-r$z;line(z[1],z[2],z[1],z[4],'black',.1);labs<-tick_labels(ticks)
 for(i in seq_along(ticks)) {yy<-r$y(ticks[i]);line(z[1]-1,yy,z[1],yy,'black',.1);if(show_labels)label(labs[i],z[1]-2,yy+.85,size,hjust=1)}
 .F$axes[[length(.F$axes)+1L]]<-list(panel=.F$current,region=r$id,orientation='y',ticks=ticks,labels=labs,labelled=show_labels)
 if(nzchar(title)) shared_title(title,r$id,title_x,'y')
}
zero_x <- function(r) if(r$xlim[1]<=0 && r$xlim[2]>=0) {
 line(r$x(0),r$z[2],r$x(0),r$z[4],connector,.12,'dashed')
 .F$zero_lines[[length(.F$zero_lines)+1L]]<-list(region=r$id,orientation='x',segments=list(c(r$x(0),r$z[2],r$x(0),r$z[4])))
}
zero_y <- function(r,gap=NULL) if(r$ylim[1]<=0 && r$ylim[2]>=0) {
 xs<-if(is.null(gap))list(r$z[c(1,3)]) else list(c(r$z[1],gap[1]),c(gap[2],r$z[3]))
 segments<-lapply(xs,function(x) {line(x[1],r$y(0),x[2],r$y(0),connector,.12,'dashed');c(x[1],r$y(0),x[2],r$y(0))})
 .F$zero_lines[[length(.F$zero_lines)+1L]]<-list(region=r$id,orientation='y',segments=segments)
}
shared_title <- function(text,ids,at,orientation='x',kind='plain',size=8) {
 rs<-.F$regions[ids];stopifnot(all(vapply(rs,function(r)r$panel==.F$current,TRUE)))
 z<-do.call(rbind,lapply(rs,function(r)r$local_mm));center<-if(orientation=='x')mean(range(z[,c(1,3)])) else mean(range(z[,c(2,4)]))
 .F$titles[[length(.F$titles)+1L]]<-list(panel=.F$current,regions=ids,orientation=orientation,center_mm=center,at_mm=at,text=text,kind=kind,size_pt=size)
 if(!is.null(.F$style) && exists('dj_label') && grepl('₂|g_av',text)) {if(orientation=='x')dj_label(text,center,at,hjust=.5) else dj_label(text,at,center,hjust=.5,rot=90)} else if(orientation=='x')label(text,center,at,size,hjust=.5) else label(text,at,center,size,hjust=.5,rot=90)
}
mark_record <- function(type,key,x,y,extra=list()) {
 p<-deviceLoc(mm(x),unit(1,'npc')-mm(y),valueOnly=TRUE)
 .F$marks[[length(.F$marks)+1L]]<-c(list(type=type,panel=.F$current,key=key,x_pt=p$x*72,y_pt=p$y*72),extra)
}
heat <- function(v,lim=c(-.8,.8),sequential=FALSE) {
 if(!is.finite(v)) return(unavailable)
 stopifnot(v>=lim[1],v<=lim[2])
 cols<-if(sequential) grDevices::colorRampPalette(c(blue,'#F1F3F5'))(257) else grDevices::colorRampPalette(c(negative,'#F1F3F5',positive))(257)
 cols[round((v-lim[1])/diff(lim)*256)+1]
}
point_xy <- function(d,r,x,y,xfield='',yfield='',dx=0,dy=0,radius=.45,fill=individual) {
 stopifnot(x>=r$xlim[1],x<=r$xlim[2],y>=r$ylim[1],y<=r$ylim[2])
 xx<-r$x(x)+dx;yy<-r$y(y)+dy;dot(xx,yy,radius,fill)
 mark_record('point',d$binding_key,xx,yy,list(region=r$id,x_value=x,y_value=y,x_field=xfield,y_field=yfield,dx_mm=dx,dy_mm=dy,fill=fill,radius_mm=radius))
}
cell_mark <- function(d,x,y,w,h,value,field,limits,kind='cell',fill=heat(value,limits),extra=list()) {
 rect(x,y,w,h,fill)
 mark_record(kind,d$binding_key,x,y,c(list(width_mm=w,height_mm=h,value=value,field=field,limits=limits,fill=fill),extra))
}
