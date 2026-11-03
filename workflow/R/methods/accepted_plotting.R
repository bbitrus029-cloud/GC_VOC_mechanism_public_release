# Accepted descriptive plotting functions; saved fits and intervals only.
accepted_draw_tissue_context <- function(x, line, sp, s) {
 # Consume the preserved 100-point fitted-line source, including its accepted CI.
 grid <- line$delta_D
 pred <- as.matrix(line[, c('fit', 'ci_lower', 'ci_upper'), drop = FALSE])
  par(mar=c(5,5,4,1)+.1)
  plot(x$delta_D,x$delta_H,type="n",xlab="Paired D7 change (one score unit; tumor - adjacent)",
   ylab="Paired accepted rank-H change (tumor - adjacent)",
   main=paste0(s," (",nrow(x)," patient pairs)\n","gastric-type mucous and chief-cell expression proxy D7"))
  polygon(c(grid,rev(grid)),c(pred[,2],rev(pred[,3])),col="#DCE5EE",border=NA)
  abline(h=0,col="#777777",lty=3);if(sp$zero_supported[sp$source==s])abline(v=0,col="#777777",lty=3)
  lines(grid,pred[,1],col="#234E70",lwd=2);points(x$delta_D,x$delta_H,pch=16,col="#234E70",cex=.65)
  mtext(sp$intercept_label[sp$source==s],side=3,line=.2,cex=.8)

}

accepted_draw_protein <- function(d) {
 par(mar=c(5,6,4,1)+.1)
 yy<-rev(seq_len(nrow(d)));plot(d$effect,yy,pch=16,xlim=range(c(d$ci_lower,d$ci_upper,0)),ylim=c(.5,nrow(d)+.5),yaxt="n",ylab="",
  xlab="Mean centered log2 tumor/adjacent protein ratio",main="PDC000214 / Mun2019\n80 patients; dependent cross-layer evidence")
 abline(v=0,col="#777777",lty=3);segments(d$ci_lower,yy,d$ci_upper,yy,col="#234E70",lwd=2)
 points(d$effect,yy,pch=16,col="#234E70");axis(2,at=yy,labels=d$family_hypothesis,las=1)

}
