# set up GenomicRanges rect plot
GRPlotBasic = function(gr, main=NULL, xlab=NULL, ylab=NULL, barheight=.8){
  plotdata=GRPlotNew(gr, main,xlab, ylab, barheight)
  # basic ranges/tracks  
  rect(start(plotdata$gr_ranges), 
       seq_along(plotdata$gr_ranges), 
       end(plotdata$gr_ranges), 
       seq_along(plotdata$gr_ranges) + barheight,col='white')
  # basic ranges/tracks
  invisible(plotdata)
}
GRPlotNew = function(gr, main=NULL, xlab=NULL, ylab=NULL, barheight=.8){
  # make the new empty plot, including axes
  if (class(gr) %in% c('CompressedGRangesList','GRangesList')) {
    gr = unlist(gr)
  }
  gr = sort(gr)
  gr_ranges = ranges(gr)
  gr_data = mcols(gr)
  chrom = seqlevelsInUse(gr)
  stopifnot(length(chrom) == 1) # handle multiple chromosomes later
  
  low_end = min(start(gr_ranges))
  high_end = max(end(gr_ranges))
  
  plot(
    x=c(low_end, 
        high_end), 
    y=c(1, 
        length(gr_ranges)+1.5), 
    type='n', 
    axes=F,
    xlab="",
    ylab=sprintf("%d genomic range%s", 
                 length(gr), 
                 ifelse(length(gr)>1,"s","")),
    main=sprintf("%s:%s-%s", 
                 chrom, 
                 commathou(low_end),
                 commathou(high_end)))
  box()
  
  ## axis
  majTx = 10^floor( log10(high_end-low_end)) # ~8000 => 10^3
  minTx = majTx/2
  
  round_factor = floor(log10(minTx)) # 500 => 2
  majTx_start = round(low_end, digits = -round_factor)
  
  if (majTx_start > low_end) { majTx_start = majTx_start - 10^round_factor }
  
  majorTicks = seq(majTx_start, high_end, by=majTx)
  minorTicks = seq(majTx_start+minTx, high_end, by=majTx)
  
  micTx = ifelse(minTx <= 5, 1, (minTx/2))
  microTicks = seq(majTx_start+micTx, 
                   high_end, 
                   by=ifelse(minTx<=5,
                             1,
                             minTx))
  
  # make sure the ticks don't stop before the data
  if (high_end > tail(majorTicks,1) ) {
    majorTicks = seq(majTx_start, high_end+majTx, by=majTx)
  }
  if (high_end > tail(minorTicks,1)) {
    minorTicks = seq(majTx_start+minTx, high_end+majTx, by=majTx)
  }
  
  # the top of the graph is reserved for the major, minor, microticks
  topBar = length(gr) + 1
  ybars = seq(topBar, topBar+.5, length.out=3)
  
  segments(majorTicks[1],ybars[3], majorTicks[1] + majTx)
  text(majorTicks[1] + majTx, ybars[3], sprintf("%d bp", majTx),pos=4,cex=.75, offset=0)
  segments(majorTicks[1],ybars[2], majorTicks[1] + minTx)
  text(majorTicks[1] + minTx, ybars[2], sprintf("%d bp", minTx),pos=4,cex=.75, offset=0)
  segments(majorTicks[1],ybars[1], majorTicks[1] + micTx)
  text(majorTicks[1] + micTx, ybars[1], sprintf("%d bp", micTx),pos=4,cex=.75, offset=0)
  
  
  axis(1,
       at= majorTicks, 
       labels = commathou(majorTicks),
       cex=.5,las=2)
  abline(v=majorTicks, lty=1, col="grey", lwd=.5)
  
  axis(1,
       at= minorTicks, 
       labels = FALSE,
       cex=.5,las=2,lwd=.5)
  abline(v=minorTicks, lty=3, col="grey")
  
  axis(1,
       at= microTicks, 
       labels = FALSE,
       cex=.5,las=2,lwd=.25)
  abline(v=microTicks, lwd=.5, lty=3, col="grey")
  
  return( 
    list('gr_data'=gr_data,
         'gr_ranges'=gr_ranges,
         'barheight'=barheight)
    )
}

plotIDR = function(gr_ranges, gr_data, barheight){

  # basic ranges/tracks  
  rect(start(gr_ranges), 
       seq_along(gr_ranges), 
       end(gr_ranges), 
       seq_along(gr_ranges) + barheight,col='white')
  # basic ranges/tracks
  
  rect(gr_data$peak_1_start, 
       seq_along(gr_ranges), 
       gr_data$peak_1_end, 
       seq_along(gr_ranges) + barheight, 
       col=rgb(0,0,1,.1),lwd=0)
  
  rect(gr_data$peak_2_start, 
       seq_along(gr_ranges), 
       gr_data$peak_2_end, 
       seq_along(gr_ranges) + barheight, 
       col=rgb(1,0,0,.1),lwd=0)
  
  text( mid(gr_ranges), seq_along(gr_ranges)+ barheight/2, 
        cex=.5, 
        labels = gr_data$called.stage )
  
  summits = start(gr_ranges) + gr_data$overall_summit_off
  segments(summits, y0=seq_along(gr_ranges), y1=seq_along(gr_ranges) + barheight)
  # summits = gr_data$peak_1_start + gr_data$peak_1_summit
  # segments(summits, y0=seq_along(gr_ranges), y1=seq_along(gr_ranges) + barheight, col='blue')
  # summits = gr_data$peak_2_start + gr_data$peak_2_summit
  # segments(summits, y0=seq_along(gr_ranges), y1=seq_along(gr_ranges) + barheight, col='red')
  # 
}