#### Run the markdown file far enough to set include3d_x and kmeans_on3d_k4
stopifnot('include3d_x' %in% ls())
stopifnot('kmeans_on3d_k4' %in% ls())

library(plotly)

source('kweights.R')
C=kmeans_on3d_k4$centers
mx = cbind(include3d_x, kmeans_on3d_k4$cluster)
kweights = apply(mx, 1, row_func)

plottable = data.frame(include3d_x, k=as.factor(kmeans_on3d_k4$cluster),
                       knames=kmeans_on3d_k4$cluster,kweights=kweights)
knames=kmeans_on3d_k4$cluster
plottable$knames[which(knames == 1)] <- 'L3'
plottable$knames[which(knames == 2)] <- 'Larval'
plottable$knames[which(knames == 3)] <- 'LE'
plottable$knames[which(knames == 4)] <- 'Increasing'
plottable$knames = factor(plottable$knames, levels=c("LE","Larval","L3","Increasing"))

trendcolors=c("#D95F02","#1B9E77","#7570B3","#E7298A")
trendcolors=c("#7570B3","#1B9E77", "#D95F02","#E7298A")
plot_ly(plottable, x=~LE, y=~L1, z=~L3, type='scatter3d',
        mode='markers',
        color=~knames, 
        size=~kweights,
        colors=trendcolors)
