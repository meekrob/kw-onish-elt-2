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

source('rotation_matrices.R')

# begin plotting on a 2D plane
flattened = flatten_std(include3d_x)
flattened_polar = xy_to_polar(flattened)
flat_polar_half = flat_polar    
flat_polar_half[,'r'] <- .3

plot(polar_to_xy(flattened_polar),col=trendcolors[kmeans_on3d_k4$cluster],asp=1,type='n')
# get the clusters individually
rings=250
spread=1.005
kweight_size = 2
transparency = .5
magenta = flat_polar_half[kmeans_on3d_k4$cluster==4,]
pointsPerArc=nrow(magenta)/rings
magenta_spread = spread_out_by_radius(magenta,pointsPerArc,spread)
clr = c(col2rgb(trendcolors[4])/255,transparency)
points(polar_to_xy(magenta_spread[,1],
                   magenta_spread[,2]),
                  col=rgb(clr[1],clr[2],clr[3],clr[4]),
                  cex=kweights[kmeans_on3d_k4$cluster==4]/kweight_size,
                  pch=15)

purple = flat_polar_half[kmeans_on3d_k4$cluster==1,]
pointsPerArc=nrow(purple)/rings
purple_spread = spread_out_by_radius(purple,pointsPerArc,spread)
clr = c(col2rgb(trendcolors[1])/255,transparency)
points(polar_to_xy(purple_spread[,1],
                   purple_spread[,2]),
       col=rgb(clr[1],clr[2],clr[3],clr[4]),
       cex=kweights[kmeans_on3d_k4$cluster==1]/kweight_size,
       pch=15)

green = flat_polar_half[kmeans_on3d_k4$cluster==2,]
pointsPerArc=nrow(green)/rings
green_spread = spread_out_by_radius(green,pointsPerArc,spread)
clr = c(col2rgb(trendcolors[2])/255,transparency)
points(polar_to_xy(green_spread[,1],
                   green_spread[,2]),
       col=rgb(clr[1],clr[2],clr[3],clr[4]),
       cex=kweights[kmeans_on3d_k4$cluster==2]/kweight_size,
       pch=15)

orange = flat_polar_half[kmeans_on3d_k4$cluster==3,]
pointsPerArc=nrow(orange)/rings
orange_spread = spread_out_by_radius(orange,pointsPerArc,spread)
clr = c(col2rgb(trendcolors[3])/255,transparency)
points(polar_to_xy(orange_spread[,1],
                   orange_spread[,2]),
       col=rgb(clr[1],clr[2],clr[3],clr[4]),
       cex=kweights[kmeans_on3d_k4$cluster==3]/kweight_size,
       pch=15)
