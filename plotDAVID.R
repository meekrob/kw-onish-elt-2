# Plot Data from DAVID similar to their web output, which is not exportable


plotDavid = function(anno_cluster.wide, 
                     main.label="Some Representative Term", 
                     chip.cluster.label="k cluster",
                     xaxis.size=6,
                     yaxis.size=6)
{
  anno_cluster.wide$X.1 = NULL # trailing tabs add a garbage column
  anno_cluster.wide$gene = factor(rownames(anno_cluster.wide), levels=rownames(anno_cluster.wide), ordered=T)
  
  anno_cluster.long = melt(anno_cluster.wide)
  green_label="corresponding gene-term association positively reported"
  black_label="corresponding gene-term association not reported yet"
  labelled_value = factor(ifelse(anno_cluster.long$value, green_label, black_label),levels=c(black_label,green_label),ordered=T)
  anno_cluster.long$value = labelled_value
  
  tileplot<-ggplot(anno_cluster.long, aes(x=gene, y=variable)) + 
    ggtitle(
      bquote(
        paste(bold(.(chip.cluster.label)), 
              ": DAVID Cluster including ", 
              italic(.(main.label)) 
        )
      )
    ) +
    xlab("Gene name") +
    ylab("Annotation term") +
    geom_tile(aes(fill=value),color="white") + 
    scale_fill_manual(values=c("black","green")) + 
    coord_equal() +
    theme(
      legend.position="top", legend.title = element_blank(), legend.direction = "vertical",
      axis.text.x = element_text(angle = 45,hjust=1,size=xaxis.size),
      axis.text.y = element_text(size = yaxis.size)) 
  return(tileplot)
}