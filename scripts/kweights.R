# unnormed gotten by cbinding 3 random vectors
# normed: mean subtracted, sd divided
# normedk = data.frame(normed,k=k$cluster)
# normed_weights = apply(normedk, 1, row_func)
# plot_ly(rbind(unnormed,rep(1,3),rep(-1,3),rep(0,3),centers), x=~x,y=~y,z=~z,type="scatter3d", mode='markers',color=as.factor(c(k$cluster,0,0,0,4,5,6)), colors=c('#666666',   '#0000FF','#00FF00','#FF0000','#000088','#008800','#880000'), size= 5*c(normed_weights,rep(1,6)))

row_func = function(r) {
  # cluster matrix C must exist outside scope of this function
  # k must be the last element of vector r
  n = length(r);
  x = r[1:(n-1)];
  k = r[n];
  weights_f = function(c) {
    return(weight_f(c,x,C[k,]));
  }
  apply(C,1,weights_f)->weights;
  return(1/sum(weights));
}

weight_f = function(c,v,ck,m=2) {
  exponent = 2/(m-1);
  quot = euclid(v,ck) / euclid(v,c);
  return(quot^exponent);
}

euclid = function(u,v) {
  return(c(dist(rbind(u,v))));
}

random_kmeans = function(n, data, k) {
 # you should set.seed before calling this function
  nrows = nrow(data);
  centers = data.frame()
  xserial = 1:ncol(data)
  for (draw in 1:n) {
     clustering = kmeans( data[sample(1:nrows),],k)
     slopes = c()
     for (j in 1:k) {
        slopes = c(slopes,lm(clustering$centers[j,] ~ xserial)$coef[2])
     }
     clustering$centers = cbind(clustering$centers, slopes)
     trial=data.frame(clustering$centers, trial=draw,k=1:k );
     centers = rbind(centers,trial);
  }
  return(centers);
}