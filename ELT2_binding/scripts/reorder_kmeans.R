# renumber kmeans output by rank of cluster size, low (1) - high (k)

reorder_kmeans = function(kobj) {
  n = kobj$size # vector of cluster sizes
  k = length(n) # scalar
  oldnumbers = seq(1,k) # 1 thru k
  newnumbers = rank(kobj$size)
  tbl = cbind(oldnumbers,newnumbers)
  # renumbered data 
  kobj$cluster = tbl[kobj$cluster,2]
  kobj$centers = kobj$centers[order(newnumbers),]
  kobj$size = kobj$size[order(newnumbers)]
  kobj$withinss = kobj$withinss[order(newnumbers)]
  return(kobj)
}