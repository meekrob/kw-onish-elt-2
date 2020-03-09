
row_func = function(r) {
  # cluster matrix C must exist outside scope of this function
  # k must be the last element of vector r
  n = length(r);
  x = r[1:(n-1)];
  k = r[n];
  return(weights_f = function(c) {
    weight_f(c,x,C[k,]);
    apply(C,2,weights_f)->weights;
    return(1/sum(weights));
  });
}

weight_f = function(c,v,ck,m=2) {
  exponent = 2/(m-1);
  quot = euclid(v,ck) / euclid(v,c);
  return(quot^exponent);
}

euclid = function(u,v) {
  return(c(dist(rbind(u,v))));
}