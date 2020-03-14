Rx = function(theta) {
  return(matrix(c(1,0,0,0,
           0, cos(theta), -sin(theta), 0,
           0, sin(theta), cos(theta),0,
           0, 0, 0, 1), nrow=4, byrow = T));
}

Ry = function(theta) {
  return(matrix(c( cos(theta),  0, sin(theta),   0,
                       0     ,  1,     0     ,   0,
                  -sin(theta),  0, cos(theta),   0,
                       0     ,  0,      0    ,   1), nrow=4, byrow = T));
}

Rz = function(theta) {
  return(matrix(c
                ( cos(theta),-sin(theta),  0, 0,
                  sin(theta), cos(theta),  0, 0,
                      0     ,     0     ,  1, 0,
                      0     ,     0     ,  0, 1), nrow=4, byrow = T));
}

ux=-sqrt(1/3)
uy=ux
uz=ux


flatten_std = function(o) {
  p = as.matrix(cbind(o,1))
  tformed = p %*% Rz(-pi/4)# %*% Ry(-pi/4) %*% Rz(-pi/4);
  colnames(tformed) <- colnames(p);
  return(as.data.frame(tformed[,1:3]));
}

polar_to_xy = function(r, theta=NULL) {
  if (is.null(theta)) {
    if (is.null(dim(r))) {
      stop("r must have at least 2 columns if theta is null.")  
    }
    else if (ncol(r) > 1) {
      theta = r[,2]
      r = r[,1]
    }
    else {
      stop("r must have at least 2 columns if theta is null.")  
    } 
  }
  x = r * cos(theta);
  y = r * sin(theta);
  return(cbind(x,y));
}

xy_to_polar = function(x,y=NULL) {
  if (is.null(y)) {
    if (is.null(dim(x))) {
      stop("x must have at least 2 columns if y is null.")  
    }
    else if (ncol(x) > 1) {
    y = x[,2]
    x = x[,1]
    }
    else {
      stop("x must have at least 2 columns if y is null.")  
    } 
  }
  r = sqrt(x^2 + y^2);
  theta = atan2(y,x);
  return(cbind(r,theta));
}

updateXY = function(df) {

    xy = polar_to_xy(df$r,df$theta)
    df$x = xy[,'x']
    df$x = xy[,'y']
    return(df);
  
}

updatePolar = function(df) {
  polar = xy_to_polar(df$x,df$y)
  df$r = polar[,'r'];
  df$theta = polar[,'theta'];
  return(df);
}

getArcSeries = function(r, startTheta, endTheta, n=10) {
  angles = seq(startTheta,endTheta,length.out=n);
  xy = polar_to_xy( rep(r,n), angles)
  return(xy)
}

ring = function(offset, thickness, startTheta, endTheta, n=10) {
  innerPoints = getArcSeries(offset, startTheta,endTheta,n);
  outerPoints = getArcSeries(offset+thickness, endTheta,startTheta,n); # in reverse order
  return(rbind(innerPoints,outerPoints));
}

radToDegrees = function(radians) {
  return((180*radians/pi ) %% 360)
}

degreesToRadians = function(degrees) {
   return(((degrees %% 360)/360) * (2*pi))
}

polar_segments = function(r0,t0,r1,t1=t0,...) {
  maxn = max(length(r0), length(t0), length(r1),length(t1))
  r0 = rep_len(r0,maxn)
  t0 = rep_len(t0,maxn)
  r1 = rep_len(r1,maxn)
  t1 = rep_len(t1,maxn)
  starts = polar_to_xy(r0,t0);
  ends = polar_to_xy(r1,t1);
  segments(starts[,1],starts[,2],ends[,1],ends[,2],...);
  
}
spread_out_by_radius<-function(polar, len, fact)
{
  i = 1:nrow(polar)
  polar[,'r'] = polar[,'r']* fact ^ round(i/len)
  return(polar);
}