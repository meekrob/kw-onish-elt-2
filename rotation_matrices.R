Rx = function(theta) {
  return(matrix(c(1,0,0,0,
           0, cos(theta), -sin(theta), 0,
           0, sin(theta), cos(theta),0,
           0, 0, 0, 1), nrow=4, byrow = T));
}

Ry = function(theta) {
  return(matrix(c(cos(theta),0,sin(theta),0,
                  0, 1, 0, 0,
                  -sin(theta), 0, cos(theta),0,
                  0, 0, 0, 1), nrow=4, byrow = T));
}

Rz = function(theta) {
  return(matrix(c
                (cos(theta),-sin(theta),0,0,
                  cos(theta),sin(theta), 0, 0,
                  0, 0, 1, 0,
                  0, 0, 0, 1), nrow=4, byrow = T));
}

flatten_std = function(p) {
  return(p %*% Rx(pi/2) %*% Rz(pi/6));
}

polar_to_xy = function(r, theta) {
  x = r * cos(theta);
  y = r * sin(theta);
  return(cbind(x,y));
}

xy_to_polar = function(x,y) {
  r = sqrt(x^2 + y^2);
  theta = atan(y/x);
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
