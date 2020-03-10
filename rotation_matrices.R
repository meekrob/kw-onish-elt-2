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
  return(p %*% Rz(-pi/2) %*% Ry(pi/2));
}