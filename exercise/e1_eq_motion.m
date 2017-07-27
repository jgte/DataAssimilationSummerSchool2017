% y(1:3)         : cartesian position in the inertial frame
% y(4:6)=dy(1:3) : cartesian velocity in the inertial frame
% dy(4:6)        : cartesian acceleration in the inertial frame
function dy = e1_eq_motion(~,y)
  % get global constants
  global G M_moon
  
  % compute the relative vector to the center of mass of the moon
  dy=y(1:3)-[0 0 0]';
  % compute distance to the center of the moon
  r=sqrt(sum(dy.^2));
  % gravity acceleration unit vector
  g_uv = - dy/r;
  % compute gravitational acceleration
  g = G*M_moon/(r^2)*g_uv;

  % assign outputs
  dy = [y(4:6);g];
end