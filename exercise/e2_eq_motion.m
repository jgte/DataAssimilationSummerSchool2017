% y(1:3)         : cartesian position in the inertial frame
% y(4:6)=dy(1:3) : cartesian velocity in the inertial frame
% dy(4:6)        : cartesian acceleration in the inertial frame
function dy = e2_eq_motion(t,y)
  % get global constants
  global G M_moon mass_anomaly_loc mass_anomaly_mag
  
  % compute the relative vector to the center of mass of the moon
  dy=y(1:3)-0;
  % compute distance to the center of the moon
  r=sqrt(sum(dy.^2));
  % gravity acceleration unit vector
  g_uv = - dy/r;
  % compute gravitational acceleration
  g = G*M_moon/(r^2)*g_uv;
  
  % rotate the planet under the orbit (faster than rotating the orbit)
  ma_y=celestial2corotating(t,mass_anomaly_loc,'planet');
  % compute relative vector to the mass anomaly
  ma_dy=y(1:3)-ma_y(:);
  % compute the distance to the mass anomaly
  ma_r=sqrt(sum(ma_dy.^2));
  % compute the unit vector of the acceleration caused by the mass anomaly
  ma_g_uv= - ma_dy/ma_r;
  % compute gravitational acceleration caused by the mass anomaly
  ma_g=G*mass_anomaly_mag/(ma_r^2)*ma_g_uv;
  
  % assign outputs
  dy = [y(4:6);g+ma_g];
end