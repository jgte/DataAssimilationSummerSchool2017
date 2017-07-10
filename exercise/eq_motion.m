% y(1:3)         : cartesian position in the inertial frame
% y(4:6)=dy(1:3) : cartesian velocity in the inertial frame
% dy(4:6)        : cartesian acceleration in the inertial frame
function dy = eq_motion(~,y)
  % get global constants
  global miu
  % spherical coordinates
  [~,~,r] = cart2sph(y(1),y(2),y(3));
  % gravity acceleration unit vector
  g_uv = - y(1:3)/r;
  % compute gravitational accelerations
  g = miu/(r^2)*g_uv;
  % assign outputs
  dy = [y(4:6);g];
end