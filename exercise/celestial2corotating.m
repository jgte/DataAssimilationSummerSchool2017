%NOTE: only position is converted
function [y_corotating,sph_corotating]=celestial2corotating(t,y_inertial,mode)

  %get relevant constants
  global omega

  %need spherical coordinates
  [lon_inertial,lat_inertial,r_inertial] = cart2sph(y_inertial(:,1),y_inertial(:,2),y_inertial(:,3));
  %co-rotating frame adds/subtracts a constant angular rate to the longitude
  switch mode
  case 'planet'
    lon_corotating=lon_inertial-omega.*t;
  case 'orbit'
    lon_corotating=lon_inertial+omega.*t;
  end
  %moke room for outputs
  y_corotating=zeros(size(y_inertial(:,1:3)));
  %convert back to cartesian
  [y_corotating(:,1),y_corotating(:,2),y_corotating(:,3)]=sph2cart(lon_corotating,lat_inertial,r_inertial);
  %assign additional outputs
  sph_corotating=[wrapToPi(lon_corotating),lat_inertial,r_inertial];

end