clear all
close all

%% define constants
global miu 
miu   =4.282837e13;              %gravitational parameter (m3/s2)
R     =1737.4e3;                 %mean radius (m)
omega =2*pi/(27.321661*24*3600); %angular speed (rad/s)
plot_scale=1000;
f_args=struct(...
  'square',{{'Position',50+[0,0,  1,  1]*plot_scale,'PaperUnits','points','PaperPosition',50+[0,0,1  ,1  ]*plot_scale}},...
  'rectg', {{'Position',50+[0,0,2.1,0.9]*plot_scale,'PaperUnits','points','PaperPosition',50+[0,0,2.1,0.9]*plot_scale}}...
);

%% define initial conditions of orbit
y0=oe_oe2cart([...
  R+100e3, ... %      a    = Semi-major axis       [m]
  0,       ... %      E    = eccentricity          [ ]
  pi/2,    ... %      i    = Inclination           [rad]
  pi/2,    ... %      Omega= Right-Ascension of the ascending node [rad]
  0,       ... %      w    = Argument of periapse (angle between the ascending node and the pericenter) [rad]
  0        ... %      nu   = True anomaly          [rad]
],miu);

%% integrate the orbit
[t,y]=orbit_integrate(y0,60,3600*24*5);

%% plot the orbit

figure(f_args.square{:})
%plot the orbit
plot3(y(:,1),y(:,2),y(:,3),'.');
hold on
%plot the moon
[mx,my,mz]=sphere(100);
mx=mx*R;
my=my*R;
mz=mz*R;
surf(mx,my,mz,'EdgeColor','none');
axis equal

%% convert from inertial to co-rotating reference frame

%need spherical coordinates
[lon,lat,r] = cart2sph(y(:,1),y(:,2),y(:,3));
%co-rotating frame adds a constant angular rate to the longitude
lon_cr=lon+omega.*t;
%convert back to cartesian
y_cr=zeros(size(y(:,1:3)));
[y_cr(:,1),y_cr(:,2),y_cr(:,3)]=sph2cart(lon_cr,lat,r);

%% plot the ground tracks

figure(f_args.square{:})
%plot the orbit
plot3(y_cr(:,1),y_cr(:,2),y_cr(:,3),'.');
hold on
%plot the moon
[mx,my,mz]=sphere(100);
mx=mx*R;
my=my*R;
mz=mz*R;
surf(mx,my,mz,'EdgeColor','none');
axis equal

%%

figure(f_args.rectg{:})
break_idx=abs(diff(lon_cr))>0.1;
lon_cr_p=lon_cr;lon_cr_p(break_idx)=nan;
lat_p=lat;lat_p(break_idx)=nan;

plot(lon_cr_p*180/pi,lat_p*180/pi), hold on
