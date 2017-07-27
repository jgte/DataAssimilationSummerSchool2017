%% Objectives:
%  - plot differences betweeen disturbed and undisturbed orbits

%% house-keeping

clear all %#ok<CLSCR>
close all

%% define constants

global G M_moon omega
G         = 6.67408e-11;           %universal gravitational constant 
M_moon    = 7.342e22;              %mass of the moon [kg]
R         = 1737.4e3;              %mean radius [m]
d_moon    = 27.321661*24*3600;     %moon day [s]
dt        = 10;                    %time step [s]
plot_scale= 600;                   %make the plots this big [pixels]
h_orb     = 100e3;                 %orbit altitude [m]

%dependent constants
omega = 2*pi/(d_moon);            %angular speed [rad/s]
t     = transpose(0:dt:d_moon/2); %time domain (why d_moon/2 for end time?)
f_args=struct(...
  'square',{{'Position',50+[0,0,  1,  1]*plot_scale,'PaperUnits','points','PaperPosition',50+[0,0,1  ,1  ]*plot_scale}},...
  'rectg', {{'Position',50+[0,0,2.1,0.9]*plot_scale,'PaperUnits','points','PaperPosition',50+[0,0,2.1,0.9]*plot_scale}}...
);

%% define the mass anomaly

global mass_anomaly_loc mass_anomaly_mag
%location of the mass anomaly (lon,lat,r)
[x,y,z]=sph2cart(pi,0,R);
mass_anomaly_loc=[x,y,z];
%magnitude of the mass anomaly [kg]
mass_anomaly_mag=1e12;

%% define initial conditions of orbit

y0=oe_oe2cart([...
  R+h_orb, ... %      a    = Semi-major axis       [m]
  0,       ... %      e    = eccentricity          [ ]
  pi/2,    ... %      i    = Inclination           [rad]
  pi/2,    ... %      Omega= Right-Ascension of the ascending node [rad]
  0,       ... %      omega= Argument of periapse (angle between the ascending node and the pericenter) [rad]
  0        ... %      nu   = True anomaly          [rad]
],G*M_moon);

%% integrate the orbit 

[~,y1]=orbit_integrate(t,y0,@e1_eq_motion);
[t,y2]=orbit_integrate(t,y0,@e2_eq_motion);

%% convert from inertial to co-rotating reference frame

[y_cr1,sph_cr1]=celestial2corotating(t,y1,'orbit');
[y_cr2,sph_cr2]=celestial2corotating(t,y2,'orbit');

%% some info plots

figure(f_args.rectg{:})
subplot(3,3,1)
plot(t,sph_cr1(:,1)*180/pi)
title('Longitude (no mass anomaly)');ylabel('[deg]')
subplot(3,3,2)
plot(t,sph_cr2(:,1)*180/pi)
title('Longitude (with mass anomaly)');ylabel('[deg]')
subplot(3,3,3)
plot(t,(sph_cr1(:,1)-sph_cr2(:,1))*180/pi*3600)
title('Longitude residual');ylabel('[arcsec]')

subplot(3,3,4)
plot(t,sph_cr1(:,2)*180/pi)
title('Latitude (no mass anomaly)');ylabel('[deg]')
subplot(3,3,5)
plot(t,sph_cr2(:,2)*180/pi)
title('Latitude (with mass anomaly)');ylabel('[deg]')
subplot(3,3,6)
plot(t,(sph_cr1(:,2)-sph_cr2(:,2))*180/pi*3600)
title('Latitude residual');ylabel('[arcsec]')

subplot(3,3,7)
plot(t,sph_cr1(:,3)-R-h_orb)
title(['Radius difference to ',num2str(h_orb*1e-3),'km (no mass anomaly)]']);ylabel('[m]')
subplot(3,3,8)
plot(t,sph_cr2(:,3)-R-h_orb)
title(['Radius difference to ',num2str(h_orb*1e-3),'km (with mass anomaly)']);ylabel('[m]')
subplot(3,3,9)
plot(t,sph_cr1(:,3)-sph_cr2(:,3))
title('Radius residual')

%% plots of the residuals in the spatial domain

figure(f_args.rectg{:})
title_str={'Longitude','Latitude','Radial'};
for i=1:3
  subplot(1,3,i)
  %compute this residuals
  res=sph_cr1(:,i)-sph_cr2(:,i);
  %plot this residuals
  scatter3(y_cr1(:,1),y_cr1(:,2),y_cr1(:,3),ones(size(i)),res,'Marker','.'); hold on
  %show the location of the mass anomaly
  plot3(mass_anomaly_loc(1),mass_anomaly_loc(2),mass_anomaly_loc(3),'r','MarkerSize',20,'Marker','.','MarkerFaceColor','r');
  %format and anotate the plot
  axis equal
  title([title_str{i},' residuals'])
  colormap jet
  colorbar
end
