%% Objectives:
%  - integrate undisturbed orbit

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

%dependent constants
omega = 2*pi/(d_moon);            %angular speed [rad/s]
t     = transpose(0:dt:d_moon/2); %time domain (why d_moon/2 for end time?)
f_args=struct(...
  'square',{{'Position',50+[0,0,  1,  1]*plot_scale,'PaperUnits','points','PaperPosition',50+[0,0,1  ,1  ]*plot_scale}},...
  'rectg', {{'Position',50+[0,0,2.1,0.9]*plot_scale,'PaperUnits','points','PaperPosition',50+[0,0,2.1,0.9]*plot_scale}}...
);

%% define initial conditions of orbit

y0=oe_oe2cart([...
  R+100e3, ... %      a    = Semi-major axis       [m]
  0,       ... %      e    = eccentricity          [ ]
  pi/2,    ... %      i    = Inclination           [rad]
  pi/2,    ... %      Omega= Right-Ascension of the ascending node [rad]
  0,       ... %      omega= Argument of periapse (angle between the ascending node and the pericenter) [rad]
  0        ... %      nu   = True anomaly          [rad]
],G*M_moon);

%% integrate the orbit 

[t,y]=orbit_integrate(t,y0,@e1_eq_motion);

%% plot the orbit (the moon rotates under the orbit!)

figure(f_args.square{:})
%plot the orbit
plot3(y(:,1),y(:,2),y(:,3));
hold on
%plot the moon
[mx,my,mz]=sphere(100);
mx=mx*R;
my=my*R;
mz=mz*R;
surf(mx,my,mz,...
  'EdgeColor','none',...
  'FaceColor','texturemap',...
  'Cdata',flipud(imread('ear1ccc2.jpg'))...
);
colormap gray
axis equal
title('The moon rotates under the orbit!','FontSize',18,'Color','red')
grid on

%% convert from inertial to co-rotating reference frame

y_cr=celestial2corotating(t,y,'orbit');

%% plot the ground tracks

figure(f_args.square{:})
%plot the orbit
plot3(y_cr(:,1),y_cr(:,2),y_cr(:,3),'b');
hold on
%plot the moon
[mx,my,mz]=sphere(100);
mx=mx*R;
my=my*R;
mz=mz*R;
h=surf(mx,my,mz,...
  'EdgeColor','none',...
  'FaceColor','texturemap',...
  'Cdata',flipud(imread('ear1ccc2.jpg'))...
);
colormap gray
axis equal
set(gca, 'visible', 'off') ;
