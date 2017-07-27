%% Objectives:
%  - 

%% house-keeping

clear all %#ok<CLSCR>
close all

%% define constants

global G M_moon omega order
G         = 6.67408e-11;           %universal gravitational constant 
M_moon    = 7.342e22;              %mass of the moon [kg]
R         = 1737.4e3;              %mean radius [m]
d_moon    = 27.321661*24*3600;     %moon day [s]
dt        = 10;                    %time step [s]
plot_scale= 600;                   %make the plots this big [pixels]
h_orb     = 100e3;                 %orbit altitude [m]
order     = 3;                     %order of the central-difference differentiation scheme

%dependent constants
omega = 2*pi/(d_moon);            %angular speed [rad/s]
t     = transpose(0:dt:d_moon/2); %time domain (why d_moon/2 for end time?)
f_args=struct(...
  'square',{{'Position',50+[0,0,  1,  1]*plot_scale,'PaperUnits','points','PaperPosition',50+[0,0,1  ,1  ]*plot_scale}},...
  'rectg', {{'Position',50+[0,0,2.1,0.9]*plot_scale,'PaperUnits','points','PaperPosition',50+[0,0,2.1,0.9]*plot_scale}}...
);

%% define the mass anomaly

global mass_anomaly_mag mass_anomaly_loc
%location of the mass anomaly (lon,lat,r)
[x,y,z]=sph2cart(-pi/2,0,R);
mass_anomaly_loc=[x,y,z];
%magnitude of the mass anomaly [kg]
mass_anomaly_mag=1e12;

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

[t,y]=orbit_integrate(t,y0,@e2_eq_motion);

%% build the design matrix

% compute observations: acceleration
[obs,e]=e4_diff(y(:,1:3),dt,'taylor');

% compute the relative vector to the center of mass of the moon
dy=y(e+1:end-e,1:3);
% compute distance to the center of the moon
r=sqrt(sum(dy.^2,2))*ones(1,3);
% gravity acceleration unit vector
g_uv = - dy./r;
% compute gravitational acceleration
g = G*M_moon./(r.^2).*g_uv;

% compute relative vector to the mass anomaly
ma_dy=y(e+1:end-e,1:3)-...
  celestial2corotating(...
    t(e+1:end-e),...
    ones(numel(t(e+1:end-e)),1)*mass_anomaly_loc,...
    'planet'...
  );
% compute the distance to the mass anomaly
ma_r=sqrt(sum(ma_dy.^2,2))*ones(1,3);
% compute the unit vector of the acceleration caused by the mass anomaly
ma_g_uv= - ma_dy./ma_r;
% compute gravitational acceleration caused by a unitary mass anomaly plus the moon 
ma_g=G./(ma_r.^2).*ma_g_uv;

%this is the design matrix
A=ma_g; 

x=A(:)\(obs(:)-g(:));
disp(['Mass of the anomaly: ',num2str(x),', e=',num2str((x-mass_anomaly_mag)/mass_anomaly_mag*100),'%'])

