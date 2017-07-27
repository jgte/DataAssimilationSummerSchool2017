function xyz = oe_oe2cart(oe,miu,method)
% OE_OE2CART transforms input Nx6 matrix of orbital elements into the
% corresponding cartesian ones.
%
%   XYZ=OE_OE2CART(OE) where OE is a 6xN matrix with the following collumns,
%      a    = Semi-major axis       [m]
%      E    = eccentricity          [ ]
%      i    = Inclination           [rad]
%      Omega= Right-Ascension of the ascending node [rad]
%      w    = Argument of periapse  [rad]
%      nu   = True anomaly          [rad]
%
%   XYZ is a 6xN matrix with the x, y and z position and velocity in [m] and
%   [m/s], respectively.
%
%   OE_OE2CART(IN,MIU) allows non-standard gravitational constant.
%
%   See also OE_CART2SPH.

% Created by J.Encarnacao <J.deTeixeiradaEncarnacao@tudelft.nl>

% TODO: Add description of the available methods.

% NOTE: This is a low level function which converts a matrix of values. For
%       orbits use the dedicated wrapper function 'orb_oe2cart'.
% TODO: Remove this error
if isorb(oe)
    error('oe_oe2cart:deprecated','Function ''oe_oe2cart'' is deprecated. Use ''orb_oe2cart'' instead')
end

global stats

%checking inputs
if ~exist('miu','var') || isempty(miu)
    %notice: the linear units (m, km, cm, etc) that this routine expects
    %and returns are the same as the ones in which <miu> is defined, so
    %ignore the commends in each subroutine regarding this
    miu=orbit_default_const;
end
if ~exist('method','var') || isempty(method)
    %Performance comparison between the methods, for converting 8k coordinates:
    %1 - 2.5563s
    %2 - 2.6030s
    %3 - 2.8473s
    %4 - 2.4484s
    %
    %Differences are negligible:
    %       x             	y             	z             	dx/dt         	dy/dt         	dz/dt
    % max  :
    % 1-2  :3.72529e-009  	2.56114e-009  	4.65661e-009  	3.63798e-012  	2.72848e-012  	4.54747e-012
    % 1-3  :1.86265e-009  	1.39698e-009  	1.86265e-009  	0             	0             	0
    % 1-4  :3.72529e-009  	2.56114e-009  	4.65661e-009  	4.09273e-012  	2.72848e-012  	5.00222e-012
    % mean :
    % 1-2  :-1.7346e-011  	4.75702e-012  	1.2643e-011   	-5.10759e-015 	4.30931e-015  	-1.07798e-014
    % 1-3  :3.8064e-012   	-1.89646e-012 	-9.90339e-012 	0             	0             	0
    % 1-4  :-1.66344e-011 	3.82037e-012  	2.38154e-011  	-5.76987e-015 	3.77928e-015  	-2.32666e-014
    % std  :
    % 1-2  :1.15854e-009  	8.60279e-010  	1.2167e-009   	1.13211e-012  	8.74506e-013  	1.50117e-012
    % 1-3  :3.93583e-010  	2.76058e-010  	4.55987e-010  	0             	0             	0
    % 1-4  :1.15908e-009  	8.60403e-010  	1.23082e-009  	1.14707e-012  	8.90041e-013  	1.52076e-012
    method = 4;
end

if method == 0
    stats=cell(4,1);
end

% check inputs
%oe=orb_dyn_input(in,'oe');
if ~isnumeric(oe) || size(oe,2) ~= 6
    error('oe_oe2cart:inputs','Input must be a Nx6 matrix of classic orbital elements.')
end

%checking orbital elements
oe = oe_fix_domain(oe);

% convert the coordinates
xyz = oe2cart_aux(oe,miu,method);

if method == 0
    disp([mfilename,': found discrepancies between the methods available:',10,...
        '      ',tablify(14,'x','y','z','dx/dt','dy/dt','dz/dt'),10,...
        'max  :',10,...
        '1-2  :',tablify(14,max(stats{1}(:,1)),max(stats{1}(:,2)),max(stats{1}(:,3)),max(stats{1}(:,4)),max(stats{1}(:,5)),max(stats{1}(:,6))),10,...
        '1-3  :',tablify(14,max(stats{2}(:,1)),max(stats{2}(:,2)),max(stats{2}(:,3)),max(stats{2}(:,4)),max(stats{2}(:,5)),max(stats{2}(:,6))),10,...
        '1-4  :',tablify(14,max(stats{3}(:,1)),max(stats{3}(:,2)),max(stats{3}(:,3)),max(stats{3}(:,4)),max(stats{3}(:,5)),max(stats{3}(:,6))),10,...
        'mean :',10,...
        '1-2  :',tablify(14,mean(stats{1}(:,1)),mean(stats{1}(:,2)),mean(stats{1}(:,3)),mean(stats{1}(:,4)),mean(stats{1}(:,5)),mean(stats{1}(:,6))),10,...
        '1-3  :',tablify(14,mean(stats{2}(:,1)),mean(stats{2}(:,2)),mean(stats{2}(:,3)),mean(stats{2}(:,4)),mean(stats{2}(:,5)),mean(stats{2}(:,6))),10,...
        '1-4  :',tablify(14,mean(stats{3}(:,1)),mean(stats{3}(:,2)),mean(stats{3}(:,3)),mean(stats{3}(:,4)),mean(stats{3}(:,5)),mean(stats{3}(:,6))),10,...
        'std  :',10,...
        '1-2  :',tablify(14,std(stats{1}(:,1)),std(stats{1}(:,2)),std(stats{1}(:,3)),std(stats{1}(:,4)),std(stats{1}(:,5)),std(stats{1}(:,6))),10,...
        '1-3  :',tablify(14,std(stats{2}(:,1)),std(stats{2}(:,2)),std(stats{2}(:,3)),std(stats{2}(:,4)),std(stats{2}(:,5)),std(stats{2}(:,6))),10,...
        '1-4  :',tablify(14,std(stats{3}(:,1)),std(stats{3}(:,2)),std(stats{3}(:,3)),std(stats{3}(:,4)),std(stats{3}(:,5)),std(stats{3}(:,6)))])
end

%% cart2oe_aux
function xyz = oe2cart_aux(oe,miu,method)

global stats

if size(oe,1) > 1
    xyz=zeros(size(oe,1),6);
    for i=1:size(oe,1)
        xyz(i,:)=oe2cart_aux(oe(i,:),miu,method);
    end
    return
end

if any(isnan(oe)) || sum(oe) == 0
    xyz = nan(1,6);
    return
end

if (oe(2) == 1)
    error([mfilename,': cannot deal with parabolic orbits.'])
end

if (method == 1) || (method == 0)
    [r1,v1] = Kepler_to_RV(oe, miu);  r=r1; v=v1;
end
if (method == 2) || (method == 0)
    [r2,v2] = orb2eci(miu,oe); r=r2; v=v2;
end
if (method == 3) || (method == 0)
    [r3,v3] = oe2rv(oe,miu); r=r3; v=v3;
end
if (method == 4) || (method == 0)
    xyz = dingo(oe,miu);
    r4 = xyz(1:3); v4=xyz(4:6);  r=r4; v=v4;
end

%checking for inconsistencies
if (method == 0)
    if any( abs([r1',v1']-[r2',v2']) > eps ) || ...
       any( abs([r1',v1']-[r3',v3']) > eps ) || ...
       any( abs([r1',v1']-[r4',v4']) > eps )

        stats{1}(end+1,:)=[r1',v1']-[r2',v2'];
        stats{2}(end+1,:)=[r1',v1']-[r3',v3'];
        stats{3}(end+1,:)=[r1',v1']-[r4',v4'];
%         disp([mfilename,': found discrepancies between the methods available; error totals ',...
%               num2str(sum([abs([r1',v1']-[r2',v2']), abs([r2',v2']-[r3',v3']),abs([r3',v3']-[r1',v1'])]))])
    end
end

xyz = [r',v'];

% disp([mfilename,':out:',num2str(size(xyz))])

return

%% Kepler_to_RV
% Name:      Kepler_to_RV.m
% Author:    Jeffrey S. Parker
% Purpose:   A function used to convert an object's
%            Keplerian elements to Cartesian position and velocity
%
% Call:      [Rxyz, Vxyz] = Kepler_to_RV(a, e, i, Omega, w, nu, miu)
%
% Inputs:    a            = Semi-major axis (km)
%            e            = Eccentricity
%            i            = Inclination (rad)
%            Omega        = Right-Ascension of the ascending node (rad)
%            w            = Argument of periapse (rad)
%            nu           = True anomaly (rad)
%            miu           = Gravitational parameter of the system (GM), km^3/s^2
%
% Outputs:   Rxyz         = (x,y,z) Cartesian Position (km)
%            Vxyz         = (vx,vy,vz) Cartesian Velocity (km/s)
%
% Required Subroutines:   N/A

function [Rxyz, Vxyz] = Kepler_to_RV(oe, miu)

a     = oe(1);
e     = oe(2);
i     = oe(3);
Omega = oe(4);
w     = oe(5);
nu    = oe(6);

p = a*(1 - e^2);                    % semiparameter (km)

%%%  Position Coordinates in Perifocal Coordinate System
x  = (p*cos(nu)) / (1 + e*cos(nu)); % x-coordinate (km)
y  = (p*sin(nu)) / (1 + e*cos(nu)); % y-coordinate (km)
z  = 0;                             % z-coordinate (km)
vx = -(miu/p)^(1/2) * sin(nu);       % velocity in x (km/s)
vy = (miu/p)^(1/2) * (e + cos(nu));  % velocity in y (km/s)
vz = 0;                             % velocity in z (km/s)

%%%  Transformation Matrix (3 Rotations)  %%%
ROT = [cos(Omega)*cos(w)-sin(Omega)*sin(w)*cos(i) ...
                (-1)*cos(Omega)*sin(w)-sin(Omega)*cos(w)*cos(i) ...
                            sin(Omega)*sin(i); ...
       sin(Omega)*cos(w)+cos(Omega)*sin(w)*cos(i) ...
                (-1)*sin(Omega)*sin(w)+cos(Omega)*cos(w)*cos(i) ...
                            (-1)*cos(Omega)*sin(i); ...
       sin(w)*sin(i)  cos(w)*sin(i)  cos(i)];

%%%  Transforming Perifocal -> xyz  %%%
Rxyz = ROT*[x y z]';
Vxyz = ROT*[vx vy vz]';

%% orb2eci
function [r, v] = orb2eci(miu, oev)

% convert classical orbital elements to eci state vector

% input

%  oev(1) = semimajor axis (kilometers)
%  oev(2) = orbital eccentricity (non-dimensional)
%           (0 <= eccentricity < 1)
%  oev(3) = orbital inclination (radians)
%           (0 <= inclination <= pi)
%  oev(4) = right ascension of ascending node (radians)
%           (0 <= raan <= 2 pi)
%  oev(5) = argument of perigee (radians)
%           (0 <= argument of perigee <= 2 pi)
%  oev(6) = true anomaly (radians)
%           (0 <= true anomaly <= 2 pi)

% output

%  r = eci position vector (kilometers)
%  v = eci velocity vector (kilometers/second)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = zeros(3, 1);
v = zeros(3, 1);

% unload orbital elements array

sma = oev(1);
ecc = oev(2);
inc = oev(3);
raan = oev(4);
argper = oev(5);
tanom = oev(6);

slr = sma * (1 - ecc * ecc);

rm = slr / (1 + ecc * cos(tanom));

arglat = argper + tanom;

sarglat = sin(arglat);
carglat = cos(arglat);

c4 = sqrt(miu / slr);
c5 = ecc * cos(argper) + carglat;
c6 = ecc * sin(argper) + sarglat;

sinc = sin(inc);
cinc = cos(inc);

sraan = sin(raan);
craan = cos(raan);

% position vector

r(1) = rm * (craan * carglat - sraan * cinc * sarglat);
r(2) = rm * (sraan * carglat + cinc * sarglat * craan);
r(3) = rm * sinc * sarglat;

% velocity vector

v(1) = -c4 * (craan * c6 + sraan * cinc * c5);
v(2) = -c4 * (sraan * c6 - craan * cinc * c5);
v(3) = c4 * c5 * sinc;

%% oe2rv
%  oe2rv.m  Orbital Elements to r,v
%
%  [r,v] = oe2rv(oe,miu)
%			oe = [a e i Om om nu]
%			r,v  expressed in  IJK  frame
%
function [ri,vi] = oe2rv(oe,miu)
	a=oe(1); e=oe(2); i=oe(3); Om=oe(4); om=oe(5); nu=oe(6);
	p = a*(1-e*e);
	r = p/(1+e*cos(nu));
	rv = [r*cos(nu); r*sin(nu); 0];			% in PQW frame
	vv = sqrt(miu/p)*[-sin(nu); e+cos(nu); 0];
%
%	now rotate
%
	cO = cos(Om);  sO = sin(Om);
	co = cos(om);  so = sin(om);
	ci = cos(i);   si = sin(i);
	R  = [cO*co-sO*so*ci  -cO*so-sO*co*ci  sO*si;
		  sO*co+cO*so*ci  -sO*so+cO*co*ci -cO*si;
		  so*si            co*si           ci];
	ri = R*rv;
	vi = R*vv;

%% dingo

function s = dingo(e,gm)

    p = e(1) * ( 1.0 - e(2) * e(2) );
    cv = cos(e(6));
    ecv = 1.0 + e(2) * cv;
    r = p / ecv;
    u = e(5) + e(6);
    cu = cos(u);
    su = sin(u);
    co = cos(e(4));
    so = sin(e(4));
    ci = cos(e(3));
    si = sin(e(3));
    cocu = co * cu;
    sosu = so * su;
    socu = so * cu;
    cosu = co * su;
    fx = cocu - sosu * ci;
    fy = socu + cosu * ci;
    fz = su * si;

    f = sqrt( gm / p );
    vr = f * e(2) * sin(e(6));
    vu = f * ecv;
    s = [r * fx;...
         r * fy;...
         r * fz;...
         vr * fx - vu * (cosu + socu * ci );...
         vr * fy - vu * (sosu - cocu * ci );...
         vr * fz + vu * cu * si];
