function oe=oe_fix_domain(oe,disp_flag)
% OUT=OE_FIX_DOMAIN(IN,DISP_FLAG) fixes the angular coordinates of the list of
% orbital elements IN to the interval [0,2pi]
%
%   IN is a Nx6 matrix with the orbital elements.
%
%   See also OE_OE2CART, OE_OE2SPH.

% Created by J.Encarnacao <J.deTeixeiradaEncarnacao@tudelft.nl>

% check inputs
if ~exist('disp_flag','var') || isempty(disp_flag)
    disp_flag=true;
end
if ~isnumeric(oe) || size(oe,2) ~= 6
  error('oe_fix_domain:input','Input in must be a Nx6 list of orbital elements.')
end

%converting to radians if needed
oe(:,3:6)=ang_get_radians(oe(:,3:6));

%checking orbital elements
if disp_flag
    out = [    fix_domain(oe(:,1),0,inf,'a'),...
               fix_domain(oe(:,2),0,1,  'e'),...
           ang_fix_domain(oe(:,3),0,pi, 'i'),...
           ang_fix_2pi(   oe(:,4)),...
           ang_fix_2pi(   oe(:,5)),...
           ang_fix_2pi(   oe(:,6))];
else
    out = [    fix_domain(oe(:,1),0,inf),...
               fix_domain(oe(:,2),0,1),...
           ang_fix_domain(oe(:,3),0,pi),...
           ang_fix_2pi(   oe(:,4)),...
           ang_fix_2pi(   oe(:,5)),...
           ang_fix_2pi(   oe(:,6))];
end

