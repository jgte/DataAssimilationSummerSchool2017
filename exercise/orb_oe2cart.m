function out = orb_oe2cart(in,miu,method)
% ORB_OE2CART(IN) transforms input orbit defined in classic orbital elements into
% the corresponding cartesian ones.
%
%   ORB_OE2CART(IN,MIU,METHOD) allows non-standard gravitational constant and
%   selecting the desired method.
%
% See low-level function 'oe_oe2cart' for details on the available methods.

% Created by J.Encarnacao <J.deTeixeiradaEncarnacao@tudelft.nl>
% Improvements by P.Inacio <p.m.g.inacio@tudelft.nl>

% load orbit if necessary
if ischar(in), in = orb_read(in); end

% check input
if ~isorb(in)
    error('orb_cart2oe:input','Input orbit must be valid orbit object')
end
% NOTE: This is a low level function used to convert the type of an orbit.
%       It is called by orb_setType. It does not make sense to accept an
%       orbit that is not in orbital elements.
if ~strcmp(in.Type,'oe')
    error('orb_cart2oe:metadata','Input must be defined in orbital elements.')
end
if ~exist('miu','var') || isempty(miu)
    %notice: the linear units (m, km, cm, etc) that this routine expects
    %and returns are the same as the ones in which <miu> is defined, so
    %ignore the commends in each subroutine regarding this
    miu=orbit_default_const;
end
if ~exist('method','var') || isempty(method)
    method = 4;
end

%% Get data
%oe=orb_dyn_input(in,'oe');
oe=in.Data;

%% Convert data
xyz = zeros(size(oe))
% copy time vector
xyz(:,1) = oe(:,1);
% convert coordinates
xyz(:,2:end) = oe_oe2cart(oe(:,2:end),miu,method);

% build output object
warning('orb_cart2oe:implementation','TODO: Handle Units field of orbit');
out = orb_constructor(xyz,'ReferenceDate',in.ReferenceDate,...
        'ReferenceFrame',in.ReferenceFrame,'Comment',in.Comment,...
        'Type','cart');
