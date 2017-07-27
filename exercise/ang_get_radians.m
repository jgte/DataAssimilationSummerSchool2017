function out=ang_get_radians(in,force_flag)
% ANG_GET_RADIANS(IN,FORCE_FLAG) if needed, converts the input IN from
% degrees into radians. The criteria used to see if IN is in degrees is if
% any of its elements are outside -2*pi and 2*pi. For this reason, if IN is
% in degrees and ranges between around -6.3 to 6.3 degrees, then it will be
% mistaken for radians and not converted.
%
% The criteria can be overriden with FORCE_FLAG

% Created by J.Encarnacao <J.deTeixeiradaEncarnacao@tudelft.nl>

if ~exist('force_flag','var') || isempty(force_flag)
    force_flag=false;
end

if  force_flag || any(abs(in(~isnan(in))) > 2*pi)
    out=in*pi/180;
else
    out=in;
end
