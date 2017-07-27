function out = ang_fix_domain(in,lower,upper,string)
% ANG_FIX_DOMAIN(IN,LOWER,UPPER) is the same as <check_domain>, except that
% the bounds are taken to be connected, much like angular quantities that
% should change smoothly betweem 0 and 360 degrees.
%
%   ANG_FIX_DOMAIN(IN,LOWER,UPPER,STRING) to pass the name of the angle
%   being checked, in which case a warning message is displayed.

% Created by J.Encarnacao <J.deTeixeiradaEncarnacao@tudelft.nl>

if ~exist('string','var')
    string='';
end
if ~isscalar(lower) || ~isscalar(upper)
    error([mfilename,': inputs <lower> and <upper> must be scalar.'])
end

while any(in(:)>upper) || any(in(:)<lower)
    idx = (in(:) > upper);
    if any(idx)
        in(idx) = in(idx) - (upper-lower);
        if ~isempty(string) && isscalar(in)
            disp([mfilename,': ',string,'=',num2str(in),' (was ',num2str(in),...
                ') needed fixing because upper circular bound is ',num2str(upper),'.'])
        end
    end
    idx = (in(:) < lower);
    if any(idx)
        in(idx) = in(idx) + (upper-lower);
        if ~isempty(string) && isscalar(in)
            disp([mfilename,': ',string,'=',num2str(in),' (was ',num2str(in),...
                ') needed fixing because lower circular bound is ',num2str(lower),'.'])
        end
    end
end
out = in;

