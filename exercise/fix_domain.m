function out = fix_domain(in,lower,upper,string)

%makes sure <in> is lower than <upper> and higher than <lower>. If not, the
%closest of the two is returned. If <in> does not violated the <upper> and
%<lower> bounds, <out = in>.
%
%The optional input <string> is used to pass the name of the quantity being
%checked, in which case a warning message is displayed.

if ~isscalar(lower) || ~isscalar(upper)
    error([mfilename,': inputs <lower> and <upper> must be scalar.'])
end

out=zeros(size(in));

idx_u = in > upper;
if any(idx_u)
    out(idx_u) = upper;
    if (nargin == 4) && ~isempty(string)
        if isscalar(in)
            disp([mfilename,': ',string,'=',num2str(in),' is ilegal because upper bound is ',num2str(upper),'.'])
        else
            disp([mfilename,': ',string,' is ilegal at ',num2str(sum(idx_u)),' entries because upper bound is ',num2str(upper),'.'])
        end
    end
end
idx_l = in < lower;
if any(idx_l)
    out(idx_l) = lower;
    if (nargin == 4) && ~isempty(string)
        if isscalar(in)
            disp([mfilename,': ',string,'=',num2str(in),' is ilegal because lower bound is ',num2str(lower),'.'])
        else
            disp([mfilename,': ',string,' is ilegal at ',num2str(sum(idx_l)),' entries because lower bound is ',num2str(lower),'.'])
        end
    end
end
out(~idx_u & ~idx_l) = in(~idx_u & ~idx_l);

