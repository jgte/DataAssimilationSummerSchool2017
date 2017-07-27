function phi=ang_fix_pi(phi)

%Makes sure phi is within [-pi,pi].
%Input <phi> can be of any size and shape.

%tranforming angles to [-pi,pi] domain
idx = abs(phi) > pi;
phi(idx) = rem(phi(idx),2*pi);

idx = abs(phi) > pi;
phi(idx) = phi(idx) - 2*sign(phi(idx))*pi;

%bug trap
if any(abs(phi(:)) > pi)
    error([mfilename,': bug trap: any(abs(phi(i)) > pi)'])
end
