function phi=ang_fix_2pi(phi)

%Makes sure phi is within [0,2*pi].
%Input <phi> can be of any size and shape.

%transforming to [-pi,pi] domain
phi=phi-pi;

%fixing phi to [-pi,pi] range
phi=ang_fix_pi(phi);

%transforming back to [0,2*pi] domain
phi=phi+pi;
