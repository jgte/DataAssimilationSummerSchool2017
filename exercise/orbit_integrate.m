function [t,y]=orbit_integrate(t,y0,eq_motion)
  % set integration options
  options = odeset('RelTol',2.5e-14,'AbsTol',1e-8);
  % call the integrator
  [t1,y] = ode113(eq_motion,t,y0,options);
  % make sure the time domain didn't change
  assert(all(t1==t),'the ode113 integrator changed the time domain.')
end
