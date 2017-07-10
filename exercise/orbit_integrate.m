function [t,y]=orbit_integrate(y0,t_step,t_end)
  % build the time domain
  tspan = 0:t_step:t_end;
  % set integration options
  options = odeset('RelTol',1e-08);
  % call the integrator
  [t,y] = ode45(@eq_motion,tspan,y0,options);
end
