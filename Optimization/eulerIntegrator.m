function error = eulerIntegrator(x,z,u,dt)
% Euler integration scheme

% (x_(k+1) - x_k)/dt = xdot ; with x_(k+1) = z and u = xdot
error = (z - x) - u*dt;