function error = eulerIntegrator(x,z,u,dt)

error = (z - x) - u*dt;