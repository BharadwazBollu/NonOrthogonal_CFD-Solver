function [error,u_curr,v_curr] = rmsError(u_curr,v_curr,u_Next,v_Next,Nx,Ny)


error = 0;

for  j = 1 : Ny+1
    for  i = 1 : Nx+1
        error = error + (u_Next(i,j) - u_curr(i,j))^2 + (v_Next(i,j) - v_curr(i,j))^2 ;     % Residue for difflux_erence in current and next time step values
        u_curr(i,j)  = u_Next(i,j);                                                              % updating values for next time step
        v_curr(i,j)  = v_Next(i,j);
    end
end

error = sqrt(error/(Nx*Ny));

end