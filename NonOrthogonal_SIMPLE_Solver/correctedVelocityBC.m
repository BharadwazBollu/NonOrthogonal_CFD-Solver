function [u_Next,v_Next] = correctedVelocityBC(u_Next,v_Next,Nx,Ny,mesh_num)
% BC for Corrected velocities

if (mesh_num == 1)
    for i = 1:Nx+1
        u_Next(i,Ny+1)  =  1 ;          % Top moving wall
        v_Next(i,Ny+1)  =  0 ;
    end
    
    for i = 1:Nx+1
        u_Next(i,1)     =  0 ;          % Bottom wall
        v_Next(i,1)     =  0 ;
    end
    
    for j = 1 : Ny+1
        u_Next(1,j)     =  0 ;          % left wall
        v_Next(1,j)     =  0 ;
    end
    
    for j = 1 : Ny+1
        u_Next(Nx+1,j)  =  0 ;          % right wall
        v_Next(Nx+1,j)  =  0 ;
    end
elseif(mesh_num == 2)
    for i = 1:Nx+1
        u_Next(i,Ny+1)  =  0 ;                      % Top no slip
        v_Next(i,Ny+1)  =  0 ;
    end
    
    for i = 1:Nx+1
        u_Next(i,1)     =  0 ;                      % Bottom no slip
        v_Next(i,1)     =  0 ;
    end
    
    for j = 1 : Ny+1
        u_Next(1,j)     =  1 ;                      % left inlet
        v_Next(1,j)     =  0 ;
    end
    
    for j = 1 : Ny+1
        u_Next(Nx+1,j)  =  0 ;                      % right outlet 
        v_Next(Nx+1,j)  = -1 ;
    end
elseif(mesh_num == 3)
    
    for j = 1:ceil(ceil(Ny/2)/4)
        u_Next(Nx+1,j)  =  1  ;
        v_Next(Nx+1,j)  =  0  ;                         % Inlet top / free stream
    end
    
    for j = (Nx+1-ceil(ceil(Ny/2)/4) +1 ):Nx+1
        u_Next(Nx+1,j)  =  1  ;
        v_Next(Nx+1,j)  =  0  ;                         % Inlet bottom / free stream
    end
    
    temp = ( Nx+1 - 4*(ceil(ceil(Ny/2)/4))-2 )/2 ;
    
    for j = ceil(ceil(Ny/2)/4)+1:(ceil(ceil(Ny/2)/4)+temp+2)
        u_Next(Ny+1,j)  =  1 ;                          % free stream top
        v_Next(Ny+1,j)  =  0 ;
    end
    
    for j = (ceil(ceil(Ny/2)/4)+temp+2 + 1):(Nx+1-ceil(ceil(Ny/2)/4)-temp -2)
        u_Next(Nx+1,j)  =  u_Next(Nx,j) ;               % Outlet nuemon
        v_Next(Nx+1,j)  =  v_Next(Nx,j) ;
    end
    
    for j = (Nx+1-ceil(ceil(Ny/2)/4)-temp-1):(Nx+1-ceil(ceil(Ny/2)/4)  )
        u_Next(Ny+1,j)  =  1 ;                          % free stream bottom
        v_Next(Ny+1,j)  =  0 ;
    end
    
    for j = 1 : Nx+1
        u_Next(1,j)     =  0    ;                    	% circle no slip
        v_Next(1,j)     =  0    ;
    end
    
    for i = 2 : Ny
        u_Next(i,1)     =  1    ;                       % cutline free stream
        v_Next(i,1)     =  0    ;
    end
    
    for i = 2 : Ny
        u_Next(i,Nx+1)  =  1    ;                       % cutline free stream
        v_Next(i,Nx+1)  =  0    ;
    end
    
end

end