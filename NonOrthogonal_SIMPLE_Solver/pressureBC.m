function [p_Next] = pressureBC(p_Next,Nx,Ny,mesh_num)
% BC for p_Next

if (mesh_num == 1)
    for i = 1:Nx+1
        p_Next(i,Ny+1)  =  p_Next(i,Ny) ;           % Top moving wall nuemon
    end
    
    for i = 1:Nx+1
        p_Next(i,1)     =  p_Next(i,2)  ;        	% Bottom wall nuemon
    end
    
    for j = 1 : Ny+1
        p_Next(1,j)     =  p_Next(2,j)  ;        	% left wall nuemon
    end
    
    for j = 1 : Ny+1
        p_Next(Nx+1,j)  =  p_Next(Nx,j) ;           % right wall nuemon
    end
elseif(mesh_num == 2)
    
    for i = 1:Nx+1
        p_Next(i,Ny+1)  =  p_Next(i,Ny) ;           % Top wall nuemon
    end
    
    for i = 1:Nx+1
        p_Next(i,1)     =  p_Next(i,2)  ;        	% Bottom wall nuemon
    end
    
    for j = 1 : Ny+1
        p_Next(1,j)     =  p_Next(2,j)  ;        	% left inlet nuemon
    end
    
    for j = 1 : Ny+1
        p_Next(Nx+1,j)  =  0 ;                      % right outlet fixed value
    end
elseif(mesh_num == 3)
    for j = 1:ceil(ceil(Ny/2)/4)
        p_Next(Nx+1,j)  =  p_Next(Nx,j) ;               % Inlet top half / free stream
    end
    
    for j = (Nx+1-ceil(ceil(Ny/2)/4) +1 ):Nx+1
        p_Next(Nx+1,j)  =  p_Next(Nx,j) ;               % Inlet bottom half/ free stream
    end
    
    temp = ( Nx+1 - 4*(ceil(ceil(Ny/2)/4))-2 )/2 ;
    
    for j = ceil(ceil(Ny/2)/4)+1:(ceil(ceil(Ny/2)/4)+temp+2)
        p_Next(Ny+1,j)  =  p_Next(Ny,j) ;               % top nuemon
    end
    
    for j = (ceil(ceil(Ny/2)/4)+temp+2 + 1):(Nx+1-ceil(ceil(Ny/2)/4)-temp -2)
        p_Next(Nx+1,j)  =  0              ;             % Outlet fixed value
    end
    
    for j = (Nx+1-ceil(ceil(Ny/2)/4)-temp-1):(Nx+1-ceil(ceil(Ny/2)/4)  )
        p_Next(Ny+1,j)  =  p_Next(Ny,j) ;               % bottom nuemon
    end
    
    for j = 2 : Nx
        p_Next(1,j)     =  p_Next(2,j)    ;             % circle nuemon
    end
    
    for i = 1 : Ny+1
        p_Next(i,1)     =  0    ;                       % cutline  
    end
    
    for i = 1 : Ny+1
        p_Next(i,Nx+1)  =  0    ;                       % cutline
    end
    
end

end