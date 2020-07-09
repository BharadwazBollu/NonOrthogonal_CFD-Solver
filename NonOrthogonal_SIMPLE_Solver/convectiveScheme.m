function [ue,uw,un,us,ve,vw,vn,vs,convecion_CentralCoeff] = convectiveScheme(convScheme,i,j,Nx,Ny,flux_e,flux_w,flux_n,flux_s,u,v)

convecion_CentralCoeff = 0 ;

if ( convScheme == 1 )
    
    
    if ( flux_e(i,j) >= 0)         % Upwind scheme for Convected term
        ue = u(i,j);
        ve = v(i,j);
        convecion_CentralCoeff = convecion_CentralCoeff + flux_e(i,j);
    else
        ue = u(i+1,j);
        ve = v(i+1,j);
    end
    
    if ( flux_w(i,j) >= 0)
        uw = u(i,j);
        vw = v(i,j);
        convecion_CentralCoeff = convecion_CentralCoeff + flux_w(i,j);
    else
        uw = u(i-1,j);
        vw = v(i-1,j);
    end
    
    if ( flux_n(i,j) >= 0)
        un = u(i,j);
        vn = v(i,j);
        convecion_CentralCoeff = convecion_CentralCoeff + flux_n(i,j);
    else
        un = u(i,j+1);
        vn = v(i,j+1);
    end
    
    if ( flux_s(i,j) >= 0)
        us = u(i,j);
        vs = v(i,j);
        convecion_CentralCoeff = convecion_CentralCoeff + flux_s(i,j);
    else
        us = u(i,j-1);
        vs = v(i,j-1);
    end
    
end

if ( convScheme == 2 )
    if ( flux_e(i,j) >= 0)
        ue = u(i,j);
        ve = v(i,j);
        convecion_CentralCoeff = convecion_CentralCoeff + flux_e(i,j);
    else
        ue = u(i+1,j);
        ve = v(i+1,j);
    end
    
    if ( flux_w(i,j) >= 0)
        uw = u(i,j);
        vw = v(i,j);
        convecion_CentralCoeff = convecion_CentralCoeff + flux_w(i,j);
    else
        uw = u(i-1,j);
        vw = v(i-1,j);
    end
    
    if ( flux_n(i,j) >= 0)
        un = u(i,j);
        vn = v(i,j);
        convecion_CentralCoeff = convecion_CentralCoeff + flux_n(i,j);
    else
        un = u(i,j+1);
        vn = v(i,j+1);
    end
    
    if ( flux_s(i,j) >= 0)
        us = u(i,j);
        vs = v(i,j);
        convecion_CentralCoeff = convecion_CentralCoeff + flux_s(i,j);
    else
        us = u(i,j-1);
        vs = v(i,j-1);
    end
    
    if ( i > 2 && j > 2 )
        if ( i < Nx && j < Ny )
            
            if ( flux_e(i,j) >= 0)
                ue = 6.0/8*u(i,j) + 3.0/8*u(i+1,j) - 1.0/8*u(i-1,j);
                ve = 6.0/8*v(i,j) + 3.0/8*v(i+1,j) - 1.0/8*v(i-1,j);
            else
                ue = 6.0/8*u(i+1,j) + 3.0/8*u(i,j) - 1.0/8*u(i+2,j);
                ve = 6.0/8*v(i+1,j) + 3.0/8*v(i,j) - 1.0/8*v(i+2,j);
            end
            if ( flux_w(i,j) >= 0)
                uw = 6.0/8*u(i-1,j) + 3.0/8*u(i,j) - 1.0/8*u(i-2,j);
                vw = 6.0/8*v(i-1,j) + 3.0/8*v(i,j) - 1.0/8*v(i-2,j);
            else
                uw = 6.0/8*u(i,j) + 3.0/8*u(i-1,j) - 1.0/8*u(i+1,j);
                vw = 6.0/8*v(i,j) + 3.0/8*v(i-1,j) - 1.0/8*v(i+1,j);
            end
            if ( flux_n(i,j) >= 0)
                un = 6.0/8*u(i,j) + 3.0/8*u(i,j+1) - 1.0/8*u(i,j-1);
                vn = 6.0/8*v(i,j) + 3.0/8*v(i,j+1) - 1.0/8*v(i,j-1);
            else
                un = 6.0/8*u(i,j+1) + 3.0/8*u(i,j) - 1.0/8*u(i,j+2);
                vn = 6.0/8*v(i,j+1) + 3.0/8*v(i,j) - 1.0/8*v(i,j+2);
            end
            if ( flux_s(i,j) >= 0)
                us = 6.0/8*u(i,j-1) + 3.0/8*u(i,j) - 1.0/8*u(i,j-2);
                vs = 6.0/8*v(i,j-1) + 3.0/8*v(i,j) - 1.0/8*v(i,j-2);
            else
                us = 6.0/8*u(i,j) + 3.0/8*u(i,j-1) - 1.0/8*u(i,j+1);
                vs = 6.0/8*v(i,j) + 3.0/8*v(i,j-1) - 1.0/8*v(i,j+1);
            end
        end
    end
    
end

end