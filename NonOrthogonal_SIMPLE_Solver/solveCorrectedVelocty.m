function [u_Next,v_Next,flux_e,flux_w,flux_n,flux_s] = solveCorrectedVelocty(Nx,Ny,Vp,Re,dt,u_curr,v_curr,u_pred,v_pred,p_Next,u_Next,v_Next,surface_e_x,surface_w_x,surface_n_x,surface_s_x,...
    surface_e_y,surface_w_y,surface_n_y,surface_s_y,normCoff_e,normCoff_w,normCoff_n,normCoff_s,crossCoff_e,crossCoff_w,crossCoff_n,crossCoff_s,...
    pdiff_e,pdiff_w,pdiff_n,pdiff_s,convScheme)

iter = 0;
residual_error  = 1 ;

while ( residual_error > 1e-07 )
    uNext_residual = 0;
    vNext_residual = 0;
    
    iter = iter + 1 ; 
    [flux_e,flux_w,flux_n,flux_s] = flux(Nx,Ny,Vp,u_pred,v_pred,surface_e_x,surface_w_x,surface_n_x,surface_s_x,surface_e_y,surface_w_y,surface_n_y,surface_s_y,dt,pdiff_e,pdiff_w,pdiff_n,pdiff_s) ;
    
    
    for j = 2 : Ny
        for i = 2 : Nx
            % upwind function
            [ue,uw,un,us,ve,vw,vn,vs,convecion_CentralCoeff] = convectiveScheme(convScheme,i,j,Nx,Ny,flux_e,flux_w,flux_n,flux_s,u_Next,v_Next) ;
            % total convection for each cell
            uConvection = flux_e(i,j) * ue + flux_w(i,j) * uw + flux_n(i,j) * un + flux_s(i,j) * us ;
            vConvection = flux_e(i,j) * ve + flux_w(i,j) * vw + flux_n(i,j) * vn + flux_s(i,j) * vs ;
            
            % volume interpolation for corner points for u
            [u_Next_ne,u_Next_nw,u_Next_se,u_Next_sw] = cornerPointsVolumeInterpolation(i,j,Vp,u_Next) ;
            % volume interpolation for corner points for v
            [v_Next_ne,v_Next_nw,v_Next_se,v_Next_sw] = cornerPointsVolumeInterpolation(i,j,Vp,v_Next) ;
            
            % calculating diffusion terms for each face for u corrected velocity
            [udiff_e,udiff_w,udiff_n,udiff_s] = diffusion(i,j,u_Next,u_Next_ne,u_Next_nw,u_Next_se,u_Next_sw,....
                normCoff_e,normCoff_w,normCoff_n,normCoff_s,crossCoff_e,crossCoff_w,crossCoff_n,crossCoff_s) ;
            % central coefficient for u corrected velocity diffusion term
            uNextDiff_CentralCoeff = normCoff_e(i,j) + normCoff_w(i,j) + normCoff_n(i,j) + normCoff_s(i,j) ;
            % total diffsion for u corrected velocity
            udiff  =  udiff_e + udiff_w + udiff_n + udiff_s  ;
            
            % calculating diffusion terms for each face for v corrected velocity
            [vdiff_e,vdiff_w,vdiff_n,vdiff_s] = diffusion(i,j,v_Next,v_Next_ne,v_Next_nw,v_Next_se,v_Next_sw,....
                normCoff_e,normCoff_w,normCoff_n,normCoff_s,crossCoff_e,crossCoff_w,crossCoff_n,crossCoff_s) ;
            % central coefficient for v corrected velocity diffusion term
            vNextDiff_CentralCoeff  = normCoff_e(i,j) + normCoff_w(i,j) + normCoff_n(i,j) + normCoff_s(i,j) ;
            % total diffsion for v corrected velocity
            vdiff  =  vdiff_e + vdiff_w + vdiff_n + vdiff_s  ;
            
            % Using gauss divergence theorem integral dp/dn*dV = Sigma(pfsfn)
            % where f is face and n is direction(x,y,z)
            PeSex = ( Vp(i+1,j)*p_Next(i,j) + Vp(i,j)*p_Next(i+1,j) )/(Vp(i,j)+Vp(i+1,j)) * surface_e_x(i,j)  ;
            PwSwx = ( Vp(i-1,j)*p_Next(i,j) + Vp(i,j)*p_Next(i-1,j) )/(Vp(i,j)+Vp(i-1,j)) * surface_w_x(i,j)  ;
            PnSnx = ( Vp(i,j+1)*p_Next(i,j) + Vp(i,j)*p_Next(i,j+1) )/(Vp(i,j)+Vp(i,j+1)) * surface_n_x(i,j)  ;
            PsSsx = ( Vp(i,j-1)*p_Next(i,j) + Vp(i,j)*p_Next(i,j-1) )/(Vp(i,j)+Vp(i,j-1)) * surface_s_x(i,j)  ;
            
            SigmaPfSfx = PeSex + PwSwx + PnSnx + PsSsx ;
            
            PeSey = ( Vp(i+1,j)*p_Next(i,j) + Vp(i,j)*p_Next(i+1,j) )/(Vp(i,j)+Vp(i+1,j)) * surface_e_y(i,j)  ;
            PwSwy = ( Vp(i-1,j)*p_Next(i,j) + Vp(i,j)*p_Next(i-1,j) )/(Vp(i,j)+Vp(i-1,j)) * surface_w_y(i,j)  ;
            PnSny = ( Vp(i,j+1)*p_Next(i,j) + Vp(i,j)*p_Next(i,j+1) )/(Vp(i,j)+Vp(i,j+1)) * surface_n_y(i,j)  ;
            PsSsy = ( Vp(i,j-1)*p_Next(i,j) + Vp(i,j)*p_Next(i,j-1) )/(Vp(i,j)+Vp(i,j-1)) * surface_s_y(i,j)  ;
            
            SigmaPfSfy = PeSey + PwSwy + PnSny + PsSsy ;
            
            % total central coefficient for u and v corrected velocities
            uNextTotal_CentralCoeff = 1 + dt/Vp(i,j) * convecion_CentralCoeff + dt/Re/Vp(i,j) * uNextDiff_CentralCoeff ;
            vNextTotal_CentralCoeff = 1 + dt/Vp(i,j) * convecion_CentralCoeff + dt/Re/Vp(i,j) * vNextDiff_CentralCoeff ;
            
            % calculating error/residue and correcting the value using Gauss Seidal
            uNext_error  	= u_curr(i,j) - u_Next(i,j) - dt/Vp(i,j) * uConvection + dt/Re * udiff/Vp(i,j) - ( dt/Vp(i,j) * SigmaPfSfx );   % residue for predicetd X velocity
            uNext_residual	= uNext_residual + uNext_error * uNext_error;
            u_Next(i,j)   	= uNext_error/uNextTotal_CentralCoeff + u_Next(i,j) ;
            
            vNext_error  	= v_curr(i,j) - v_Next(i,j) - dt/Vp(i,j) * vConvection + dt/Re * vdiff/Vp(i,j) - ( dt/Vp(i,j) * SigmaPfSfy );   % residue for predicetd Y velocity
            vNext_residual	= vNext_residual + vNext_error * vNext_error;
            v_Next(i,j)   	= vNext_error/vNextTotal_CentralCoeff + v_Next(i,j) ;
            
        end
    end
    residual_error = sqrt( (uNext_residual+vNext_residual)/(Nx*Ny));
end

fprintf('GaussSeidalSolver:  Solving for Ux, Final residual = %e, No Iterations %d\n',sqrt(uNext_residual/Nx/Ny),iter) ;
fprintf('GaussSeidalSolver:  Solving for Uy, Final residual = %e, No Iterations %d\n\n',sqrt(uNext_residual/Nx/Ny),iter) ;

end