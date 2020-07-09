function [p_Next,pdiff_e,pdiff_w,pdiff_n,pdiff_s] = solvePressure(Nx,Ny,Vp,dt,u_pred,v_pred,p_Next,surface_e_x,surface_w_x,surface_n_x,surface_s_x,...
    surface_e_y,surface_w_y,surface_n_y,surface_s_y,normCoff_e,normCoff_w,normCoff_n,normCoff_s,crossCoff_e,crossCoff_w,crossCoff_n,crossCoff_s,...
    pdiff_e,pdiff_w,pdiff_n,pdiff_s)

iter = 0 ;
residual_error  = 1 ;

while ( residual_error > 1e-07 )
    
    iter = iter + 1 ;
    pressure_residual = 0;
    
    for  j = 2 : Ny
        for i = 2 : Nx
            % volume interpolation for corner points for pressure
            [p_Next_ne,p_Next_nw,p_Next_se,p_Next_sw] = cornerPointsVolumeInterpolation(i,j,Vp,p_Next) ;
            % calculating diffusion terms for each face for pressure
            [pdiff_e(i,j),pdiff_w(i,j),pdiff_n(i,j),pdiff_s(i,j)] = diffusion(i,j,p_Next,p_Next_ne,p_Next_nw,p_Next_se,p_Next_sw,....
                normCoff_e,normCoff_w,normCoff_n,normCoff_s,crossCoff_e,crossCoff_w,crossCoff_n,crossCoff_s) ;
            
            % central coefficient for pressure poisson equation
            pNext_CentralCoeff = normCoff_e(i,j) + normCoff_w(i,j) + normCoff_n(i,j) + normCoff_s(i,j) ;
            
            % total diffsion for presuure
            pdiff  =  pdiff_e(i,j) + pdiff_w(i,j) + pdiff_n(i,j) + pdiff_s(i,j) ;
            
            % calculating flux without pressure correction for each face
            flux_noPressureCorrection_e =   ( Vp(i+1,j)*u_pred(i,j) + Vp(i,j)*u_pred(i+1,j) )/(Vp(i,j)+Vp(i+1,j)) * surface_e_x(i,j)...
                +  ( Vp(i+1,j)*v_pred(i,j) + Vp(i,j)*v_pred(i+1,j) )/(Vp(i,j)+Vp(i+1,j)) * surface_e_y(i,j)  ;
            flux_noPressureCorrection_w =   ( Vp(i-1,j)*u_pred(i,j) + Vp(i,j)*u_pred(i-1,j) )/(Vp(i,j)+Vp(i-1,j)) * surface_w_x(i,j)...
                +  ( Vp(i-1,j)*v_pred(i,j) + Vp(i,j)*v_pred(i-1,j) )/(Vp(i,j)+Vp(i-1,j)) * surface_w_y(i,j)  ;
            flux_noPressureCorrection_n =   ( Vp(i,j+1)*u_pred(i,j) + Vp(i,j)*u_pred(i,j+1) )/(Vp(i,j)+Vp(i,j+1)) * surface_n_x(i,j)...
                +  ( Vp(i,j+1)*v_pred(i,j) + Vp(i,j)*v_pred(i,j+1) )/(Vp(i,j)+Vp(i,j+1)) * surface_n_y(i,j)  ;
            flux_noPressureCorrection_s =   ( Vp(i,j-1)*u_pred(i,j) + Vp(i,j)*u_pred(i,j-1) )/(Vp(i,j)+Vp(i,j-1)) * surface_s_x(i,j)...
                +  ( Vp(i,j-1)*v_pred(i,j) + Vp(i,j)*v_pred(i,j-1) )/(Vp(i,j)+Vp(i,j-1)) * surface_s_y(i,j)  ;
            
            totalflux_noPressureCorrection = flux_noPressureCorrection_e + flux_noPressureCorrection_w...
                + flux_noPressureCorrection_n + flux_noPressureCorrection_s ;
            
            p_Next_error        =  totalflux_noPressureCorrection - dt * pdiff ;
            pressure_residual   =  pressure_residual + p_Next_error * p_Next_error;
            p_Next(i,j)         = -p_Next_error/( dt*pNext_CentralCoeff ) + p_Next(i,j);
        end
    end
    residual_error = sqrt(pressure_residual/(Nx*Ny));                    
end

fprintf('GaussSeidalSolver:  Solving for P, Final residual = %e, No Iterations %d\n',residual_error,iter) ;

end