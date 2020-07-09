function [u_pred,v_pred] = solvePredictedVelocty(Nx,Ny,Vp,Re,dt,flux_e,flux_w,flux_n,flux_s,u_curr,v_curr,u_pred,v_pred,...
    normCoff_e,normCoff_w,normCoff_n,normCoff_s,crossCoff_e,crossCoff_w,crossCoff_n,crossCoff_s,convScheme)

iter = 0 ;
residual_error = 1 ;

while ( residual_error > 1e-07 )
    uPred_residual = 0;
    vPred_residual = 0;
    
    iter = iter + 1 ;
    
    for  j = 2 : Ny
        for  i = 2 : Nx
            % upwind function
            [ue,uw,un,us,ve,vw,vn,vs,convecion_CentralCoeff] = convectiveScheme(convScheme,i,j,Nx,Ny,flux_e,flux_w,flux_n,flux_s,u_pred,v_pred) ;
            
            % total convection for each cell
            uConvection = flux_e(i,j) * ue + flux_w(i,j) * uw + flux_n(i,j) * un + flux_s(i,j) * us ;
            vConvection = flux_e(i,j) * ve + flux_w(i,j) * vw + flux_n(i,j) * vn + flux_s(i,j) * vs ;
            
            % volume interpolation for corner points for u
            [u_pred_ne,u_pred_nw,u_pred_se,u_pred_sw] = cornerPointsVolumeInterpolation(i,j,Vp,u_pred) ;
            % volume interpolation for corner points for v
            [v_pred_ne,v_pred_nw,v_pred_se,v_pred_sw] = cornerPointsVolumeInterpolation(i,j,Vp,v_pred) ;
            
            % calculating diffusion terms for each face for u predicted velocity
            [udiff_e,udiff_w,udiff_n,udiff_s] = diffusion(i,j,u_pred,u_pred_ne,u_pred_nw,u_pred_se,u_pred_sw,....
                normCoff_e,normCoff_w,normCoff_n,normCoff_s,crossCoff_e,crossCoff_w,crossCoff_n,crossCoff_s) ;
            % central coefficient for u predicted velocity diffusion term
            uPredDiff_CentralCoeff = normCoff_e(i,j) + normCoff_w(i,j) + normCoff_n(i,j) + normCoff_s(i,j) ;
            % total diffsion for u predicted velocity
            udiff  =  udiff_e + udiff_w + udiff_n + udiff_s  ;
            
            % calculating diffusion terms for each face for v predicted velocity
            [vdiff_e,vdiff_w,vdiff_n,vdiff_s] = diffusion(i,j,v_pred,v_pred_ne,v_pred_nw,v_pred_se,v_pred_sw,....
                normCoff_e,normCoff_w,normCoff_n,normCoff_s,crossCoff_e,crossCoff_w,crossCoff_n,crossCoff_s) ;
            % central coefficient for v predicted velocity diffusion term
            vPredDiff_CentralCoeff = normCoff_e(i,j) + normCoff_w(i,j) + normCoff_n(i,j) + normCoff_s(i,j) ;
            % total diffsion for v predicted velocity
            vdiff  =  vdiff_e + vdiff_w + vdiff_n + vdiff_s  ;
            
            % total central coefficient for u and v predicted velocities
            uPredTotal_CentralCoeff	= 1 + dt/Vp(i,j) * convecion_CentralCoeff + dt/Re/Vp(i,j) * uPredDiff_CentralCoeff  ;
            vPredTotal_CentralCoeff = 1 + dt/Vp(i,j) * convecion_CentralCoeff + dt/Re/Vp(i,j) * vPredDiff_CentralCoeff ;
            
            % calculating error/residue and correcting the value using Gauss Seidal
            uPred_error   	= u_curr(i,j) - u_pred(i,j) - dt/Vp(i,j) * uConvection + dt/Re * udiff/Vp(i,j) ;                              % residue for predicetd X velocity
            uPred_residual 	= uPred_residual + uPred_error * uPred_error;
            u_pred(i,j)  	= uPred_error/uPredTotal_CentralCoeff + u_pred(i,j) ;
            
            vPred_error  	= v_curr(i,j) - v_pred(i,j) - dt/Vp(i,j) * vConvection + dt/Re * vdiff/Vp(i,j) ;                              % residue for predicetd Y velocity
            vPred_residual	= vPred_residual + vPred_error * vPred_error;
            v_pred(i,j)    	= vPred_error/vPredTotal_CentralCoeff + v_pred(i,j) ;
            
        end
    end
    residual_error = sqrt( ( uPred_residual + vPred_residual )/(Nx*Ny)) ;                                                       % RMS for predicted velocity
end

fprintf('GaussSeidalSolver:  Solving for Ux predicted, Final residual = %e, No Iterations %d\n',sqrt(uPred_residual/Nx/Ny),iter) ;
fprintf('GaussSeidalSolver:  Solving for Uy predicted, Final residual = %e, No Iterations %d\n',sqrt(vPred_residual/Nx/Ny),iter) ;


end