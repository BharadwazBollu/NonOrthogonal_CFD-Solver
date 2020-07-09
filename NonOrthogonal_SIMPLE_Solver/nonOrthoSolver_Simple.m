% *Solver for Non-Orthogonal/Non-Rectangular Mesh using Non-Dimensional Governing Equations
% with SIMPLE Algorithm and predictor corrector method
% Volumetric interpolation used for better accuracy
% Incluces functions for genrating Mesh, PreProcess and also postProcess
% Solver is powerful enough to handle Max Non_Orthogonal angle > 80 degrees
% Can handle highly skewed meshes also with high accuracy*
%%
clear
clc
close all
%%  Initialize Grid size and other parameters
Nx = 41 ;               % number of points in x-direction(should be odd number)
Ny = Nx ;               % number of points in y-direction(should be odd number)
Re = 100 ;              % Reynolds number
dt = 0.01 ;             % time step
export2tecplot = 0 ;  	% 1 for exporting and other number for ignoring
%%  Generation of Mesh/Grid
mesh_num = 1 ;
fprintf('  Creating Mesh \n\n  ');
if (mesh_num == 1)
    [X,Y] = lidDrivenMesh(Nx,Ny) ;          
elseif(mesh_num == 2)
    [X,Y] = pipeBend(Nx,Ny) ; 
elseif(mesh_num == 3)
    [X,Y] = bluffCircleMesh(Nx,Ny) ;
end


%%  PreProcess of Mesh
fprintf('  prePrcessing the Mesh \n\n  ');
preProcess                          % script for preprocessing
%%

u_curr = zeros(Nx+1,Ny+1) ; u_pred = zeros(Nx+1,Ny+1) ; u_Next = zeros(Nx+1,Ny+1) ;
v_curr = zeros(Nx+1,Ny+1) ; v_pred = zeros(Nx+1,Ny+1) ; v_Next = zeros(Nx+1,Ny+1) ;
p_Next = zeros(Nx+1,Ny+1) ;

flux_e = zeros(Nx+1,Ny+1) ; flux_w = zeros(Nx+1,Ny+1) ;
flux_n = zeros(Nx+1,Ny+1) ; flux_s = zeros(Nx+1,Ny+1) ;

pdiff_e = zeros(Nx+1,Ny+1) ; pdiff_w = zeros(Nx+1,Ny+1) ;
pdiff_n = zeros(Nx+1,Ny+1) ; pdiff_s = zeros(Nx+1,Ny+1) ;

iter    = 0;
error   = 1;
convScheme = 2 ; % 1 for upwind and 2 for QUICK

fprintf('  Starting time loop \n\n  ');

while error > 1e-3                                                  % Main error loop for Checking steady stat
    
    iter = iter + 1 ;                                                             % calculating number of iterations
    fprintf(' Time = %f \n', iter * dt );
    %% Predictor
    
    [flux_e,flux_w,flux_n,flux_s] = flux(Nx,Ny,Vp,u_pred,v_pred,surface_e_x,surface_w_x,surface_n_x,surface_s_x,surface_e_y,surface_w_y,surface_n_y,surface_s_y,dt,pdiff_e,pdiff_w,pdiff_n,pdiff_s) ;
    
    
    [u_pred,v_pred] = solvePredictedVelocty(Nx,Ny,Vp,Re,dt,flux_e,flux_w,flux_n,flux_s,u_curr,v_curr,u_pred,v_pred,...
    normCoff_e,normCoff_w,normCoff_n,normCoff_s,crossCoff_e,crossCoff_w,crossCoff_n,crossCoff_s,convScheme) ;
    
    % BC for Predicted velocities
    
    [u_pred,v_pred] = predictedVelocityBC(u_pred,v_pred,Nx,Ny,mesh_num) ;
    
    %% Pressure
    
    [p_Next,pdiff_e,pdiff_w,pdiff_n,pdiff_s] = solvePressure(Nx,Ny,Vp,dt,u_pred,v_pred,p_Next,surface_e_x,surface_w_x,surface_n_x,surface_s_x,...
    surface_e_y,surface_w_y,surface_n_y,surface_s_y,normCoff_e,normCoff_w,normCoff_n,normCoff_s,crossCoff_e,crossCoff_w,crossCoff_n,crossCoff_s,...
    pdiff_e,pdiff_w,pdiff_n,pdiff_s) ;
    % BC for p_Next
    
    [p_Next] = pressureBC(p_Next,Nx,Ny,mesh_num) ;
    
    %% Corerctor 
    
    [u_Next,v_Next,flux_e,flux_w,flux_n,flux_s] = solveCorrectedVelocty(Nx,Ny,Vp,Re,dt,u_curr,v_curr,u_pred,v_pred,p_Next,u_Next,v_Next,surface_e_x,surface_w_x,surface_n_x,surface_s_x,...
    surface_e_y,surface_w_y,surface_n_y,surface_s_y,normCoff_e,normCoff_w,normCoff_n,normCoff_s,crossCoff_e,crossCoff_w,crossCoff_n,crossCoff_s,...
    pdiff_e,pdiff_w,pdiff_n,pdiff_s,convScheme) ;
    
    % BC for Corrected velocities
    
    [u_Next,v_Next] = correctedVelocityBC(u_Next,v_Next,Nx,Ny,mesh_num) ;
    
    %% calculating the error values
    
    [error,u_curr,v_curr] = rmsError(u_curr,v_curr,u_Next,v_Next,Nx,Ny) ;
    
    error = error/dt;                                               % checking error convergence for steady state
    if ( error < 1e-3)
        break
    end
    
end

%% postProcessing 
% function for postProcessing of contours and also exporting fields to
% tecplot/ParaView with file extension of .plt
postProcess(x,y,X,Y,Nx,Ny,Vp,u_curr,v_curr,export2tecplot,mesh_num)  ;  
