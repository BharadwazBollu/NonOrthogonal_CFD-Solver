function [flux_e,flux_w,flux_n,flux_s] = flux(Nx,Ny,Vp,u_pred,v_pred,surface_e_x,surface_w_x,surface_n_x,surface_s_x,surface_e_y,surface_w_y,surface_n_y,surface_s_y,dt,diff_e,diff_w,diff_n,diff_s)

flux_e = zeros(Nx+1,Ny+1) ;     flux_w = zeros(Nx+1,Ny+1) ;
flux_n = zeros(Nx+1,Ny+1) ;     flux_s = zeros(Nx+1,Ny+1) ;

for  j = 2 : Ny
    for  i = 2 : Nx
        flux_e(i,j) =   ( Vp(i+1,j)*u_pred(i,j) + Vp(i,j)*u_pred(i+1,j) )/(Vp(i,j)+Vp(i+1,j)) * surface_e_x(i,j)  +  ( Vp(i+1,j)*v_pred(i,j) + Vp(i,j)*v_pred(i+1,j) )/(Vp(i,j)+Vp(i+1,j)) * surface_e_y(i,j)  - dt * diff_e(i,j) ;              % Calculating volume flux on all faces
        flux_w(i,j) =   ( Vp(i-1,j)*u_pred(i,j) + Vp(i,j)*u_pred(i-1,j) )/(Vp(i,j)+Vp(i-1,j)) * surface_w_x(i,j)  +  ( Vp(i-1,j)*v_pred(i,j) + Vp(i,j)*v_pred(i-1,j) )/(Vp(i,j)+Vp(i-1,j)) * surface_w_y(i,j)  - dt * diff_w(i,j) ;
        flux_n(i,j) =   ( Vp(i,j+1)*u_pred(i,j) + Vp(i,j)*u_pred(i,j+1) )/(Vp(i,j)+Vp(i,j+1)) * surface_n_x(i,j)  +  ( Vp(i,j+1)*v_pred(i,j) + Vp(i,j)*v_pred(i,j+1) )/(Vp(i,j)+Vp(i,j+1)) * surface_n_y(i,j)  - dt * diff_n(i,j) ;
        flux_s(i,j) =   ( Vp(i,j-1)*u_pred(i,j) + Vp(i,j)*u_pred(i,j-1) )/(Vp(i,j)+Vp(i,j-1)) * surface_s_x(i,j)  +  ( Vp(i,j-1)*v_pred(i,j) + Vp(i,j)*v_pred(i,j-1) )/(Vp(i,j)+Vp(i,j-1)) * surface_s_y(i,j)  - dt * diff_s(i,j) ;
    end
end

end