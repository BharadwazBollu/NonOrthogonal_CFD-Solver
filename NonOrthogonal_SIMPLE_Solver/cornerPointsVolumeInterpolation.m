function [phi_ne,phi_nw,phi_se,phi_sw] = cornerPointsVolumeInterpolation(i,j,Vp,phi)

phi_ne = ( Vp(i,j)*phi(i+1,j+1) + Vp(i+1,j+1)*phi(i,j) + Vp(i,j+1)*phi(i+1,j) + Vp(i+1,j)*phi(i,j+1)  )/( Vp(i,j) + Vp(i+1,j) + Vp(i,j+1) + Vp(i+1,j+1) )   ;
phi_nw = ( Vp(i,j)*phi(i-1,j+1) + Vp(i-1,j+1)*phi(i,j) + Vp(i,j+1)*phi(i-1,j) + Vp(i-1,j)*phi(i,j+1)  )/( Vp(i,j) + Vp(i-1,j) + Vp(i,j+1) + Vp(i-1,j+1) )   ;
phi_se = ( Vp(i,j)*phi(i+1,j-1) + Vp(i+1,j-1)*phi(i,j) + Vp(i,j-1)*phi(i+1,j) + Vp(i+1,j)*phi(i,j-1)  )/( Vp(i,j) + Vp(i+1,j) + Vp(i,j-1) + Vp(i+1,j-1) )   ;
phi_sw = ( Vp(i,j)*phi(i-1,j-1) + Vp(i-1,j-1)*phi(i,j) + Vp(i,j-1)*phi(i-1,j) + Vp(i-1,j)*phi(i,j-1)  )/( Vp(i,j) + Vp(i-1,j) + Vp(i,j-1) + Vp(i-1,j-1) )   ;


end