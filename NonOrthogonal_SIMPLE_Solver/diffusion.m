function [diff_e,diff_w,diff_n,diff_s] = diffusion(i,j,phi,phi_ne,phi_nw,phi_se,phi_sw,....
    normCoff_e,normCoff_w,normCoff_n,normCoff_s,crossCoff_e,crossCoff_w,crossCoff_n,crossCoff_s)

diff_e = normCoff_e(i,j)*( phi(i+1,j) - phi(i,j)) + crossCoff_e(i,j)*( phi_ne - phi_se  ) ;
diff_w = normCoff_w(i,j)*( phi(i-1,j) - phi(i,j)) + crossCoff_w(i,j)*( phi_sw - phi_nw  ) ;
diff_n = normCoff_n(i,j)*( phi(i,j+1) - phi(i,j)) + crossCoff_n(i,j)*( phi_nw - phi_ne  ) ;
diff_s = normCoff_s(i,j)*( phi(i,j-1) - phi(i,j)) + crossCoff_s(i,j)*( phi_se - phi_sw  ) ;


end