%% Script for PreProcessing of Mesh
Vp = zeros(Nx+1,Ny+1) ;     % Volume array
x  = zeros(Nx+1,Ny+1) ;     % Cell centers for x-cordinate
y  = zeros(Nx+1,Ny+1) ;     % Cell centers for y-cordinate

surface_e_x = zeros(Nx+1,Ny+1) ;    surface_e_y = zeros(Nx+1,Ny+1) ;    % surface_e_x gives x-cordinate of normal surface projected outward from cell center
surface_w_x = zeros(Nx+1,Ny+1) ;    surface_w_y = zeros(Nx+1,Ny+1) ;    % surface_w_y gives y-cordinate of normal surface projected outward from cell center
surface_n_x = zeros(Nx+1,Ny+1) ;    surface_n_y = zeros(Nx+1,Ny+1) ;
surface_s_x = zeros(Nx+1,Ny+1) ;    surface_s_y = zeros(Nx+1,Ny+1) ;
% Normal diffussion coefficients of all faces
normCoff_e = zeros(Nx+1,Ny+1) ;   normCoff_w = zeros(Nx+1,Ny+1) ;
normCoff_n = zeros(Nx+1,Ny+1) ;   normCoff_s = zeros(Nx+1,Ny+1) ;
% Cross diffussion coefficients of all faces
crossCoff_e = zeros(Nx+1,Ny+1) ;   crossCoff_w = zeros(Nx+1,Ny+1) ;
crossCoff_n = zeros(Nx+1,Ny+1) ;   crossCoff_s = zeros(Nx+1,Ny+1) ;
% Calculating non orthogoanlity angle for all faces
nonOrthoAngle_e = zeros(Nx+1,Ny+1) ;    nonOrthoAngle_n = zeros(Nx+1,Ny+1) ;
nonOrthoAngle_w = zeros(Nx+1,Ny+1) ;    nonOrthoAngle_s = zeros(Nx+1,Ny+1) ;

%% Creating cell corners and face centers which are used as boundaries
x(1,1) = X(1,1) ;               y(1,1) = Y(1,1) ;           % left bottom corner
x(1,Nx+1) = X(1,Nx) ;           y(1,Nx+1) = Y(1,Nx) ;       % left top corner
x(Nx+1,Nx+1) = X(Nx,Nx) ;       y(Nx+1,Nx+1) = Y(Nx,Nx) ;   % right top corner
x(Nx+1,1) = X(Nx,1) ;           y(Nx+1,1) = Y(Nx,1) ;       % right bottom corner

for i=2:Ny
    x(i,1) = ( X(i,1) + X(i-1,1) )/2 ;          % Bottom boundary 
    y(i,1) = ( Y(i,1) + Y(i-1,1) )/2 ;
    
    x(i,Nx+1) = ( X(i,Nx) + X(i-1,Nx) )/2 ;     % Top boundary
    y(i,Nx+1) = ( Y(i,Nx) + Y(i-1,Nx) )/2 ;
end

for j=2:Nx
    x(1,j) = ( X(1,j) + X(1,j-1) )/2 ;          % Left boundary
    y(1,j) = ( Y(1,j) + Y(1,j-1) )/2 ;
    
    x(Nx+1,j) = ( X(Nx,j) + X(Nx,j-1) )/2 ;     % Right boundary
    y(Nx+1,j) = ( Y(Nx,j) + Y(Nx,j-1) )/2 ;
end

for j=2:Ny
    for i=2:Nx
        % taking corner points for each cell for both co-ordinates
        x_sw_vertex = X(i-1,j-1)     ;
        x_se_vertex = X(i,j-1)       ;
        x_nw_vertex = X(i-1,j)       ;
        x_ne_vertex = X(i,j)         ;
        
        y_sw_vertex = Y(i-1,j-1)     ;
        y_se_vertex = Y(i,j-1)       ;
        y_nw_vertex = Y(i-1,j)       ;
        y_ne_vertex = Y(i,j)         ;
        % calculating centroids/cell enters
        x(i,j) = ( x_sw_vertex + x_se_vertex + x_nw_vertex + x_ne_vertex )/4 ;
        y(i,j) = ( y_sw_vertex + y_se_vertex + y_nw_vertex + y_ne_vertex )/4 ;
        % calculating surface co-ordinates for all faces
        surface_e_x(i,j) = y_ne_vertex - y_se_vertex 	;
        surface_w_x(i,j) = y_sw_vertex - y_nw_vertex 	;
        surface_n_x(i,j) = y_nw_vertex - y_ne_vertex  	;
        surface_s_x(i,j) = y_se_vertex - y_sw_vertex  	;
        
        surface_e_y(i,j) = x_se_vertex - x_ne_vertex  	;
        surface_w_y(i,j) = x_nw_vertex - x_sw_vertex   	;
        surface_n_y(i,j) = x_ne_vertex - x_nw_vertex  	;
        surface_s_y(i,j) = x_sw_vertex - x_se_vertex  	;
        
        % x-component of vector joining cell center and face center
        dpe_x = ( x_se_vertex + x_ne_vertex )/2 - x(i,j)  ;
        dpw_x = ( x_nw_vertex + x_sw_vertex )/2 - x(i,j)  ;
        dpn_x = ( x_ne_vertex + x_nw_vertex )/2 - x(i,j)  ;
        dps_x = ( x_sw_vertex + x_se_vertex )/2 - x(i,j)  ;
        % y-component of vector joining cell center and face center  
        dpe_y = ( y_ne_vertex + y_se_vertex )/2 - y(i,j)  ;
        dpw_y = ( y_sw_vertex + y_nw_vertex )/2 - y(i,j)  ;
        dpn_y = ( y_nw_vertex + y_ne_vertex )/2 - y(i,j)  ;
        dps_y = ( y_se_vertex + y_sw_vertex )/2 - y(i,j)  ;
        
        % Performing dot product of surface normal vector and vector joining cell
        % center and face center for all faces
        dot_e = surface_e_x(i,j) * dpe_x + surface_e_y(i,j) * dpe_y  ;
        dot_w = surface_w_x(i,j) * dpw_x + surface_w_y(i,j) * dpw_y  ;
        dot_n = surface_n_x(i,j) * dpn_x + surface_n_y(i,j) * dpn_y  ;
        dot_s = surface_s_x(i,j) * dps_x + surface_s_y(i,j) * dps_y  ;
        % Making sure that computed surafce vectors are pointing outward from cell centers
        if ( dot_e < 0 )
            surface_e_x(i,j) = -surface_e_x(i,j) ;
            surface_e_y(i,j) = -surface_e_y(i,j) ;
        end
        if ( dot_w < 0 )
            surface_w_x(i,j) = -surface_w_x(i,j) ;
            surface_w_y(i,j) = -surface_w_y(i,j) ;
        end
        if ( dot_n < 0 )
            surface_n_x(i,j) = -surface_n_x(i,j) ;
            surface_n_y(i,j) = -surface_n_y(i,j) ;
        end
        if ( dot_s < 0 )
            surface_s_x(i,j) = -surface_s_x(i,j) ;
            surface_s_y(i,j) = -surface_s_y(i,j) ;
        end
        % Computing volume for each cell
        Vp(i,j) = abs(0.5*( x_sw_vertex*y_se_vertex - x_se_vertex*y_sw_vertex  +  x_se_vertex*y_ne_vertex - x_ne_vertex*y_se_vertex...
            +  x_ne_vertex*y_nw_vertex -  x_nw_vertex*y_ne_vertex  +  x_nw_vertex*y_sw_vertex - x_sw_vertex*y_nw_vertex  ))  ;
        
    end
end

for j=2:Ny
    for i=2:Nx
        % calculating distance between cell center and other faces i,e delta value
        dpe = sqrt(( x(i+1,j) - x(i,j) )^2 + ( y(i+1,j) - y(i,j) )^2) ;
        dpw = sqrt(( x(i-1,j) - x(i,j) )^2 + ( y(i-1,j) - y(i,j) )^2) ;
        dpn = sqrt(( x(i,j+1) - x(i,j) )^2 + ( y(i,j+1) - y(i,j) )^2) ;
        dps = sqrt(( x(i,j-1) - x(i,j) )^2 + ( y(i,j-1) - y(i,j) )^2) ;
        % calculating surface magnitude for all faces
        spe = sqrt(( surface_e_x(i,j) )^2 + ( surface_e_y(i,j) )^2) ;
        spw = sqrt(( surface_w_x(i,j) )^2 + ( surface_w_y(i,j) )^2) ;
        spn = sqrt(( surface_n_x(i,j) )^2 + ( surface_n_y(i,j) )^2) ;
        sps = sqrt(( surface_s_x(i,j) )^2 + ( surface_s_y(i,j) )^2) ;
        % x-component of unit vector along line joing cell centers
        e1_x = (x(i+1,j)-x(i,j))/( (x(i+1,j)-x(i,j))^2 + (y(i+1,j)-y(i,j))^2  )^0.5 ;
        w1_x = (x(i-1,j)-x(i,j))/( (x(i-1,j)-x(i,j))^2 + (y(i-1,j)-y(i,j))^2  )^0.5  ;
        n1_x = (x(i,j+1)-x(i,j))/( (x(i,j+1)-x(i,j))^2 + (y(i,j+1)-y(i,j))^2  )^0.5  ;
        s1_x = (x(i,j-1)-x(i,j))/( (x(i,j-1)-x(i,j))^2 + (y(i,j-1)-y(i,j))^2  )^0.5  ;
        % y-component of unit vector along line joing cell centers
        e1_y = (y(i+1,j)-y(i,j))/( (x(i+1,j)-x(i,j))^2 + (y(i+1,j)-y(i,j))^2  )^0.5  ;
        w1_y = (y(i-1,j)-y(i,j))/( (x(i-1,j)-x(i,j))^2 + (y(i-1,j)-y(i,j))^2  )^0.5  ;
        n1_y = (y(i,j+1)-y(i,j))/( (x(i,j+1)-x(i,j))^2 + (y(i,j+1)-y(i,j))^2  )^0.5  ;
        s1_y = (y(i,j-1)-y(i,j))/( (x(i,j-1)-x(i,j))^2 + (y(i,j-1)-y(i,j))^2  )^0.5  ;
        
        % x-component of Normal diffussion 
        Ee1_x = ( (surface_e_x(i,j))^2 + (surface_e_y(i,j))^2 )* e1_x/( surface_e_x(i,j)*e1_x + surface_e_y(i,j)*e1_y) ;
        Ew1_x = ( (surface_w_x(i,j))^2 + (surface_w_y(i,j))^2 )* w1_x/( surface_w_x(i,j)*w1_x + surface_w_y(i,j)*w1_y) ;
        En1_x = ( (surface_n_x(i,j))^2 + (surface_n_y(i,j))^2 )* n1_x/( surface_n_x(i,j)*n1_x + surface_n_y(i,j)*n1_y) ;
        Es1_x = ( (surface_s_x(i,j))^2 + (surface_s_y(i,j))^2 )* s1_x/( surface_s_x(i,j)*s1_x + surface_s_y(i,j)*s1_y) ;
        % y-component of Normal diffussion
        Ee1_y = ( (surface_e_x(i,j))^2 + (surface_e_y(i,j))^2 )* e1_y/( surface_e_x(i,j)*e1_x + surface_e_y(i,j)*e1_y) ;
        Ew1_y = ( (surface_w_x(i,j))^2 + (surface_w_y(i,j))^2 )* w1_y/( surface_w_x(i,j)*w1_x + surface_w_y(i,j)*w1_y) ;
        En1_y = ( (surface_n_x(i,j))^2 + (surface_n_y(i,j))^2 )* n1_y/( surface_n_x(i,j)*n1_x + surface_n_y(i,j)*n1_y) ;
        Es1_y = ( (surface_s_x(i,j))^2 + (surface_s_y(i,j))^2 )* s1_y/( surface_s_x(i,j)*s1_x + surface_s_y(i,j)*s1_y) ;
        % Normal diffussion coefficient 
        normCoff_e(i,j) = (sqrt( (Ee1_x)^2 + (Ee1_y)^2 ))/dpe ;
        normCoff_w(i,j) = (sqrt( (Ew1_x)^2 + (Ew1_y)^2 ))/dpw ;
        normCoff_n(i,j) = (sqrt( (En1_x)^2 + (En1_y)^2 ))/dpn ;
        normCoff_s(i,j) = (sqrt( (Es1_x)^2 + (Es1_y)^2 ))/dps ;
        
        % x-component of Cross diffussion
        Te1_x = surface_e_x(i,j) - Ee1_x ;
        Tw1_x = surface_w_x(i,j) - Ew1_x ;
        Tn1_x = surface_n_x(i,j) - En1_x ;
        Ts1_x = surface_s_x(i,j) - Es1_x ;
        % y-component of Cross diffussion
        Te1_y = surface_e_y(i,j) - Ee1_y ;
        Tw1_y = surface_w_y(i,j) - Ew1_y ;
        Tn1_y = surface_n_y(i,j) - En1_y ;
        Ts1_y = surface_s_y(i,j) - Es1_y ;
        % Cross diffussion coefficient
        crossCoff_e(i,j) = (sqrt( (Te1_x)^2 + (Te1_y)^2 ))/spe ;
        crossCoff_w(i,j) = (sqrt( (Tw1_x)^2 + (Tw1_y)^2 ))/spw ;
        crossCoff_n(i,j) = (sqrt( (Tn1_x)^2 + (Tn1_y)^2 ))/spn ;
        crossCoff_s(i,j) = (sqrt( (Ts1_x)^2 + (Ts1_y)^2 ))/sps ;
        % Calculating non orthonality angles
        nonOrthoAngle_e(i,j) = acosd( ( e1_x*surface_e_x(i,j) + e1_y*surface_e_y(i,j) )/spe ) ;
        nonOrthoAngle_w(i,j) = acosd( ( w1_x*surface_w_x(i,j) + w1_y*surface_w_y(i,j) )/spw ) ;
        nonOrthoAngle_n(i,j) = acosd( ( n1_x*surface_n_x(i,j) + n1_y*surface_n_y(i,j) )/spn ) ;
        nonOrthoAngle_s(i,j) = acosd( ( s1_x*surface_s_x(i,j) + s1_y*surface_s_y(i,j) )/sps ) ;    
            
    end
end

max_nonOrthoAngle_e = max(max(nonOrthoAngle_e)) ;
max_nonOrthoAngle_w = max(max(nonOrthoAngle_w)) ;
max_nonOrthoAngle_n = max(max(nonOrthoAngle_n)) ;
max_nonOrthoAngle_s = max(max(nonOrthoAngle_s)) ;

overAll_max_nonOrthoAngle = max( [ max_nonOrthoAngle_e  max_nonOrthoAngle_w  max_nonOrthoAngle_n  max_nonOrthoAngle_s ] ) ;


% hold on
% plot(X,Y,'k*')
% plot(x,y,'r*')

% clf
% hold on
% axis equal
% for m=1:Nx
% plot(X(m,:),Y(m,:),'b');
% end
% for m=1:Ny
% plot(X(:,m),Y(:,m),'Color',[0 0 0]);
% end
% title('Mesh Grid')
% xlim([ 0 1 ] ); ylim([-0.3393 0.3393])
% plot(x,y,'r*')
% print(gcf,'Mesh_Re100_d0_25.jpg','-dpng','-r300');


