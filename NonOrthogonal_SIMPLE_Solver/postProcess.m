function [] = postProcess(x,y,X,Y,Nx,Ny,Vp,u_curr,v_curr,export2tecplot,mesh_num)

x_cord  = zeros(Nx,1)  ;
y_cord  = zeros(1,Ny)  ;
u       = zeros(Nx,Ny) ;
v       = zeros(Nx,Ny) ;
%%

if ( mesh_num == 1 )
ghiaU=dlmread('ghiaU.dat');
choiU=dlmread('choiU.dat');
ghiaV=dlmread('ghiaV.dat');
choiV=dlmread('choiV.dat');

U = zeros(1,Ny+1) ;     V = zeros(1,Ny+1) ;

U(1) = 0.5 * ( u_curr((Nx+1)/2,1) + u_curr(ceil((Nx)/2),1)  );
x_cord(1) = 0 ;
for i = 2 : Ny 
    U(i) = ( u_curr((Nx+1)/2,i)*Vp((Nx+1)/2,i) + u_curr(ceil((Nx)/2),i)*Vp(ceil((Nx)/2),i) )/( Vp((Nx+1)/2,i) + Vp(ceil((Nx)/2),i) ) ;
   y_cord(i) = ( y((Nx+1)/2,i)*Vp((Nx+1)/2,i) + y(ceil((Nx)/2),i)*Vp(ceil((Nx)/2),i) )/( Vp((Nx+1)/2,i) + Vp(ceil((Nx)/2),i) )  ;
end
U(Ny+1) = 0.5 * ( ( u_curr((Nx+1)/2,Ny+1) + u_curr(ceil((Nx)/2),Ny+1) )/2  + ( u_curr((Nx+1)/2,Ny) + u_curr(ceil((Nx)/2),Ny) )/2 );
x_cord(Ny+1) = 1 ;

V(1) = 0.5 * ( ( v_curr(2,(Ny+1)/2)*Vp(2,(Ny+1)/2) + v_curr(2,ceil((Ny)/2))*Vp(2,ceil((Ny)/2)) )/( Vp(2,(Ny+1)/2) + Vp(2,ceil((Ny)/2)) ) + ( v_curr(1,(Ny+1)/2)*Vp(1,(Ny+1)/2) + v_curr(1,ceil((Ny)/2))*Vp(1,ceil((Ny)/2)) )/( Vp(1,(Ny+1)/2) + Vp(1,ceil((Ny)/2)) ) );
y_cord(1) = 0 ;
for i = 2 : Nx 
    V(i) = ( v_curr(i,(Ny+1)/2)*Vp(i,(Ny+1)/2) + v_curr(i,ceil((Ny)/2))*Vp(i,ceil((Ny)/2)) )/( Vp(i,(Ny+1)/2) + Vp(i,ceil((Ny)/2)) ) ;
    x_cord(i) = ( x(i,(Ny+1)/2)*Vp(i,(Ny+1)/2) + x(i,ceil((Ny)/2))*Vp(i,ceil((Ny)/2)) )/( Vp(i,(Ny+1)/2) + Vp(i,ceil((Ny)/2)) ) ;
end
V(Nx+1) = 0.5 * ( ( v_curr(Nx+1,(Ny+1)/2)*Vp(Nx+1,(Ny+1)/2) + v_curr(Nx+1,ceil((Ny)/2))*Vp(Nx+1,ceil((Ny)/2)) )/( Vp(Nx+1,(Ny+1)/2) + Vp(Nx+1,ceil((Ny)/2)) ) + ( v_curr(Nx,(Ny+1)/2)*Vp(Nx,(Ny+1)/2) + v_curr(Nx,ceil((Ny)/2))*Vp(Nx,ceil((Ny)/2)) )/( Vp(Nx,(Ny+1)/2) + Vp(Nx,ceil((Ny)/2)) ) );
y_cord(Ny+1) = 1 ;


figure(1)
plot(U,y_cord,'r*',ghiaU(:,1),ghiaU(:,2),'kd',choiU(:,1),choiU(:,2),'ko')
xlabel(' U Velocity ');
ylabel(' Y - distance ');
title('U Velocity Center line')
legend('UPWIND','GHIA','CHOI','Location','northwest')
% print(gcf,'U_UPWIND_IMPLICIT.jpg','-dpng','-r300');

figure(2)
plot(x_cord,V,'r*',ghiaV(:,1),ghiaV(:,2),'kd',choiV(:,1),choiV(:,2),'ko')
xlim([ 0 1 ] );
grid on
xlabel(' X - distance ');
ylabel(' V Velocity Center line');
title('V Velocity')
legend('UPWIND','GHIA','CHOI','Location','northeast')
% print(gcf,'V_UPWIND_IMPLICIT.jpg','-dpng','-r300');

end

% interior points
for j = 2:Ny-1
    for i = 2:Nx-1
        u(i,j) = ( Vp(i,j)*u_curr(i+1,j+1) + Vp(i+1,j+1)*u_curr(i,j) + Vp(i,j+1)*u_curr(i+1,j) + Vp(i+1,j)*u_curr(i,j+1)  )/( Vp(i,j) + Vp(i+1,j) + Vp(i,j+1) + Vp(i+1,j+1) )   ;
        v(i,j) = ( Vp(i,j)*v_curr(i+1,j+1) + Vp(i+1,j+1)*v_curr(i,j) + Vp(i,j+1)*v_curr(i+1,j) + Vp(i+1,j)*v_curr(i,j+1)  )/( Vp(i,j) + Vp(i+1,j) + Vp(i,j+1) + Vp(i+1,j+1) )   ;
    end
end

for j = 1 : Ny
    u(1,j)   =  0.5*( u_curr(1,j+1)    + u_curr(1,j)    )  ;            % left 
    v(1,j)   =  0.5*( v_curr(1,j+1)    + v_curr(1,j)    )  ;
    u(Nx,j)  =  0.5*( u_curr(Nx+1,j+1) + u_curr(Nx+1,j) )  ;            % right
    v(Nx,j)  =  0.5*( v_curr(Nx+1,j+1) + v_curr(Nx+1,j) )  ;
end

for i = 1:Nx
    u(i,Ny)  =  0.5*( u_curr(i+1,Ny+1) + u_curr(i,Ny+1) )  ;            % Top
    v(i,Ny)  =  0.5*( v_curr(i+1,Ny+1) + v_curr(i,Ny+1) )  ;
    u(i,1)   =  0.5*( u_curr(i+1,1)    + u_curr(i,1) )     ;            % Bottom
    v(i,1)   =  0.5*( v_curr(i+1,1)    + v_curr(i,1) )     ;
end

figure(3)
contourf(X,Y,u,25)
xlabel(' X ');
ylabel(' Y ');
title('U velocity')
colorbar
%print(gcf,'Ucontour.jpg','-dpng','-r300');

figure(4)
contourf(X,Y,v,25)
xlabel(' X ');
ylabel(' Y ');
title('V velocity')
colorbar
%print(gcf,'Vcontour.jpg','-dpng','-r300');

figure(5)
quiver(x,y,u_curr,v_curr)
xlabel(' X ');
ylabel(' Y ');
title(' Quiver ')
% xlim([0 1])
% ylim([-0.35 0.35])
%print(gcf,'Quiver.jpg','-dpng','-r300');

%%  Exporting to Tecplot/ParaView (BLUFF.plt):
if ( export2tecplot == 1 )
    fid=fopen('BLUFF.plt','w');             % file name
    fprintf(fid,'VARIABLES= X,Y,U,V\n');    % parameters to export
    fprintf(fid,'ZONE I=');
    fprintf(fid,'%d',Nx);
    fprintf(fid,',J=');
    fprintf(fid,'%d',Ny);
    fprintf(fid,' F=POINT\n');
    for j=1:Ny
        for i=1:Nx
            fprintf(fid,'%f %f %f %f \n',x(i,j),y(i,j),u_curr(i,j),v_curr(i,j));
        end
    end
    fclose(fid);
end

end