function [X,Y] = pipeBend(nx,ny)

X=zeros(nx,ny);         Y=zeros(nx,ny);

temp = linspace(6,10,nx);

for j=1:ny
    X(:,j)= temp(j)*cos(linspace(1/2*pi,0,nx)) ;
    Y(:,j)= temp(j)*sin(linspace(1/2*pi,0,ny)) ;
end

%% uncomment below lines if you want to see grid points and mesh
% clf
% plot(X,Y,'k*')
% hold on
% axis equal
% for m=1:nx
% plot(X(m,:),Y(m,:),'b');
% end
% for m=1:ny
% plot(X(:,m),Y(:,m),'Color',[0 0 0]);
% end

end