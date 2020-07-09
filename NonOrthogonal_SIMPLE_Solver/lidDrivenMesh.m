function [X,Y] = lidDrivenMesh(nx,ny)
% creates 2D non-orthogonal mesh by generating structured mesh and giving some
% random +/- error for all inside points which makes it non-rectangualar
x=linspace(0,1,nx);
y=linspace(0,1,ny);

[Y,X] = meshgrid(x,y) ;

for j=2:ny-1
    for i=2:nx-1
        X(i,j) = X(i,j) + (-1)^(randi(9))*rand/(2.5*(nx-1)) ;    % random +/- error
        Y(i,j) = Y(i,j) + (-1)^(randi(9))*rand/(2.5*(ny-1)) ;
    end
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