% % x = [1 2 3];
% % y = [1 2 4];
% % z = [1 2 5];
% x = ones(3,3,3) + 1;
% y = ones(3,3,3) + 2;
% z = ones(3,3,3) + 3;
% 
% x1 = 3
% y1 = 2
% z1 = 1
% 
% 
% x
% a = [x; y; z]
% b = repmat([x1; y1; z1],1, length([x y z]), length(x))
% 
% dot(a, b)

% v1 =  [1 2 3] 
% [X,Y,Z] = meshgrid( linspace(-5,5,5), linspace(-5,5,5),  linspace(-5,5,5))
% x = X(:) ; y = Y(:) ; z = Z(:);
% v2 = [x y z];
% d = sum(v1.*v2,2) ;
% A = reshape(d,size(X))

% Q_gas = 5e10;
% r_comet = 2000;
% q_gas = @(phi, theta) (Q_gas ./ ( 4 * pi * r_comet.^2 )) .* sin(theta) .* r_comet.^2; 
% 
% 
% 
% integral2(q_gas , 0, 2*pi, 0, pi )


% theta = linspace(0, 2*pi);
% polarplot(theta, 1 - 0.8 *cos(theta + pi), "LineWidth", 3)
% hold on;
% polarplot(theta,0.02+zeros(size(theta)), "LineWidth", 6)

% [X,Y,Z] = meshgrid(-5:0.2:5);
% V = X.*exp(-X.^2-Y.^2-Z.^2);
% 
% [xsurf,ysurf] = meshgrid(-2:0.2:2);
% zsurf = -xsurf.^2-ysurf.^2;
% slice(X,Y,Z,V,xsurf,ysurf,zsurf)

% clc,clear
%     % generate some data
% theta1 = linspace(-1,1,60)*pi;
% phi1 = linspace(-1,1,20)*pi/2;
% r1 = 0:10;
% [theta, phi, r] = meshgrid(theta1,phi1,r1);
% f = r + cos(10*theta);
%     % boundaries of volume
% [X,Y,Z] = sphere(20);
% X = r1(end)*X;
% Y = r1(end)*Y;
% Z = r1(end)*Z;
%     % create isosurface where f=5
% p = patch(isosurface(theta,phi,r,f,5));
%     % extract spherical data
% fc = get(p,'Faces');
% vc = get(p,'Vertices');
%     % convert spherical data to cartesian
% [x1,y1,z1] = sph2cart(vc(:,1),vc(:,2),vc(:,3));
% vc = [x1 y1 z1];
%     % plot boundaries
% cla
% surf(X,Y,Z,'Facecolor','none','edgeColor',[1 1 1]*0.5)
% hold on
%     % plot data f=5 in cartesian
% h = patch('Faces',fc,'Vertices',vc);
% hold off
% set(h,'FaceColor','r','EdgeColor','none');
% camlight 
% lighting gouraud
% 

[x,y,z] = meshgrid(-10:0.1:10);
data = x.^2 + y.^2 + z.^2;
cdata = smooth3(rand(size(data)),'box',7);
p = patch(isosurface(x,y,z,data));
isonormals(x,y,z,data,p)
isocolors(x,y,z,cdata,p)
p.FaceColor = 'interp';
p.EdgeColor = 'none';
view(150,30)
daspect([1 1 1])
axis tight
camlight
lighting gouraud
