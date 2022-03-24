a = 0.8;
% 
% U = linspace(0,0.1,2*pi);
% V = linspace(0,0.1,2*pi);
% [u, v] = meshgrid(U,V);
% 
% xt = 
% yt = 
% zt = 



funx = @(u,v) cos(u).*((1+a*cos(u))/(1+a));
funy = @(u,v) sin(u).*((1+a*cos(u))/(1+a)).*cos(v);
funz = @(u,v) sin(u).*((1+a*cos(u))/(1+a)).*sin(v); 
fsurf(funx,funy,funz, [0 2*pi 0 pi])
