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


