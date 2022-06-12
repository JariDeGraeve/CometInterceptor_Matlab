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

% [x,y,z] = meshgrid(-10:0.1:10);
% data = x.^2 + y.^2 + z.^2;
% cdata = smooth3(rand(size(data)),'box',7);
% p = patch(isosurface(x,y,z,data));
% isonormals(x,y,z,data,p)
% isocolors(x,y,z,cdata,p)
% p.FaceColor = 'interp';
% p.EdgeColor = 'none';
% view(150,30)
% daspect([1 1 1])
% axis tight
% camlight
% lighting gouraud
%% Problem 1: Q/r^gamma
function [ gamma, dgamma, Q, dQ ] = FindSolution_QOverRgamma( debug, Q0, gamma0, t, r, phi, theta, f, df, tol, v_gas )
    nr_of_sets = length(t);
    n = zeros( nr_of_sets, 1 );
    for k = 1:nr_of_sets
        n(k) = length(t{k});
    end
    data = ...
        struct( ...
            'debug', debug, ...
            'nr_of_sets', length(t), ...
            'n', n, ...
            't', {t}, ...
            'r', {r}, ...
            'phi', {phi}, ...
            'theta', {theta}, ...
            'f', {f}, ...
            'df', {df}, ...
            'v_gas', v_gas ...
        );
    x0 = [ gamma0; log(Q0) ];   % log is natural logarithm (ln(Q0))
    sx = [ 0.1; 1 ];
    [ x, dx ] = ...
        expl_optimization( [], 'FindSolution_QOverRgamma', x0, sx, tol, @TargetFunction_QOverRgamma, data, 'enhanced', 'silent' );
    gamma = x(1);
    Q = exp(x(2));
    dgamma = dx(1);
    dQ = Q * dx(2);
end
    
function F = TargetFunction_QOverRgamma( x, data )
    %% Target function
    % Least-squares formulation given measurements at points
    gamma = x(1);
%     fprintf("%f\n",gamma )
    if ( gamma < 0 ) 
        F = 1e20 * (-gamma);
    elseif ( gamma > 10 )
        F = 1e20 * (gamma-10);
    else
        Q = exp(x(2));
        F = 0;
        f_model = cell( size(data.f) );
        for k = 1:data.nr_of_sets
            f_model{k} = CometModel_QOverRgamma( data.r{k}, data.phi{k}, data.theta{k}, Q, gamma, data.v_gas );
            F = F + ...
                sum( ...
                    ( ( data.f{k} - f_model{k} ) ./ data.df{k} ).^2 ...
                )/data.n(k);
        end
        F = F / data.nr_of_sets;
        if data.debug
            DebugPlot( data.nr_of_sets, data.t, data.f, f_model, sprintf( 'gamma %.4f Q %.4e', gamma, Q ) );
        end
    end
end


if algo.do_QOverRgamma_sza
        t_ls = ( -t0 : dt_ls : t0 )';
        n_t_ls = length(t_ls);
        gamma_13 = zeros( n_t_ls, 1 );
        dgamma_13 = zeros( n_t_ls, 1 );
        Q_13 = zeros( n_t_ls, 1 );
        dQ_13 = zeros( n_t_ls, 1 );
        fprintf( 'Problem 1.3 (Q/r^gamma at given sza)\n' );
        fprintf( '\n' );
        fprintf( '      Q                             gamma\n' );
        sza_A = SolarZenithAngle(phi_A, theta_A);
        sza_B2 = SolarZenithAngle(phi_B2, theta_B2);
        sza_B1 = SolarZenithAngle(phi_B1, theta_B1);
        for i = 1:n_t_ls
            t_ref = t_ls(i);
            sza_ref = interp1( t, sza_A, t_ref );
            [ factor_A, selection_A ] = HomogeneitySelection( sza_A, sza_ref, dsza_homogeneity );           
            [ factor_B2, selection_B2 ] = HomogeneitySelection( sza_B2, sza_ref, dsza_homogeneity );           
            [ factor_B1, selection_B1 ] = HomogeneitySelection( sza_B1, sza_ref, dsza_homogeneity );           
            set_t = { t(selection_A); t(selection_B2); t(selection_B1) };
            set_r = { r_A(selection_A); r_B2(selection_B2); r_B1(selection_B1) };
            set_theta = { theta_A(selection_A); theta_B2(selection_B2); theta_B1(selection_B1) };
            set_phi = { phi_A(selection_A); phi_B2(selection_B2); phi_B1(selection_B1) };
            set_f = { f_A_obs(selection_A); f_B2_obs(selection_B2); f_B1_obs(selection_B1) };
            set_df = ...
                { ...
                    rel_error * f_A_obs(selection_A) .* factor_A(selection_A); ...
                    rel_error * f_B2_obs(selection_B2) .* factor_B2(selection_B2);  ...
                    rel_error * f_B1_obs(selection_B1) .* factor_B1(selection_B1)  ...
                };
            [ gamma_13(i), dgamma_13(i), Q_13(i), dQ_13(i) ] = ...
                FindSolution_QOverRgamma( ...
                    false, Q_start, gamma_start, ...
                    set_t, set_r, set_phi, set_theta, set_f, set_df, ...
                    tolerance, v_gas ...
                );
            fprintf( '    : %12.6e +/- %12.6e  %10.6f +/- %10.6f\n', Q_13(i), dQ_13(i), gamma_13(i), dgamma_13(i) );
        end
        fprintf('mean: %12.6e                   %10.6f\n', exp(mean(log(Q_13), 'omitnan')), mean(gamma_13, 'omitnan'));
        fprintf( '\n' );
    end

if algo.do_QOverRSquared_sza
        t_ls = ( -t0 : dt_ls : t0 )';
        n_t_ls = length(t_ls);
        Q_33 = zeros( n_t_ls, 1 );
        dQ_33 = zeros( n_t_ls, 1 );
        fprintf( 'Problem 3.3 (Q/r^2 at given sza)\n' );
        fprintf( '\n' );
        fprintf( '      Q\n' );
        sza_A = SolarZenithAngle(phi_A, theta_A);
        sza_B2 = SolarZenithAngle(phi_B2, theta_B2);
        sza_B1 = SolarZenithAngle(phi_B1, theta_B1);
        for i = 1:n_t_ls
            t_ref = t_ls(i);
            sza_ref = interp1( t, sza_A, t_ref );
            [ factor_A, selection_A ] = HomogeneitySelection( sza_A, sza_ref, dsza_homogeneity );           
            [ factor_B2, selection_B2 ] = HomogeneitySelection( sza_B2, sza_ref, dsza_homogeneity );           
            [ factor_B1, selection_B1 ] = HomogeneitySelection( sza_B1, sza_ref, dsza_homogeneity );           
            set_t = { t(selection_A); t(selection_B2); t(selection_B1) };
            set_r = { r_A(selection_A); r_B2(selection_B2); r_B1(selection_B1) };
            set_theta = { theta_A(selection_A); theta_B2(selection_B2); theta_B1(selection_B1) };
            set_phi = { phi_A(selection_A); phi_B2(selection_B2); phi_B1(selection_B1) };
            set_f = { f_A_obs(selection_A); f_B2_obs(selection_B2); f_B1_obs(selection_B1) };
            set_df = ...
                { ...
                    rel_error * f_A_obs(selection_A) .* factor_A(selection_A); ...
                    rel_error * f_B2_obs(selection_B2) .* factor_B2(selection_B2);  ...
                    rel_error * f_B1_obs(selection_B1) .* factor_B1(selection_B1)  ...
                };
            [ Q_33(i), dQ_33(i) ] = ...
                FindSolution_QOverRSquared( ...
                    false, Q_start, ...
                    set_t, set_r, set_phi, set_theta, set_f, set_df, ...
                    tolerance, v_gas ...
                );
            fprintf( '    : %12.6e +/- %12.6e\n', Q_33(i), dQ_33(i));
        end
        fprintf('mean: %12.6e\n', exp(mean(log(Q_33), 'omitnan')));
        fprintf( '\n' );
    end
