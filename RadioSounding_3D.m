function [] = RadioSounding_3D()

    % Algorithm selection
%     algo = ...
%         struct( ...
%             'do_OneOverRgamma_A', true ... 
%         );
    
    algo = ...
        struct( ...
              'do_OneOverRgamma_A', true ... 
        );
    RadioSounding_Core( algo );
    
end

function [] = RadioSounding_Core( algo )
    
    %% Parameters
    % Physical parameters
    r_comet = 2000;            % comet radius, in m
    Q_gas_kg = 1000;           % gas production rate, in kg/s
    m_ave = 18 * 1.6e-27;      % average neutral gas particle mass, in kg/particle
    v_gas = 1000;              % neutral gas speed, in m/s
    D_ca_A = 1000000;          % distance at closest approach of A, in m
    t_ca_A = 0;                % time at closest approach of A, in s

    delta_B1 = [-183000, 150000, -766000];   % coordinates spacecraft B1 relative to A, in m
    delta_B2 = [-233000, -350000, -466000];  % coordinates spacecraft B2 relative to A, in m

    V_flyby = [42000, 42000, 0];  % flyby speedvector, in m/s
    phi_ca = 30;               % azimuth angle at closest approach, in degrees
    theta_ca = 30;             % polar angle at closest approach, in degrees
    gamma = 0.8;               % ionization rate exponent, dimensionless
    alpha = 0.00002;           % ionization rate, in 1/s
    a = 0.8;                   % gas coma asymmetry factor, dimensionless
    theta_jet = [ 90,60];       % polar angle of center of gas jet, in degrees
    phi_jet = [50, 0];        % azimuth angle of center of gas jet, in degrees
    dangle_jet = [10, 20];     % angular half-width of gas jet, in degrees
    f_jet = [3,5];             % density contrast in jet, dimensionless
%    sc_inclination = 0;        % inclination of orbit
    rel_error = 0.05;          % relative measurement error, dimensionless
    t0 = 150;                  % limits of time interval, in s
    dt = 0.2;                  % measurement time resolution, in s
    
    % Numerical parameters
    tolerance = 1e-6;          % tolerance for solving optimization problem

    % Plot parameters
    font_size = 14;
    marker_size = 8;
    line_width = 2;
    
    % Derive Q_gas in particles per second
    Q_gas = (Q_gas_kg/m_ave)/(4*pi*r_comet^2*v_gas);
    gamma_gas = 2;
    
    %% Generate all trajectories and conditions in the ionosphere
    t = ( -t0 : dt : t0 )';
    % Trajectory of A (spherical coordinates)
    [r_A, phi_A, theta_A] = FlybyTrajectory( V_flyby, D_ca_A, phi_ca, theta_ca, t_ca_A, t );
    % Cartesian coordinates of A at closest approach
    [x_ca_A, y_ca_A, z_ca_A] = SphericalToCartesian(D_ca_A, phi_ca, theta_ca); 
    % Cartesian coordinates of probes at closest approach of A
    [r_ca_B1, phi_ca_B1, theta_ca_B1] = CartesianToSpherical(x_ca_A + delta_B1(1), y_ca_A + delta_B1(2), z_ca_A + delta_B1(3));
    [r_ca_B2, phi_ca_B2, theta_ca_B2] = CartesianToSpherical(x_ca_A + delta_B2(1), y_ca_A + delta_B2(2), z_ca_A + delta_B2(3));
    % Trajectories of probes (spherical coordinates)
    [r_B1, phi_B1, theta_B1] = FlybyTrajectory( V_flyby, r_ca_B1, phi_ca_B1, theta_ca_B1, t_ca_A, t );
    [r_B2, phi_B2, theta_B2] = FlybyTrajectory( V_flyby, r_ca_B2, phi_ca_B2, theta_ca_B2, t_ca_A, t );
    
    % Conditions in ionosphere along trajectories
%     f_A = CometReal_OneOverRgamma( r_A, phi_A, theta_A, Q_gas, gamma_gas );
%     f_B1 = CometReal_OneOverRgamma( r_B1, phi_B1, theta_B1, Q_gas, gamma_gas );
%     f_B2 = CometReal_OneOverRgamma( r_B2, phi_B2, theta_B2, Q_gas, gamma_gas );

    jet = struct( 'theta', theta_jet, 'phi', phi_jet, 'dangle', dangle_jet ,'f', f_jet );
%     % Conditions in ionosphere along trajectories
    f_A = CometIonosphere( r_A, phi_A, theta_A, r_comet, Q_gas, v_gas, a, alpha, gamma, jet, t );
    f_B1 = CometIonosphere( r_B1, phi_B1, theta_B1, r_comet, Q_gas, v_gas, a, alpha, gamma, jet, t );
    f_B2 = CometIonosphere( r_B2, phi_B2, theta_B2, r_comet, Q_gas, v_gas, a, alpha, gamma, jet, t );
    
    n_t = length(r_A);
    
    %% Generate measurements with noise
    rng(1);
    f_A_obs = f_A .* ( 1 + rel_error * randn(size(f_A)) );
    f_B1_obs = f_B1 .* ( 1 + rel_error * randn(size(f_A)) );
    f_B2_obs = f_B2 .* ( 1 + rel_error * randn(size(f_A)) );
    
    %% Plot measurements time sequence
    hf1 = findobj( 0, 'Tag', 'RadioSounding - plasma density' );
    if isempty(hf1)
        hf1 = figure( 'Color', 'white', 'Tag', 'RadioSounding - plasma density' );
    else
        clf(hf1);
    end
    ha1 = ...
        axes( ...
            'Parent', hf1, ...
            'Units', 'normalized', 'Position', [ 0.15, 0.15, 0.80, 0.80 ], ...
            'Box', 'on', 'Layer', 'top', ...
            'LineWidth', line_width, ...
            'FontSize', font_size ...
        );
    hold( ha1, 'on' );
    plot( ha1, t, f_A_obs/1e6, 'o', 'Color', [ 0 0 1 ], 'MarkerSize', marker_size );
    plot( ha1, t, f_A/1e6, '-', 'Color', [ 0 0 1 ], 'LineWidth', line_width );
    plot( ha1, t, f_B1_obs/1e6, 'o', 'Color', [ 0 1 0 ], 'MarkerSize', marker_size );
    plot( ha1, t, f_B1/1e6, '-', 'Color', [ 0 1 0 ], 'LineWidth', line_width );
    plot( ha1, t, f_B2_obs/1e6, 'o', 'Color', [ 1 0 0 ], 'MarkerSize', marker_size );
    plot( ha1, t, f_B2/1e6, '-', 'Color', [ 1 0 0 ], 'LineWidth', line_width );
    hxl = get( ha1, 'XLabel' );
    set( hxl, 'String', '$t$ [s]', 'Interpreter', 'latex', 'FontSize', font_size );
    hyl = get( ha1, 'YLabel' );
    set( hyl, 'String', '$f$ [cm$^{-3}$]', 'Interpreter', 'latex', 'FontSize', font_size );
    
    %% Plot problem geometry
    hf2 = findobj( 0, 'Tag', 'RadioSounding - geometry' );
    if isempty(hf2)
        hf2 = figure( 'Color', 'white', 'Tag', 'RadioSounding - geometry' );
    else
        clf(hf2);
    end
    ha2 = ...
        axes( ...
            'Parent', hf2, ...
            'Units', 'normalized', 'Position', [ 0.15, 0.15, 0.80, 0.80 ], ...
            'Box', 'on', 'Layer', 'top', ...
            'LineWidth', line_width, ...
            'FontSize', font_size, ...
            'XLimMode', 'manual', ...
            'XLim', [ -5*D_ca_A, 5*D_ca_A ]/1000, ...
            'YLim', [ -5*D_ca_A, 5*D_ca_A ]/1000,  ...
            'ZLim', [ -5*D_ca_A, 5*D_ca_A ]/1000 ...
        );
    [ xx, yy, zz] = meshgrid( linspace(-5*D_ca_A,5*D_ca_A,100), linspace(-5*D_ca_A,5*D_ca_A,100),  linspace(-5*D_ca_A,5*D_ca_A,100));
    plot_rr = sqrt( xx.^2 + yy.^2 + zz.^2);
    plot_phi = atan2d( yy, xx );
    plot_theta = atan2d(sqrt(xx.^2 + yy.^2), zz);
    plot_ne = CometIonosphere( plot_rr, plot_phi, plot_theta, r_comet, Q_gas, v_gas, a, alpha, gamma, jet );
%     plot_ne = CometReal_OneOverRgamma( plot_rr, plot_phi, plot_theta, Q_gas, gamma_gas );
    xslice = [];   
    yslice = 0;
    zslice = 0;
    slices = slice(ha2, xx/1000, yy/1000, zz/1000,log10(plot_ne/1e6),xslice,yslice,zslice);
    %, 'EdgeColor', 'none'      % zwarte kotjes weg?
    set(slices,'EdgeColor','none');
    view(3)
    hold( ha2, 'on' );
    idx = round(n_t/2);
    plot3( ha2, ...
        -r_A.*sind(theta_A).*cosd(phi_A)/1000, r_A.*sind(theta_A).*sind(phi_A)/1000, r_A.*cosd(theta_A)/1000, 'b-', ...
        -r_B1.*sind(theta_B1).*cosd(phi_B1)/1000, r_B1.*sind(theta_B1).*sind(phi_B1)/1000, r_B1.*cosd(theta_B1)/1000,'g-', ...
        -r_B2.*sind(theta_B2).*cosd(phi_B2)/1000, r_B2.*sind(theta_B2).*sind(phi_B2)/1000, r_B2.*cosd(theta_B2)/1000,'r-', ...
        0, 0, 0, 'ko', ...
        -r_A(idx).*sind(theta_A(idx)).*cosd(phi_A(idx))/1000, r_A(idx).*sind(theta_A(idx)).*sind(phi_A(idx))/1000, r_A(idx).*cosd(theta_A(idx))/1000, 'bo', ...
        -r_B1(idx).*sind(theta_B1(idx)).*cosd(phi_B1(idx))/1000, r_B1(idx).*sind(theta_B1(idx)).*sind(phi_B1(idx))/1000, r_B1(idx).*cosd(theta_B1(idx))/1000, 'go', ...
        -r_B2(idx).*sind(theta_B2(idx)).*cosd(phi_B2(idx))/1000, r_B2(idx).*sind(theta_B2(idx)).*sind(phi_B2(idx))/1000, r_B2(idx).*cosd(theta_B2(idx))/1000, 'ro', ...
        'LineWidth', line_width, 'MarkerSize', marker_size...
    );
    hxl = get( ha2, 'XLabel' );
    set( hxl, 'String', '$x$ [km]', 'Interpreter', 'latex', 'FontSize', font_size );
    hyl = get( ha2, 'YLabel' );
    set( hyl, 'String', '$y$ [km]', 'Interpreter', 'latex', 'FontSize', font_size );
    hzl = get( ha2, 'ZLabel' );
    set( hzl, 'String', '$z$ [km]', 'Interpreter', 'latex', 'FontSize', font_size );
    
    hc = colorbar( ha2 );
    set( get(hc, 'Title' ), 'String', '$f$ [cm$^{-3}$]', 'Interpreter', 'latex', 'FontSize', font_size );
    set( ...
        hc, ...
        'YTick', ...
            log10( ...
                [ ...
                    0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, ...
                    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, ...
                    1, 2, 3, 4, 5, 6, 7, 8, 9, ...
                    10, 20, 30, 40, 50, 60, 70, 80, 90, ...
                    100, 200, 300, 400, 500, 600, 700, 800, 900, ...
                    1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, ...
                    10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, ...
                    100000 ...
                ] ...
            ), ...
        'YTickLabel', ...
            { ...
                '', '', '', '', '', '', '', '', '', ...
                '10^{-1}', '', '', '', '', '', '', '', '', ...
                '10^0', '', '', '', '', '', '', '', '', ...
                '10^1', '', '', '', '', '', '', '', '', ...
                '10^2', '', '', '', '', '', '', '', '', ...
                '10^3', '', '', '', '', '', '', '', '', ...
                '10^4', '', '', '', '', '', '', '', '', ...
                '' ...
            }, ...
        'FontSize', font_size ...
    );

    drawnow;
         
    Q_start = Q_gas*1.33;
    gamma_start = gamma_gas*1.25;
    
    %% Problem 0A:
    % Derive Q, gamma from measurements by A global
    if algo.do_OneOverRgamma_A
        set_t = { t };
        set_r = { r_A };
        set_phi = { phi_A };
        set_f = { f_A_obs };
        set_df = { f_A_obs * rel_error };    
        [ gamma_0A, dgamma_0A, Q_0A, dQ_0A ] = ...
            FindSolution_OneOverRgamma( ...
                false, Q_start, gamma_start, ...
                set_t, set_r, set_phi, set_f, set_df, ...
                tolerance ...
            );
        fprintf( 'Problem 0A\n' );
        fprintf( '\n' );
        fprintf( '      Q                             gamma\n' );
        fprintf( '    : %12.6e +/- %12.6e  %10.6f +/- %10.6f\n', Q_0A, dQ_0A, gamma_0A, dgamma_0A );
        fprintf( '\n' );
    end
    
    %% Plot parameters time sequence
    hf3 = findobj( 0, 'Tag', 'RadioSounding - parameters' );
    if isempty(hf3)
        hf3 = figure( 'Color', 'white', 'Tag', 'RadioSounding - parameters' );
    else
        clf(hf3);
    end
    ha3A = ...
        axes( ...
            'Parent', hf3, ...
            'Units', 'normalized', 'Position', [ 0.15, 0.55, 0.80, 0.35 ], ...
            'Box', 'on', 'Layer', 'top', ...
            'LineWidth', line_width, ...
            'FontSize', font_size ...
        );
    ha3B = ...
        axes( ...
            'Parent', hf3, ...
            'Units', 'normalized', 'Position', [ 0.15, 0.15, 0.80, 0.35 ], ...
            'Box', 'on', 'Layer', 'top', ...
            'LineWidth', line_width, ...
            'FontSize', font_size ...
        );
    hold( ha3A, 'on' );
    plot( ha3A, t, ones(size(t))*gamma_gas, '-', 'Color', [ 0 0 0 ], 'LineWidth', line_width );
    if algo.do_OneOverRgamma_A
        plot( ha3A, [ -1, 1 ]*t0, [ 1, 1]*gamma_0A, '-', 'Color', [ 0.6 0.6 0.6 ], 'LineWidth', 0.5*line_width );
    end
%     hxl = get( ha3A, 'XLabel' );
%     set( hxl, 'String', '$t$ [s]', 'Interpreter', 'latex', 'FontSize', font_size );
    hyl = get( ha3A, 'YLabel' );
    set( hyl, 'String', '$\gamma$', 'Interpreter', 'latex', 'FontSize', font_size );
    set( ha3A, 'YLim', [ 0.5, 1.5 ]*gamma_gas );

    set( ha3B, ...
        'YLim', [ 0.005, 2000 ]*Q_gas, ...
        'YScale', 'log', ...
        'YTickMode', 'manual', ...
        'YTick', [ 10, 100, 1000, 10000, 100000, 1000000 ], ...
        'YTickLabel',{ '10^{1}', '10^{2}', '10^{3}', '10^{4}', '10^{5}', '10^{6}' } ...
    );
    hold( ha3B, 'on' );
    plot( ha3B, t, Q_gas*ones(size(t)), '-', 'Color', [ 0 0 0 ], 'LineWidth', line_width );
    if algo.do_OneOverRgamma_A
        plot( ha3B, [ -1, 1 ]*t0, [ 1, 1]*Q_0A, '-', 'Color', [ 0.6 0.6 0.6 ], 'LineWidth', 0.5*line_width );
    end
    hxl = get( ha3B, 'XLabel' );
    set( hxl, 'String', '$t$ [s]', 'Interpreter', 'latex', 'FontSize', font_size );
    hyl = get( ha3B, 'YLabel' );
    set( hyl, 'String', '$Q$ [kg/s$^{-1}$]', 'Interpreter', 'latex', 'FontSize', font_size );
    set( ha3B, 'YLim', [ 0.5, 1.5 ]*Q_gas );
end

function [ gamma, dgamma, Q, dQ ] = FindSolution_OneOverRgamma( debug, Q0, gamma0, t, r, phi, f, df, tol )
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
            'f', {f}, ...
            'df', {df} ...
        );
    x0 = [ gamma0; log(Q0) ];
    sx = [ 0.1; 1 ];
    [ x, dx ] = ...
        expl_optimization( [], 'FindSolution_OneOverRgamma', x0, sx, tol, @TargetFunction_OneOverRgamma, data, 'enhanced', 'silent' );
    gamma = x(1);
    Q = exp(x(2));
    dgamma = dx(1);
    dQ = Q * dx(2);
end
    
function F = TargetFunction_OneOverRgamma( x, data )
    %% Target function
    % Least-squares formulation given measurements at points
    gamma = x(1);
    if ( gamma < 0 ) 
        F = 1e20 * (-gamma);
    elseif ( gamma > 10 )
        F = 1e20 * (gamma-10);
    else
        Q = exp(x(2));
        F = 0;
        f_model = cell( size(data.f) );
        for k = 1:data.nr_of_sets
            f_model{k} = CometModel_OneOverRgamma( data.r{k}, data.phi{k}, Q, gamma );
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

function f = CometReal_OneOverRgamma( r, phi, theta, Q, gamma ) 
    %% Routine that computes the power law model
    % f = Q / r^gamma
    f = Q ./ r.^(gamma);
end

function sza = SolarZenithAngle( phi, theta) 
    x = -sind(theta) .* cosd(phi);
    y = sind(theta) .* sind(phi);
    z = cosd(theta);
    % Cosine rule: a^2 = b^2 + c^2 - 2*a*c*cos(alpha)
    sza = acosd((sqrt(y.^2 + z.^2).^2 - x.^2 - 1) ./ (-2 * x));
end

function q_factor = JetEnhancement( phi, theta, jet )
    % Unpack jet-related data
    f_jet = jet.f; 
    theta_jet = jet.theta;
    phi_jet = jet.phi;
    dangle_jet = jet.dangle; %half width
    
    q_factor = ones(size(phi));

    [x,y,z] = SphericalToCartesian(1, phi, theta);
    
    n_jet = length(jet.f);
    for i=1:n_jet
        % Compute angular distance
        [xjet,yjet,zjet] = SphericalToCartesian(1, phi_jet(i), theta_jet(i));

        % Angle between [x,y,z] and [xjet, yjet, zjet] ( A.B / |A||B| = cos(A,B) )
        % Dot-product
        X = x(:) ; Y = y(:) ; Z = z(:);
        v = [X Y Z];
        d = sum([xjet,yjet,zjet].*v,2);
        d = reshape(d,size(x));
        % Angle
        dangle = acosd(d ./ (sqrt(x.^2 + y.^2 + z.^2) * sqrt(xjet^2 + yjet^2 + zjet^2)));

        % Compute jet density enhancement factor
        q_factor = q_factor + f_jet(i) * exp( -(dangle/dangle_jet(i)).^2 );

    end
end

function q_gas = GasProduction( phi, theta, r_comet, a, Q_gas, jet )
    % Gas production at surface
    
    % Solar Zenith Angle (calculate with phi and theta) => angle between x-axis and r
    sza = SolarZenithAngle(phi, theta);

    % Compute the flux on the nucleus surface at the subsolar point
    q_factor = JetEnhancement( phi, theta, jet );
    q_gas0 = ( 1 + a ) * q_factor .* Q_gas ./ ( 4 * pi * r_comet.^2 ); %Total gas production devided by surface

    % Compute the flux on the nucleus surface at the base of the streamline
    % as a function of theta
    q_gas = q_gas0 .* ( 1 + a*cosd(sza+180) ) / ( 1 + a ); %sza ipv phi
end

function q_gas = GasProductionNormalized( phi, theta, r_comet, a, Q_gas, jet )
    
    % Unnormalized gas production
    q_gas0 = GasProduction(phi,theta,r_comet,a,Q_gas,jet);
    % Total gas production in model
    integr = integral2(@(x,y) GasProduction(x,y,r_comet,a,Q_gas,jet), 0, 2*pi, 0, pi);
    % Correct gas production so total gas production is Q_gas;
    q_gas = q_gas0 * (Q_gas/integr);


    % Cardioid in 3D parameter curve (asymmetry) => kan gebruiken als
    % visualisatie van de asymmetrie (isodensity oppervlak)
%     funx = @(u,v) cos(u).*((1+a*cos(u))/(1+a));
%     funy = @(u,v) sin(u).*((1+a*cos(u))/(1+a)).*cos(v);
%     funz = @(u,v) sin(u).*((1+a*cos(u))/(1+a)).*sin(v); 
%     fsurf(funx,funy,funz, [0 2*pi 0 pi])

end

function n_e = CometIonosphere( r, phi, theta, r_comet, Q_gas, v_gas, a, alpha, gamma, jet, t ) 
    % Routine that computes total electron density as a function of theta and r

    % Gas production at the surface
    q_gas = GasProductionNormalized( phi, theta, r_comet, a, Q_gas, jet );
    
    %n_n = q_gas / v_gas;
    n_e = ...
        ( r_comet^(2*gamma) * alpha / ( (3-2*gamma) .* v_gas^(gamma+1) ) ) ...
        .* q_gas.^gamma .* ( r - r_comet ).^(3-2*gamma) ./ r.^2;
end

function f = CometModel_OneOverRgamma( r, phi, Q, gamma ) 
    % Routine that computes the power law model
    % f = Q / r^gamma
    f = Q ./ r.^(gamma);
end

function [x, y, z] = SphericalToCartesian(r, phi, theta)
    % Routine that converts spherical coordinates to cartesian coordinates
    x = r .* sind(theta) .* cosd(phi);
    y = r .* sind(theta) .* sind(phi);
    z = r .* cosd(theta);
end

function [r, phi, theta] = CartesianToSpherical(x, y, z)
    % Routine that converts cartesian coordinates to spherical coordinates
    r = sqrt( x.^2 + y.^2 + z.^2);
    phi = atan2d( y, x );
    theta = atan2d(sqrt(x.^2 + y.^2), z);
end

function [ r, phi, theta ] = FlybyTrajectory( Vflyby, D_ca, phi_ca, theta_ca, t_ca, t )
    % Trajectory as a function of time
    
    % Cartesian coordinates closest approach
    [x_ca, y_ca, z_ca] = SphericalToCartesian(D_ca, phi_ca, theta_ca);

    % Cartesian coordinates in time
    x_sc = (t-t_ca) * Vflyby(1) - x_ca;
    y_sc = (t-t_ca) * Vflyby(2) + y_ca;
    z_sc = (t-t_ca) * Vflyby(3) + z_ca;
    
    % Spherical coordinates
    r = sqrt( x_sc.^2 + y_sc.^2 + z_sc.^2);
    phi = atan2d( y_sc, x_sc );
    theta = atan2d(sqrt(x_sc.^2 + y_sc.^2), z_sc);
    
    % Make sure there are no jumps in phi
    n = length(phi);
    for i = 2:n
        if phi(i) > phi(i-1)
            % apparently ascending theta; check whether it wasn't a jump
            if abs(phi(i) - 360 - phi(i-1)) < (phi(i) - phi(i-1))
                % then it was a jump; correct all following theta
                phi(i:end) = phi(i:end) - 360;
            % else
                % smoothly increasing
            end
        else % phi(i) < phi(i-1)
            % apparently descending phi; check whether it wasn't a jump
            if abs(phi(i-1) - 360 - phi(i)) < (phi(i-1) - phi(i))
                % then it was a jump; correct all following phi
                phi(i:end) = phi(i:end) + 360;
            % else
                % smoothly decreasing
            end
        end
    end

    % Make sure there are no jumps in theta
    n = length(theta);
    for i = 2:n
        if theta(i) > theta(i-1)
            % apparently ascending theta; check whether it wasn't a jump
            if abs(theta(i) - 360 - theta(i-1)) < (theta(i) - theta(i-1))
                % then it was a jump; correct all following theta
                theta(i:end) = theta(i:end) - 360;
            % else
                % smoothly increasing
            end
        else % theta(i) < theta(i-1)
            % apparently descending theta; check whether it wasn't a jump
            if abs(theta(i-1) - 360 - theta(i)) < (theta(i-1) - theta(i))
                % then it was a jump; correct all following theta
                theta(i:end) = theta(i:end) + 360;
            % else
                % smoothly decreasing
            end
        end
    end
end
   
function [] = DebugPlot( nr_of_sets, t, f, f_model, txt )
     hf = findobj(0,'Tag','Debug');
     if isempty(hf)
         hf = figure( 'Tag', 'Debug' );
     else
         clf(hf);
     end
     ha = axes(hf);
     for k = 1:nr_of_sets
         semilogy( ha, t{k}, f{k}, 'b.', t{k}, f_model{k}, 'k-' );
         hold( ha, 'on' );
     end
     title( txt );
     pause( 0.001 );
 end


    
    