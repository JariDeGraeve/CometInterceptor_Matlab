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
    
    % Physical parameters
    r_comet = 2000;            % comet radius, in m
    Q_gas_kg = 1000;           % gas production rate, in kg/s
    m_ave = 18 * 1.6e-27;      % average neutral gas particle mass, in kg/particle
    v_gas = 1000;              % neutral gas speed, in m/s
    D_ca_A = 1000000;          % distance at closest approach of A, in m
    t_ca_A = 0;                % time at closest approach of A, in s
    V_flyby = 60000;           % flyby speed, in m/s
    phi_ca = 30;               % azimuth angle at closest approach, in degrees
%    theta_ca = 30;             % elevation angle at closest approach, in degrees
%    sc_inclination = 0;        % inclination of orbit
    a = 0.8;                   % gas coma asymmetry factor, dimensionless
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
    
    % Generate all trajectories and conditions in the ionosphere
    t = ( -t0 : dt : t0 )';
    [ r_A, phi_A ] = FlybyTrajectory( V_flyby, D_ca_A, phi_ca, t_ca_A, t );
    f_A = CometReal_OneOverRgamma( r_A, phi_A, Q_gas, gamma_gas );
    n_t = length(r_A);
    
    % Generate measurements with noise
    rng(1);
    f_A_obs = f_A .* ( 1 + rel_error * randn(size(f_A)) );
    
    % Plot measurements time sequence
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
    hxl = get( ha1, 'XLabel' );
    set( hxl, 'String', '$t$ [s]', 'Interpreter', 'latex', 'FontSize', font_size );
    hyl = get( ha1, 'YLabel' );
    set( hyl, 'String', '$f$ [cm$^{-3}$]', 'Interpreter', 'latex', 'FontSize', font_size );
    
    % Plot problem geometry
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
            'XLim', [ -3*D_ca_A, 3*D_ca_A ]/1000, ...
            'YLim', [ -3*D_ca_A, 3*D_ca_A ]/1000  ...
        );
    [ xx, yy ] = meshgrid( linspace(-3*D_ca_A,3*D_ca_A,100), linspace(-3*D_ca_A,3*D_ca_A,100) );
    plot_rr = sqrt( xx.^2 + yy.^2 );
    plot_phi = atan2d( yy, xx );
    plot_ne = CometModel_OneOverRgamma( plot_rr, plot_phi, Q_gas, gamma_gas );
    surface( ha2, xx/1000, yy/1000, 0*plot_ne, log10(plot_ne/1e6), 'EdgeColor', 'none', 'FaceColor', 'interp' );
    view( 0, 90 );
    hold( ha2, 'on' );
    idx = round(n_t/2);
    plot( ha2, ...
        r_A.*cosd(phi_A)/1000, r_A.*sind(phi_A)/1000, 'b-', ...
        0, 0, 'ko', ...
        r_A(idx).*cosd(phi_A(idx))/1000, r_A(idx).*sind(phi_A(idx))/1000, 'bo', ...
        'LineWidth', line_width, 'MarkerSize', marker_size ...
    );
    hxl = get( ha2, 'XLabel' );
    set( hxl, 'String', '$x$ [km]', 'Interpreter', 'latex', 'FontSize', font_size );
    hyl = get( ha2, 'YLabel' );
    set( hyl, 'String', '$y$ [km]', 'Interpreter', 'latex', 'FontSize', font_size );
    
    hc = colorbar( ha2 );
    set( get(hc, 'Title' ), 'String', '$f$ [cm$^{-3}$]', 'Interpreter', 'latex', 'FontSize', font_size );
    set( ...
        hc, ...
        'YTick', ...
            log10( ...
                [ ...
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
    
    % Problem 0A:
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
    
    % Plot parameters time sequence
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
    % Target function
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

function f = CometReal_OneOverRgamma( r, phi, Q, gamma ) 
    % Routine that computes the power law model
    % f = Q / r^gamma
    f = Q ./ r.^(gamma);
end

function f = CometModel_OneOverRgamma( r, phi, Q, gamma ) 
    % Routine that computes the power law model
    % f = Q / r^gamma
    f = Q ./ r.^(gamma);
end

function [ r, phi ] = FlybyTrajectory( Vflyby, D_ca, phi_ca, t_ca, t )
    % Trajectory as a function of time
    
    % Cartesian coordinates
    x_sc = Vflyby * sind(phi_ca) * ( t - t_ca ) - D_ca * cosd( phi_ca );
    y_sc = Vflyby * cosd(phi_ca) * ( t - t_ca ) + D_ca * sind( phi_ca );
    
    % Polar coordinates
    r = sqrt( x_sc.^2 + y_sc.^2 );
    phi = atan2d( y_sc, x_sc );
    
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


    
    