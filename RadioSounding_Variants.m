function [] = RadioSounding_Variants()

    % Algorithm selection
%     algo = ...
%         struct( ...
%             'do_GammaGlobalQGlobal_A', true, ...
%             'do_GammaGlobalQGlobal_AB', true, ...
%             'do_GammaGlobalQGlobal_ABLOS', true, ...
%             'do_GammaGlobalQLocal_A', true, ...
%             'do_GammaGlobalQLocal_AB', true, ...
%             'do_GammaGlobalQLocal_ABLOS', true, ...
%             'do_GammaLocalQLocal_Time_A', true, ...
%             'do_GammaLocalQLocal_Time_AB', true, ...
%             'do_GammaLocalQLocal_Time_ABLOS', true, ...
%             'do_GammaLocalQLocal_Theta_A', true, ...
%             'do_GammaLocalQLocal_Theta_AB', true, ...
%             'do_GammaLocalQLocal_Theta_ABLOS', true, ...
%             'do_GammaGlobalQiGlobal_A', true, ...
%             'do_GammaGlobalQiGlobal_AB', true, ...
%             'do_GammaGlobalQiGlobal_ABLOS', true ... 
%         );
%     RadioSounding_Core( algo );
    
    algo = ...
        struct( ...
            'do_GammaGlobalQGlobal_A', true, ...
            'do_GammaGlobalQGlobal_AB', true, ...
            'do_GammaGlobalQGlobal_ABLOS', true, ...
            'do_GammaGlobalQLocal_A', false, ...
            'do_GammaGlobalQLocal_AB', false, ...
            'do_GammaGlobalQLocal_ABLOS', false, ...
            'do_GammaLocalQLocal_Time_A', true, ...
            'do_GammaLocalQLocal_Time_AB', true, ...
            'do_GammaLocalQLocal_Time_ABLOS', true, ...
            'do_GammaLocalQLocal_Theta_A', true, ...
            'do_GammaLocalQLocal_Theta_AB', true, ...
            'do_GammaLocalQLocal_Theta_ABLOS', true, ...
            'do_GammaGlobalQiGlobal_A', false, ...
            'do_GammaGlobalQiGlobal_AB', false, ...
            'do_GammaGlobalQiGlobal_ABLOS', false ... 
        );
    RadioSounding_Core( algo );
    
end

function [] = RadioSounding_Core( algo )
    
    % Physical parameters
    r_comet = 2000;            % comet radius, in m
    Q_gas_kg = 1200;           % gas production rate, in kg/s
    m_ave = 18 * 1.6e-27;      % average neutral gas particle mass, in kg/particle
    v_gas = 1000;              % neutral gas speed, in m/s
    D_ca_A = 1000000;          % distance at closest approach of A, in m
    D_ca_B = 800000;           % distance at closest approach of B, in m
    t_ca_A = 0;                % time at closest approach of A, in s
    t_ca_B = 5;                % time at closest approach of B, in s
    V_flyby = 60000;           % flyby speed, in m/s
    theta_ca = 30;             % solar zenith angle at closest approach, in degrees
    gamma = 0.8;               % ionization rate exponent, dimensionless
    alpha = 0.00002;           % ionization rate, in 1/s
    a = 0.8;                   % gas coma asymmetry factor, dimensionless
    theta_jet = [ -30,80];     % solar zenith angle of center of gas jet, in degrees
    dtheta_jet = [ 15, 10 ];   % angular half-width of gas jet, in degrees
    f_jet = 1*[3,5];           % density contrast in jet, dimensionless
    dn_e_rel = 0.05;           % relative measurement error, dimensionless
    dn_e_LOS_rel = 0.05;       % relative LOS measurement error, dimensionless
    
    % Numerical parameters
    t0 = 150;                  % limits of time interval, in s
    dt = 0.2;                  % measurement time resolution, in s
    N_sample = 10;             % number of sampling points along LOS
    dt_homogeneity = 5;        % time homogeneity scale, in s
    dtheta_homogeneity = 1.0;  % theta homogeneity scale, in degrees
    tolerance = 1e-6;          % tolerance for solving optimization problem
    dt_ls = 5;                 % time interval for which to compute least-squares solutions
    n_theta_i = 72;              

    % Plot parameters
    font_size = 14;
    marker_size = 8;
    line_width = 2;
    
    % Derive Q_gas in particles per second
    Q_gas = Q_gas_kg/m_ave;
    
    % Prepare structures with jet properties
    jet = struct( 'theta', theta_jet, 'dtheta', dtheta_jet, 'f', f_jet );
    
    % Generate all trajectories and conditions in the ionosphere
    t = ( -t0 : dt : t0 )';
    [ r_A, theta_A ] = FlybyTrajectory( V_flyby, D_ca_A, theta_ca, t_ca_A, t );
    n_e_A = CometIonosphere( r_A, theta_A, r_comet, Q_gas, v_gas, a, alpha, gamma, jet );
    [ r_B, theta_B ] = FlybyTrajectory( V_flyby, D_ca_B, theta_ca, t_ca_B, t );
    n_e_B = CometIonosphere( r_B, theta_B, r_comet, Q_gas, v_gas, a, alpha, gamma, jet );
    n_t = length(r_A);
    n_e_LOS = zeros( n_t, 1 );
    n_e_ave_LOS = zeros( n_t, 1 );
    scale = linspace( 0, 1, N_sample )';
    for i = 1:n_t
        x_A = r_A(i) .* cosd(theta_A(i));
        y_A = r_A(i) .* sind(theta_A(i));
        x_B = r_B(i) .* cosd(theta_B(i));
        y_B = r_B(i) .* sind(theta_B(i));
        x_i = x_A + scale*(x_B-x_A);
        y_i = y_A + scale*(y_B-y_A);
        r_i = sqrt( x_i.^2 + y_i.^2 );
        theta_i = atan2d( y_i, x_i );
        d_AB = sqrt( ( x_A - x_B ).^2 + ( y_A - y_B ).^2 );
        n_e_i = CometIonosphere( r_i, theta_i, r_comet, Q_gas, v_gas, a, alpha, gamma, jet );
        n_e_ave_LOS(i) = sum(n_e_i)/N_sample;
        n_e_LOS(i) = d_AB * n_e_ave_LOS(i);
    end
    
    % Generate measurements with noise
    rng(1);
    n_e_A_obs = n_e_A .* ( 1 + dn_e_rel * randn(size(n_e_A)) );
    n_e_B_obs = n_e_B .* ( 1 + dn_e_rel * randn(size(n_e_B)) );
    the_rand = randn(size(n_e_LOS));
    n_e_ave_LOS_obs = n_e_ave_LOS .* ( 1 + dn_e_LOS_rel * the_rand );
    n_e_LOS_obs = n_e_LOS .* ( 1 + dn_e_LOS_rel * the_rand );
    
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
    plot( ha1, t, n_e_A_obs/1e6, 'o', 'Color', [ 0 0 1 ], 'MarkerSize', marker_size );
    plot( ha1, t, n_e_B_obs/1e6, 'd', 'Color', [ 0 0.5 0 ], 'MarkerSize', marker_size );
    plot( ha1, t, n_e_ave_LOS_obs/1e6, 's', 'Color', [ 1 0 0 ], 'MarkerSize', marker_size );
    plot( ha1, t, n_e_A/1e6, '-', 'Color', [ 0 0 1 ], 'LineWidth', line_width );
    plot( ha1, t, n_e_B/1e6, '-', 'Color', [ 0 0.5 0 ], 'LineWidth', line_width );
    plot( ha1, t, n_e_ave_LOS/1e6, '-', 'Color', [ 1 0 0 ], 'LineWidth', line_width );
    hxl = get( ha1, 'XLabel' );
    set( hxl, 'String', '$t$ [s]', 'Interpreter', 'latex', 'FontSize', font_size );
    hyl = get( ha1, 'YLabel' );
    set( hyl, 'String', '$n_e$ [cm$^{-3}$]', 'Interpreter', 'latex', 'FontSize', font_size );
    
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
    plot_theta = atan2d( yy, xx );
    plot_ne = CometIonosphere( plot_rr, plot_theta, r_comet, Q_gas, v_gas, a, alpha, gamma, jet );
    surface( ha2, xx/1000, yy/1000, 0*plot_ne, log10(plot_ne/1e6), 'EdgeColor', 'none', 'FaceColor', 'interp' );
    view( 0, 90 );
    hold( ha2, 'on' );
    idx = round(n_t/2);
    plot( ha2, ...
        r_A.*cosd(theta_A)/1000, r_A.*sind(theta_A)/1000, 'b-', ...
        r_B.*cosd(theta_B)/1000, r_B.*sind(theta_B)/1000, 'g-', ...
        0, 0, 'ko', ...
        r_A(idx).*cosd(theta_A(idx))/1000, r_A(idx).*sind(theta_A(idx))/1000, 'bo', ...
        r_B(idx).*cosd(theta_B(idx))/1000, r_B(idx).*sind(theta_B(idx))/1000, 'go', ...
        'LineWidth', line_width, 'MarkerSize', marker_size ...
    );
    hxl = get( ha2, 'XLabel' );
    set( hxl, 'String', '$x$ [km]', 'Interpreter', 'latex', 'FontSize', font_size );
    hyl = get( ha2, 'YLabel' );
    set( hyl, 'String', '$y$ [km]', 'Interpreter', 'latex', 'FontSize', font_size );
    
    hc = colorbar( ha2 );
    set( get(hc, 'Title' ), 'String', '$n_e$ [cm$^{-3}$]', 'Interpreter', 'latex', 'FontSize', font_size );
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
                    
    gamma_start = 1.0;
    Q_start = 1000/m_ave;
    
    % Problem 0A:
    % Derive Q, gamma from measurements by A global
    if algo.do_GammaGlobalQGlobal_A
        set_t = { t };
        set_r = { r_A };
        set_theta = { theta_A };
        set_n_e = { n_e_A_obs };
        set_dn_e = { n_e_A_obs * dn_e_rel };    
        [ gamma_0A, dgamma_0A, Q_0A, dQ_0A ] = ...
            FindSolutionGammaGlobalQGlobal( ...
                false, gamma_start, Q_start, ...
                r_comet, v_gas, a, alpha, m_ave, ...
                set_t, set_r, set_theta, set_n_e, set_dn_e, ...
                tolerance ...
            );
        fprintf( 'Problem 0A\n' );
        fprintf( '\n' );
        fprintf( '      Q                             gamma\n' );
        fprintf( '    : %12.6e +/- %12.6e  %10.6f +/- %10.6f\n', Q_0A*m_ave, dQ_0A*m_ave, gamma_0A, dgamma_0A );
        fprintf( '\n' );
    end
        
    % Problem 0AB:
    % Derive Q, gamma from measurements by A and B global
    if algo.do_GammaGlobalQGlobal_AB
        set_t = { t; t };
        set_r = { r_A; r_B };
        set_theta = { theta_A; theta_B };
        set_n_e = { n_e_A_obs; n_e_B_obs };
        set_dn_e = { dn_e_rel * n_e_A_obs; dn_e_rel * n_e_B_obs };    
        [ gamma_0AB, dgamma_0AB, Q_0AB, dQ_0AB ] = ...
            FindSolutionGammaGlobalQGlobal( ...
                false, gamma_start, Q_start, ...
                r_comet, v_gas, a, alpha, m_ave, ...
                set_t, set_r, set_theta, set_n_e, set_dn_e, ...
                tolerance ...
            );
        fprintf( 'Problem 0AB\n' );
        fprintf( '\n' );
        fprintf( '      Q                             gamma\n' );
        fprintf( '    : %12.6e +/- %12.6e  %10.6f +/- %10.6f\n', Q_0AB*m_ave, dQ_0AB*m_ave, gamma_0AB, dgamma_0AB );
        fprintf( '\n' );
    end
        
    % Problem 0ABLOS:
    % Derive Q, gamma from measurements by A and A-B global
    if algo.do_GammaGlobalQGlobal_ABLOS
        set_t = { t; t };
        set_r = { r_A; r_B };
        set_theta = { theta_A; theta_B };
        set_n_e = { n_e_A_obs; n_e_LOS_obs };
        set_dn_e = { dn_e_rel * n_e_A_obs; dn_e_LOS_rel * n_e_LOS_obs };    
        [ gamma_0ABLOS, dgamma_0ABLOS, Q_0ABLOS, dQ_0ABLOS ] = ...
            FindSolutionGammaGlobalQGlobalLOS( ...
                false, gamma_start, Q_start, ...
                r_comet, v_gas, a, alpha, m_ave, ...
                set_t, set_r, set_theta, set_n_e, set_dn_e, ...
                tolerance ...
            );
        fprintf( 'Problem 0ABLOS\n' );
        fprintf( '\n' );
        fprintf( '      Q                             gamma\n' );
        fprintf( '    : %12.6e +/- %12.6e  %10.6f +/- %10.6f\n', Q_0ABLOS*m_ave, dQ_0ABLOS*m_ave, gamma_0ABLOS, dgamma_0ABLOS );
        fprintf( '\n' );
    end
        
    % Problem 1A:
    % Derive gamma_ave, Q_ave from measurements by A at a given time
    if algo.do_GammaLocalQLocal_Time_A
        t_ls = ( -t0 : dt_ls : t0 )';
        n_t_ls = length(t_ls);
        gamma_1A = zeros( n_t_ls, 1 );
        dgamma_1A = zeros( n_t_ls, 1 );
        Q_1A = zeros( n_t_ls, 1 );
        dQ_1A = zeros( n_t_ls, 1 );
        fprintf( 'Problem 1A\n' );
        fprintf( '\n' );
        fprintf( '      Q                             gamma\n' );
        for i = 1:n_t_ls
            t_ref = t_ls(i);
            [ factor, selection ] = HomogeneitySelection( t, t_ref, dt_homogeneity );           
            set_t = { t(selection) };
            set_r = { r_A(selection) };
            set_theta = { theta_A(selection) };
            set_n_e = { n_e_A_obs(selection) };
            set_dn_e = ...
                    { dn_e_rel * n_e_A_obs(selection) .* factor(selection) };
            [ gamma_1A(i), dgamma_1A(i), Q_1A(i), dQ_1A(i) ] = ...
                FindSolutionGammaGlobalQGlobal( ...
                    false, gamma_start, Q_start, ...
                    r_comet, v_gas, a, alpha, m_ave, ...
                    set_t, set_r, set_theta, set_n_e, set_dn_e, ...
                    tolerance ...
                );
            fprintf( '    : %12.6e +/- %12.6e  %10.6f +/- %10.6f\n', Q_1A(i)*m_ave, dQ_1A(i)*m_ave, gamma_1A(i), dgamma_1A(i) );
        end
        fprintf( '\n' );
    end

    % Problem 1AB:
    % Derive gamma_ave, Q_ave from measurements by A and B at a given time
    if algo.do_GammaLocalQLocal_Time_AB
        t_ls = ( -t0 : dt_ls : t0 )';
        n_t_ls = length(t_ls);
        gamma_1AB = zeros( n_t_ls, 1 );
        dgamma_1AB = zeros( n_t_ls, 1 );
        Q_1AB = zeros( n_t_ls, 1 );
        dQ_1AB = zeros( n_t_ls, 1 );
        fprintf( 'Problem 1AB\n' );
        fprintf( '\n' );
        fprintf( '      Q                             gamma\n' );
        for i = 1:n_t_ls
            t_ref = t_ls(i);
            [ factor, selection ] = HomogeneitySelection( t, t_ref, dt_homogeneity );           
            set_t = { t(selection); t(selection) };
            set_r = { r_A(selection); r_B(selection) };
            set_theta = { theta_A(selection); theta_B(selection) };
            set_n_e = { n_e_A_obs(selection); n_e_B_obs(selection) };
            set_dn_e = ...
                { ...
                    dn_e_rel * n_e_A_obs(selection) .* factor(selection); ...
                    dn_e_rel * n_e_B_obs(selection) .* factor(selection)  ...
                };
            [ gamma_1AB(i), dgamma_1AB(i), Q_1AB(i), dQ_1AB(i) ] = ...
                FindSolutionGammaGlobalQGlobal( ...
                    false, gamma_start, Q_start, ...
                    r_comet, v_gas, a, alpha, m_ave, ...
                    set_t, set_r, set_theta, set_n_e, set_dn_e, ...
                    tolerance ...
                );
            fprintf( '    : %12.6e +/- %12.6e  %10.6f +/- %10.6f\n', Q_1AB(i)*m_ave, dQ_1AB(i)*m_ave, gamma_1AB(i), dgamma_1AB(i) );
        end
        fprintf( '\n' );
    end

    % Problem 1ABLOS:
    % Derive gamma_ave, Q_ave from measurements by A and A-B at a given time
    if algo.do_GammaLocalQLocal_Time_ABLOS
        t_ls = ( -t0 : dt_ls : t0 )';
        n_t_ls = length(t_ls);
        gamma_1ABLOS = zeros( n_t_ls, 1 );
        dgamma_1ABLOS = zeros( n_t_ls, 1 );
        Q_1ABLOS = zeros( n_t_ls, 1 );
        dQ_1ABLOS = zeros( n_t_ls, 1 );
        fprintf( 'Problem 1ABLOS\n' );
        fprintf( '\n' );
        fprintf( '      Q                             gamma\n' );
        for i = 1:n_t_ls
            t_ref = t_ls(i);
            [ factor, selection ] = HomogeneitySelection( t, t_ref, dt_homogeneity );           
            set_t = { t(selection); t(selection) };
            set_r = { r_A(selection); r_B(selection) };
            set_theta = { theta_A(selection); theta_B(selection) };
            set_n_e = { n_e_A_obs(selection); n_e_LOS_obs(selection) };
            set_dn_e = ...
                { ...
                    dn_e_rel * n_e_A_obs(selection) .* factor(selection); ...
                    dn_e_LOS_rel * n_e_LOS_obs(selection) .* factor(selection)  ...
                };
            [ gamma_1ABLOS(i), dgamma_1ABLOS(i), Q_1ABLOS(i), dQ_1ABLOS(i) ] = ...
                FindSolutionGammaGlobalQGlobalLOS( ...
                    false, gamma_start, Q_start, ...
                    r_comet, v_gas, a, alpha, m_ave, ...
                    set_t, set_r, set_theta, set_n_e, set_dn_e, ...
                    tolerance ...
                );
            fprintf( '    : %12.6e +/- %12.6e  %10.6f +/- %10.6f\n', Q_1ABLOS(i)*m_ave, dQ_1ABLOS(i)*m_ave, gamma_1ABLOS(i), dgamma_1ABLOS(i) );
        end
        fprintf( '\n' );
    end

    % Problem 4A:
    % Derive gamma_ave, Q_ave from measurements by A at a given theta
    if algo.do_GammaLocalQLocal_Theta_A
        t_ls = ( -t0 : dt_ls : t0 )';
        n_t_ls = length(t_ls);
        gamma_4A = zeros( n_t_ls, 1 );
        dgamma_4A = zeros( n_t_ls, 1 );
        Q_4A = zeros( n_t_ls, 1 );
        dQ_4A = zeros( n_t_ls, 1 );
        fprintf( 'Problem 4A\n' );
        fprintf( '\n' );
        fprintf( '      Q                             gamma\n' );
        for i = 1:n_t_ls
            t_ref = t_ls(i);
            theta_ref = interp1( t, theta_A, t_ref );
            [ factor, selection ] = HomogeneitySelection( theta_A, theta_ref, dtheta_homogeneity );           
            set_t = { t(selection) };
            set_r = { r_A(selection) };
            set_theta = { theta_A(selection) };
            set_n_e = { n_e_A_obs(selection) };
            set_dn_e = ...
                    { dn_e_rel * n_e_A_obs(selection) .* factor(selection) };
            [ gamma_4A(i), dgamma_4A(i), Q_4A(i), dQ_4A(i) ] = ...
                FindSolutionGammaGlobalQGlobal( ...
                    false, gamma_start, Q_start, ...
                    r_comet, v_gas, a, alpha, m_ave, ...
                    set_t, set_r, set_theta, set_n_e, set_dn_e, ...
                    tolerance ...
                );
            fprintf( '    : %12.6e +/- %12.6e  %10.6f +/- %10.6f\n', Q_4A(i)*m_ave, dQ_4A(i)*m_ave, gamma_4A(i), dgamma_4A(i) );
        end
        fprintf( '\n' );
    end

    % Problem 4AB:
    % Derive gamma_ave, Q_ave from measurements by A and B at a given theta
    if algo.do_GammaLocalQLocal_Theta_AB
        t_ls = ( -t0 : dt_ls : t0 )';
        n_t_ls = length(t_ls);
        gamma_4AB = zeros( n_t_ls, 1 );
        dgamma_4AB = zeros( n_t_ls, 1 );
        Q_4AB = zeros( n_t_ls, 1 );
        dQ_4AB = zeros( n_t_ls, 1 );
        fprintf( 'Problem 4AB\n' );
        fprintf( '\n' );
        fprintf( '      Q                             gamma\n' );
        for i = 1:n_t_ls
            t_ref = t_ls(i);
            theta_ref = interp1( t, theta_A, t_ref );
            [ factor_A, selection_A ] = HomogeneitySelection( theta_A, theta_ref, dtheta_homogeneity );           
            [ factor_B, selection_B ] = HomogeneitySelection( theta_B, theta_ref, dtheta_homogeneity );           
            set_t = { t(selection_A); t(selection_B) };
            set_r = { r_A(selection_A); r_B(selection_B) };
            set_theta = { theta_A(selection_A); theta_B(selection_B) };
            set_n_e = { n_e_A_obs(selection_A); n_e_B_obs(selection_B) };
            set_dn_e = ...
                { ...
                    dn_e_rel * n_e_A_obs(selection_A) .* factor_A(selection_A); ...
                    dn_e_rel * n_e_B_obs(selection_B) .* factor_B(selection_B)  ...
                };
            [ gamma_4AB(i), dgamma_4AB(i), Q_4AB(i), dQ_4AB(i) ] = ...
                FindSolutionGammaGlobalQGlobal( ...
                    false, gamma_start, Q_start, ...
                    r_comet, v_gas, a, alpha, m_ave, ...
                    set_t, set_r, set_theta, set_n_e, set_dn_e, ...
                    tolerance ...
                );
            fprintf( '    : %12.6e +/- %12.6e  %10.6f +/- %10.6f\n', Q_4AB(i)*m_ave, dQ_4AB(i)*m_ave, gamma_4AB(i), dgamma_4AB(i) );
        end
        fprintf( '\n' );
    end

    % Problem 4ABLOS:
    % Derive gamma_ave, Q_ave from measurements by A and A-B at a given time
    if algo.do_GammaLocalQLocal_Theta_ABLOS
        t_ls = ( -t0 : dt_ls : t0 )';
        n_t_ls = length(t_ls);
        gamma_4ABLOS = zeros( n_t_ls, 1 );
        dgamma_4ABLOS = zeros( n_t_ls, 1 );
        Q_4ABLOS = zeros( n_t_ls, 1 );
        dQ_4ABLOS = zeros( n_t_ls, 1 );
        fprintf( 'Problem 4ABLOS\n' );
        fprintf( '\n' );
        fprintf( '      Q                             gamma\n' );
        for i = 1:n_t_ls
            t_ref = t_ls(i);
            theta_ref = interp1( t, theta_A, t_ref );
            [ factor, selection ] = HomogeneitySelection( theta_A, theta_ref, dtheta_homogeneity );           
            set_t = { t(selection); t(selection) };
            set_r = { r_A(selection); r_B(selection) };
            set_theta = { theta_A(selection); theta_B(selection) };
            set_n_e = { n_e_A_obs(selection); n_e_LOS_obs(selection) };
            set_dn_e = ...
                { ...
                    dn_e_rel * n_e_A_obs(selection) .* factor(selection); ...
                    dn_e_LOS_rel * n_e_LOS_obs(selection) .* factor(selection)  ...
                };
            [ gamma_4ABLOS(i), dgamma_4ABLOS(i), Q_4ABLOS(i), dQ_4ABLOS(i) ] = ...
                FindSolutionGammaGlobalQGlobalLOS( ...
                    false, gamma_start, Q_start, ...
                    r_comet, v_gas, a, alpha, m_ave, ...
                    set_t, set_r, set_theta, set_n_e, set_dn_e, ...
                    tolerance ...
                );
            fprintf( '    : %12.6e +/- %12.6e  %10.6f +/- %10.6f\n', Q_4ABLOS(i)*m_ave, dQ_4ABLOS(i)*m_ave, gamma_4ABLOS(i), dgamma_4ABLOS(i) );
        end
        fprintf( '\n' );
    end

    % Problem 2A:
    % Derive gamma, Q_ave from measurements by A at a given theta
    if algo.do_GammaGlobalQLocal_A
        t_ls = ( -t0 : dt_ls : t0 )';
        n_t_ls = length(t_ls);
        set_t = { t };
        set_r = { r_A };
        set_theta = { theta_A };
        set_n_e = { n_e_A_obs };
        set_dn_e = { n_e_A_obs * dn_e_rel };
        [ gamma_2A, dgamma_2A, Q_2A, dQ_2A ] = ...
            FindSolutionGammaGlobalQLocal( ...
                false, gamma_start, Q_start, t_ls, dtheta_homogeneity, r_comet, v_gas, a, alpha, m_ave, ...
                set_t, set_r, set_theta, set_n_e, set_dn_e, tolerance ...
            );    
        fprintf( 'Problem 2A\n' );
        fprintf( '\n' );
        fprintf( '      Q                             gamma\n' );
        for i = 1:n_t_ls
            fprintf( '    : %12.6e +/- %12.6e  %10.6f +/- %10.6f\n', Q_2A(i)*m_ave, dQ_2A(i)*m_ave, gamma_2A, dgamma_2A );
        end
        fprintf( '\n' );
    end

    % Problem 2AB:
    % Derive gamma, Q_ave from measurements by A and B at a given theta
    if algo.do_GammaGlobalQLocal_AB
        t_ls = ( -t0 : dt_ls : t0 )';
        n_t_ls = length(t_ls);
        set_t = { t; t };
        set_r = { r_A; r_B };
        set_theta = { theta_A; theta_B };
        set_n_e = { n_e_A_obs; n_e_B_obs };
        set_dn_e = { n_e_A_obs * dn_e_rel; n_e_B_obs * dn_e_rel };
        [ gamma_2AB, dgamma_2AB, Q_2AB, dQ_2AB ] = ...
            FindSolutionGammaGlobalQLocal( ...
                false,  gamma_start, Q_start, t_ls, dtheta_homogeneity, r_comet, v_gas, a, alpha, m_ave, ...
                set_t, set_r, set_theta, set_n_e, set_dn_e, tolerance ...
            );    
        fprintf( 'Problem 2AB\n' );
        fprintf( '\n' );
        fprintf( '      Q                             gamma\n' );
        for i = 1:n_t_ls
            fprintf( '    : %12.6e +/- %12.6e  %10.6f +/- %10.6f\n', Q_2AB(i)*m_ave, dQ_2AB(i)*m_ave, gamma_2AB, dgamma_2AB );
        end
        fprintf( '\n' );
    end

    % Problem 2ABLOS:
    % Derive gamma, Q_ave from measurements by A and B at a given theta
    if algo.do_GammaGlobalQLocal_ABLOS
        t_ls = ( -t0 : dt_ls : t0 )';
        n_t_ls = length(t_ls);
        set_t = { t; t };
        set_r = { r_A; r_B };
        set_theta = { theta_A; theta_B };
        set_n_e = { n_e_A_obs; n_e_LOS_obs };
        set_dn_e = { n_e_A_obs * dn_e_rel; n_e_LOS_obs * dn_e_LOS_rel };
        [ gamma_2ABLOS, dgamma_2ABLOS, Q_2ABLOS, dQ_2ABLOS ] = ...
            FindSolutionGammaGlobalQLocalLOS( ...
                false, gamma_start, Q_start, t_ls, dtheta_homogeneity, r_comet, v_gas, a, alpha, m_ave, ...
                set_t, set_r, set_theta, set_n_e, set_dn_e, tolerance ...
            );    
        fprintf( 'Problem 2ABLOS\n' );
        fprintf( '\n' );
        fprintf( '      Q                             gamma\n' );
        for i = 1:n_t_ls
            fprintf( '    : %12.6e +/- %12.6e  %10.6f +/- %10.6f\n', Q_2ABLOS(i)*m_ave, dQ_2ABLOS(i)*m_ave, gamma_2ABLOS, dgamma_2ABLOS );
        end
        fprintf( '\n' );
    end

    % Problem 3A:
    % Derive gamma and Qi from measurements by A
    if algo.do_GammaGlobalQiGlobal_A
        theta_i = linspace(0,360,n_theta_i+1)';
        theta_i = theta_i(1:n_theta_i,1);
        n_i = length(theta_i);
        set_t = { t };
        set_r = { r_A };
        set_theta = { theta_A };
        set_n_e = { n_e_A_obs };
        set_dn_e = { n_e_A_obs * dn_e_rel };
        fprintf( 'Problem 3A\n' );
        fprintf( '\n' );
        fprintf( '      Q                             gamma\n' );
        [ gamma_3A, dgamma_3A, Qi_3A, dQi_3A ] = ...
            FindSolutionGammaGlobalQiGlobal( ...
                false, gamma_start, Q_start * ones( n_i, 1 ), ...
                theta_i, ...
                r_comet, v_gas, a, alpha, m_ave, ...
                set_t, set_r, set_theta, set_n_e, set_dn_e, ...
                tolerance ...            
            );
        for i = 1:n_i
            fprintf( '%03.0f : %12.6e +/- %12.6e  %10.6f +/- %10.6f\n', theta_i(i), Qi_3A(i)*m_ave, dQi_3A(i)*m_ave, gamma_3A, dgamma_3A );
        end
        fprintf( '\n' );
        t_ls = ( -t0 : dt_ls : t0 )';
        ipl_Qi_3A = InterpolateQiFromThetaToTime( theta_i, Qi_3A, interp1( t, theta_A, t_ls ) );
    end
    
    % Problem 3AB:
    % Derive gamma and Qi from measurements by A and B
    if algo.do_GammaGlobalQiGlobal_AB
        theta_i = linspace(0,360,n_theta_i+1)';
        theta_i = theta_i(1:n_theta_i,1);
        n_i = length(theta_i);
        set_t = { t; t };
        set_r = { r_A; r_B };
        set_theta = { theta_A; theta_B };
        set_n_e = { n_e_A_obs; n_e_B_obs };
        set_dn_e = { n_e_A_obs * dn_e_rel; n_e_B_obs * dn_e_rel };
        fprintf( 'Problem 3AB\n' );
        fprintf( '\n' );
        fprintf( '      Q                             gamma\n' );
        [ gamma_3AB, dgamma_3AB, Qi_3AB, dQi_3AB ] = ...
            FindSolutionGammaGlobalQiGlobal( ...
                false, gamma_start, Q_start * ones( n_i, 1 ), ...
                theta_i, ...
                r_comet, v_gas, a, alpha, m_ave, ...
                set_t, set_r, set_theta, set_n_e, set_dn_e, ...
                tolerance ...            
            );
        for i = 1:n_i
            fprintf( '%03.0f : %12.6e +/- %12.6e  %10.6f +/- %10.6f\n', theta_i(i), Qi_3AB(i)*m_ave, dQi_3AB(i)*m_ave, gamma_3AB, dgamma_3AB );
        end
        fprintf( '\n' );
        t_ls = ( -t0 : dt_ls : t0 )';
        ipl_Qi_3AB_A = InterpolateQiFromThetaToTime( theta_i, Qi_3AB, interp1( t, theta_A, t_ls ) );    
        ipl_Qi_3AB_B = InterpolateQiFromThetaToTime( theta_i, Qi_3AB, interp1( t, theta_B, t_ls ) );    
    end
    
    % Problem 3ABLOS:
    % Derive gamma and Qi from measurements by A and A-B
    if algo.do_GammaGlobalQiGlobal_ABLOS
        theta_i = linspace(0,360,n_theta_i+1)';
        theta_i = theta_i(1:n_theta_i,1);
        n_i = length(theta_i);
        set_t = { t; t };
        set_r = { r_A; r_B };
        set_theta = { theta_A; theta_B };
        set_n_e = { n_e_A_obs; n_e_LOS_obs };
        set_dn_e = { n_e_A_obs * dn_e_rel; n_e_LOS_obs * dn_e_LOS_rel };
        fprintf( 'Problem 3ABLOS\n' );
        fprintf( '\n' );
        fprintf( '      Q                             gamma\n' );
        [ gamma_3ABLOS, dgamma_3ABLOS, Qi_3ABLOS, dQi_3ABLOS ] = ...
            FindSolutionGammaGlobalQiGlobalLOS( ...
                false, gamma_start, Q_start * ones( n_i, 1 ), ...
                theta_i, ...
                r_comet, v_gas, a, alpha, m_ave, ...
                set_t, set_r, set_theta, set_n_e, set_dn_e, ...
                tolerance ...            
            );
        for i = 1:n_i
            fprintf( '%03.0f : %12.6e +/- %12.6e  %10.6f +/- %10.6f\n', theta_i(i), Qi_3ABLOS(i)*m_ave, dQi_3ABLOS(i)*m_ave, gamma_3ABLOS, dgamma_3ABLOS );
        end
        fprintf( '\n' );
        t_ls = ( -t0 : dt_ls : t0 )';
        ipl_Qi_3ABLOS_A = InterpolateQiFromThetaToTime( theta_i, Qi_3ABLOS, interp1( t, theta_A, t_ls ) );    
        ipl_Qi_3ABLOS_B = InterpolateQiFromThetaToTime( theta_i, Qi_3ABLOS, interp1( t, theta_B, t_ls ) );    
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
    plot( ha3A, t, ones(size(t))*gamma, '-', 'Color', [ 0 0 0 ], 'LineWidth', line_width );
    if algo.do_GammaGlobalQGlobal_A
        plot( ha3A, [ -1, 1 ]*t0, [ 1, 1]*gamma_0A, '-', 'Color', [ 0.6 0.6 0.6 ], 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaGlobalQGlobal_AB
        plot( ha3A, [ -1, 1 ]*t0, [ 1, 1]*gamma_0AB, '--', 'Color', [ 0.6 0.6 0.6 ], 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaGlobalQGlobal_ABLOS
        plot( ha3A, [ -1, 1 ]*t0, [ 1, 1]*gamma_0ABLOS, '.-', 'Color', [ 0.6 0.6 0.6 ], 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaLocalQLocal_Time_A
        plot( ha3A, t_ls, gamma_1A, 'o-', 'Color', [ 0 0 1 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaLocalQLocal_Time_AB
        plot( ha3A, t_ls, gamma_1AB, 'd-', 'Color', [ 0 0 1 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaLocalQLocal_Time_ABLOS
        plot( ha3A, t_ls, gamma_1ABLOS, 's-', 'Color', [ 0 0 1 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaLocalQLocal_Theta_A
        plot( ha3A, t_ls, gamma_4A, 'o-', 'Color', [ 0.5 0 0.5 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaLocalQLocal_Theta_AB
        plot( ha3A, t_ls, gamma_4AB, 'd-', 'Color', [ 0.5 0 0.5 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaLocalQLocal_Theta_ABLOS
        plot( ha3A, t_ls, gamma_4ABLOS, 's-', 'Color', [ 0.5 0 0.5 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaGlobalQLocal_A
        plot( ha3A, [ -1, 1 ]*t0, [ 1, 1]*gamma_2A, '-', 'Color', [ 0 0.5 0 ], 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaGlobalQLocal_AB
        plot( ha3A, [ -1, 1 ]*t0, [ 1, 1]*gamma_2AB, '--', 'Color', [ 0 0.5 0 ], 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaGlobalQLocal_ABLOS
        plot( ha3A, [ -1, 1 ]*t0, [ 1, 1]*gamma_2ABLOS, '.-', 'Color', [ 0 0.5 0 ], 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaGlobalQiGlobal_A
        plot( ha3A, [ -1, 1 ]*t0, [ 1, 1]*gamma_3A, '-', 'Color', [ 1 0 0 ], 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaGlobalQiGlobal_AB
        plot( ha3A, [ -1, 1 ]*t0, [ 1, 1]*gamma_3AB, '--', 'Color', [ 1 0 0 ], 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaGlobalQiGlobal_ABLOS
        plot( ha3A, [ -1, 1 ]*t0, [ 1, 1]*gamma_3ABLOS, '.-', 'Color', [ 1 0 0 ], 'LineWidth', 0.5*line_width );
    end
%     hxl = get( ha3A, 'XLabel' );
%     set( hxl, 'String', '$t$ [s]', 'Interpreter', 'latex', 'FontSize', font_size );
    hyl = get( ha3A, 'YLabel' );
    set( hyl, 'String', '$\gamma$', 'Interpreter', 'latex', 'FontSize', font_size );
    set( ha3A, 'YLim', [ 0.5, 1.2 ] );

    hold( ha3B, 'on' );
    plot( ha3B, t, Q_gas_kg*JetEnhancement( theta_A, jet ), '-', 'Color', [ 0 0 0 ], 'LineWidth', line_width );
    if algo.do_GammaGlobalQGlobal_A
        plot( ha3B, [ -1, 1 ]*t0, [ 1, 1]*Q_0A*m_ave, '-', 'Color', [ 0.6 0.6 0.6 ], 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaGlobalQGlobal_AB
        plot( ha3B, [ -1, 1 ]*t0, [ 1, 1]*Q_0AB*m_ave, '--', 'Color', [ 0.6 0.6 0.6 ], 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaGlobalQGlobal_ABLOS
        plot( ha3B, [ -1, 1 ]*t0, [ 1, 1]*Q_0ABLOS*m_ave, '.-', 'Color', [ 0.6 0.6 0.6 ], 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaLocalQLocal_Time_A
        plot( ha3B, t_ls, Q_1A*m_ave, 'o-', 'Color', [ 0 0 1 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaLocalQLocal_Time_AB
        plot( ha3B, t_ls, Q_1AB*m_ave, 'd-', 'Color', [ 0 0 1 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaLocalQLocal_Time_ABLOS
        plot( ha3B, t_ls, Q_1ABLOS*m_ave, 's-', 'Color', [ 0 0 1 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaLocalQLocal_Theta_A
        plot( ha3B, t_ls, Q_4A*m_ave, 'o-', 'Color', [ 0.5 0 0.5 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaLocalQLocal_Theta_AB
        plot( ha3B, t_ls, Q_4AB*m_ave, 'd-', 'Color', [ 0.5 0 0.5 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaLocalQLocal_Theta_ABLOS
        plot( ha3B, t_ls, Q_4ABLOS*m_ave, 's-', 'Color', [ 0.5 0 0.5 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaGlobalQLocal_A
        plot( ha3B, t_ls, Q_2A*m_ave, 'o-', 'Color', [ 0 0.5 0 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaGlobalQLocal_AB
        plot( ha3B, t_ls, Q_2AB*m_ave, 'd-', 'Color', [ 0 0.5 0 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaGlobalQLocal_ABLOS
        plot( ha3B, t_ls, Q_2ABLOS*m_ave, 's-', 'Color', [ 0 0.5 0 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaGlobalQiGlobal_A
        plot( ha3B, t_ls, ipl_Qi_3A*m_ave, 'o-', 'Color', [ 1 0 0 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaGlobalQiGlobal_AB
        plot( ha3B, t_ls, ipl_Qi_3AB_A*m_ave, 'd-', 'Color', [ 1 0 0 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
        % plot( ha3B, t_ls, ipl_Qi_3AB_B*m_ave, 'd-', 'Color', [ 1 0 0 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    if algo.do_GammaGlobalQiGlobal_ABLOS
        plot( ha3B, t_ls, ipl_Qi_3ABLOS_A*m_ave, 's-', 'Color', [ 1 0 0 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
        % plot( ha3B, t_ls, ipl_Qi_3ABLOS_B*m_ave, 's-', 'Color', [ 1 0 0 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    set( ha3B, ...
        'YLim', [ 0.005, 2000 ]*Q_gas_kg, ...
        'YScale', 'log', ...
        'YTickMode', 'manual', ...
        'YTick', [ 10, 100, 1000, 10000, 100000, 1000000 ], ...
        'YTickLabel',{ '10^{1}', '10^{2}', '10^{3}', '10^{4}', '10^{5}', '10^{6}' } ...
    );
    hxl = get( ha3B, 'XLabel' );
    set( hxl, 'String', '$t$ [s]', 'Interpreter', 'latex', 'FontSize', font_size );
    hyl = get( ha3B, 'YLabel' );
    set( hyl, 'String', '$Q$ [kg/s$^{-1}$]', 'Interpreter', 'latex', 'FontSize', font_size );
end

function [ gamma, dgamma, Q, dQ ] = FindSolutionGammaGlobalQGlobal( debug, gamma0, Q0, r_comet, v_gas, a, alpha, m_ave, t, r, theta, n_e, dn_e, tol )
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
            'm_ave', m_ave, ...
            'a', a, ...
            'alpha', alpha, ...
            'r_comet', r_comet, ...
            'v_gas', v_gas, ...
            't', {t}, ...
            'r', {r}, ...
            'theta', {theta}, ...
            'n_e', {n_e}, ...
            'dn_e', {dn_e} ...
        );
    x0 = [ gamma0; log(Q0) ];
    sx = [ 0.1; 1 ];
    [ x, dx ] = ...
        expl_optimization( [], 'FindSolutionGammaGlobalQGlobal', x0, sx, tol, @TargetFunctionGammaGlobalQGlobal, data, 'enhanced', 'silent' );
    gamma = x(1);
    Q = exp(x(2));
    dgamma = dx(1);
    dQ = Q * dx(2);
end
    
function F = TargetFunctionGammaGlobalQGlobal( x, data )
    % Target function
    % Least-squares formulation given measurements at points
    gamma = x(1);
    if ( gamma < 0 ) 
        F = 1e20 * (1-gamma);
    elseif ( gamma > 1.5 )
        F = 1e20 * (gamma-1.5);
    else
        Q = exp(x(2));
        F = 0;
        n_e_model = cell( size(data.n_e) );
        for k = 1:data.nr_of_sets
            n_e_model{k} = ModelCometIonosphere( data.r{k}, data.theta{k}, data.r_comet, data.v_gas, data.a, data.alpha, gamma, Q );
            F = F + ...
                sum( ...
                    ( ( data.n_e{k} - n_e_model{k} ) ./ data.dn_e{k} ).^2 ...
                )/data.n(k);
        end
        F = F / data.nr_of_sets;
        if data.debug
            DebugPlot( data.nr_of_sets, data.t, data.n_e, n_e_model, sprintf( 'gamma %.4f Q %.4e', gamma, Q*data.m_ave ) );
        end
    end
end

function [ gamma, dgamma, Q, dQ ] = FindSolutionGammaGlobalQGlobalLOS( debug, gamma0, Q0, r_comet, v_gas, a, alpha, m_ave, t, r, theta, n_e, dn_e, tol )
    nr_of_sets = length(t);
    n = zeros( nr_of_sets, 1 );
    for k = 1:nr_of_sets
        n(k) = length(t{k});
    end
    data = ...
        struct( ...
            'debug', debug, ...
            'nr_of_sets', nr_of_sets, ...
            'n', n, ...
            'm_ave', m_ave, ...
            'a', a, ...
            'alpha', alpha, ...
            'r_comet', r_comet, ...
            'v_gas', v_gas, ...
            't', {t}, ...
            'r', {r}, ...
            'theta', {theta}, ...
            'n_e', {n_e}, ...
            'dn_e', {dn_e} ...
        );
    x0 = [ gamma0; log(Q0) ];
    sx = [ 0.1; 1 ];
    [ x, dx ] = ...
        expl_optimization( [], 'FindSolutionGammaGlobalQGlobalLOS', x0, sx, tol, @TargetFunctionGammaGlobalQGlobalLOS, data, 'enhanced', 'silent' );
    gamma = x(1);
    Q = exp(x(2));
    dgamma = dx(1);
    dQ = Q * dx(2);
end
    
function F = TargetFunctionGammaGlobalQGlobalLOS( x, data )
    % Target function
    % Least-squares formulation given measurements at points
    gamma = x(1);
    if ( gamma < 0 ) 
        F = 1e20 * (1-gamma);
    elseif ( gamma > 1.5 )
        F = 1e20 * (gamma-1.5);
    else
        N_sample = 10;
        Q = exp(x(2));
        F = 0;
        scale = linspace( 0, 1, N_sample )';
        n_e_model = cell( size(data.n_e) );
        for k = 1:data.nr_of_sets
            if k == 1
                n_e_model{k} = ...
                    ModelCometIonosphere( data.r{k}, data.theta{k}, data.r_comet, data.v_gas, data.a, data.alpha, gamma, Q );
            else
                n_e_LOS_model = zeros( data.n(k), 1 );
                the_r_A = data.r{1};
                the_theta_A = data.theta{1};
                the_r_B = data.r{k};
                the_theta_B = data.theta{k};
                for i = 1:data.n(k)
                    x_A = the_r_A(i) .* cosd(the_theta_A(i));
                    y_A = the_r_A(i) .* sind(the_theta_A(i));
                    x_B = the_r_B(i) .* cosd(the_theta_B(i));
                    y_B = the_r_B(i) .* sind(the_theta_B(i));
                    x_i = x_A + (x_B-x_A).*scale;
                    y_i = y_A + (y_B-y_A).*scale;
                    r_i = sqrt( x_i.^2 + y_i.^2 );
                    theta_i = atan2d( y_i, x_i );
                    d_AB = sqrt( ( x_A - x_B ).^2 + ( y_A - y_B ).^2 );
                    n_e_model_i = ModelCometIonosphere( r_i, theta_i, data.r_comet, data.v_gas, data.a, data.alpha, gamma, Q );
                    n_e_LOS_model(i) = d_AB * sum(n_e_model_i)/N_sample;
                end
                n_e_model{k} = n_e_LOS_model;
            end
            F = F + ...
                sum( ...
                    ( ( data.n_e{k} - n_e_model{k} ) ./ data.dn_e{k} ).^2 ...
                )/data.n(k);
        end
        F = F / data.nr_of_sets;
        if data.debug
            DebugPlotLOS( data.nr_of_sets, data.t, data.n_e, n_e_model, sprintf( 'gamma %.4f Q %.4e', gamma, Q*data.m_ave ) );
        end
    end
end

function [ gamma, dgamma, Q, dQ ] = FindSolutionGammaGlobalQLocal( debug, gamma0, Q0, t_ls, dtheta_homogeneity, r_comet, v_gas, a, alpha, m_ave, t, r, theta, n_e, dn_e, tol )
    % Optimize for gamma using optimization for Q as a subproblem
    nr_of_sets = length(t);
    n = zeros( nr_of_sets, 1 );
    for k = 1:nr_of_sets
        n(k) = length(t{k});
    end
    data = ...
        struct( ...
            'debug', debug, ...
            'nr_of_sets', nr_of_sets, ...
            'n', n, ...
            'm_ave', m_ave, ...
            'Q0', Q0, ...
            't_ls', t_ls, ...
            'dtheta_homogeneity', dtheta_homogeneity, ...
            'tolerance', tol, ...
            'a', a, ...
            'alpha', alpha, ...
            'r_comet', r_comet, ...
            'v_gas', v_gas, ...
            't', { t }, ...
            'r', { r }, ...
            'theta', { theta }, ...
            'n_e', { n_e }, ...
            'dn_e', { dn_e } ...
        );
    x0 = gamma0;
    sx = 0.1;
    [ x, dx ] = ...
        expl_optimization( [], 'FindSolutionGammaGlobalQLocal', x0, sx, tol, @TargetFunctionGammaGlobalQLocal, data, 'enhanced', 'silent' );
    gamma = x(1);
    dgamma = dx(1);
    % Then rerun for the subproblem with the optimal gamma
    [ Q, dQ, ~ ] = SolutionGammaGlobalQLocal( x, data );
end

function F = TargetFunctionGammaGlobalQLocal( x, data )
    % Target function
    % Least-squares formulation given measurements at points
    gamma = x(1);
    if ( gamma <= 0.5 ) 
        F = 1.e30 * ( 1 + 0.5 - gamma );
    elseif ( gamma > 1.5 )
        F = 1.e30 * ( 1 + gamma - 0.5 );
    else
        [ ~, ~, F ] = SolutionGammaGlobalQLocal( x, data );
        % Compute the target function value
        F = sum(F)/length(F);
    end
    % fprintf( 'Trying gamma %10.8f => F = %20.8f\n', gamma, F );
end

function [ Q, dQ, F ] = SolutionGammaGlobalQLocal( x, data )
    % Set up optimization of Q for the given gamma
    % set up an optimization for all local Q
    gamma = x(1);
    t_ls = data.t_ls;
    n_t_ls = length(t_ls);
    Q = zeros( n_t_ls, 1 );
    dQ = zeros( n_t_ls, 1 );
    F = zeros( n_t_ls, 1 );
    for i = 1:n_t_ls
        t_ref = t_ls(i);
        theta_ref = interp1( data.t{1}, data.theta{1}, t_ref );
        set_t = cell( data.nr_of_sets, 1 );
        set_r = cell( data.nr_of_sets, 1 );
        set_theta = cell( data.nr_of_sets, 1 );
        set_n_e = cell( data.nr_of_sets, 1 );
        set_dn_e = cell( data.nr_of_sets, 1 );
        for k = 1:data.nr_of_sets
            [ factor, selection ] = HomogeneitySelection( data.theta{k}, theta_ref, data.dtheta_homogeneity ); 
            set_t{k} = data.t{k}(selection);
            set_r{k} = data.r{k}(selection);
            set_theta{k} = data.theta{k}(selection);
            set_n_e{k} = data.n_e{k}(selection);
            set_dn_e{k} = data.dn_e{k}(selection) .* factor(selection);
        end
        [ Q(i), dQ(i), F(i) ] = ...
            FindSolutionGammaGlobalQLocal_sub( ...
                data.debug, ...
                gamma, data.Q0, ...
                data.r_comet, data.v_gas, ...
                data.a, data.alpha, data.m_ave, ...
                set_t, set_r, set_theta, set_n_e, set_dn_e, ...
                data.tolerance ...            
            );
    end
end

function [ Q, dQ, F ] = FindSolutionGammaGlobalQLocal_sub( debug, gamma, Q0, r_comet, v_gas, a, alpha, m_ave, t, r, theta, n_e, dn_e, tol )
    nr_of_sets = length(t);
    n = zeros( nr_of_sets, 1 );
    for k = 1:nr_of_sets
        n(k) = length(t{k});
    end
    data = ...
        struct( ...
            'debug', debug, ...
            'nr_of_sets', nr_of_sets, ...
            'n', n, ...
            'm_ave', m_ave, ...
            'gamma', gamma, ...
            'a', a, ...
            'alpha', alpha, ...
            'r_comet', r_comet, ...
            'v_gas', v_gas, ...
            't', { t }, ...
            'r', { r }, ...
            'theta', { theta }, ...
            'n_e', { n_e }, ...
            'dn_e', { dn_e } ...
        );
    x0 = log(Q0);
    sx = 1;
    [ x, dx ] = ...
        expl_optimization( [], 'FindSolutionGammaGlobalQLocal_sub', x0, sx, tol, @TargetFunctionGammaGlobalQLocal_sub, data, 'enhanced', 'silent' );
    Q = exp(x(1));
    dQ = Q * dx(1);
    % Get target function value at the optimum
    F = TargetFunctionGammaGlobalQLocal_sub( x, data );
end
    
function F = TargetFunctionGammaGlobalQLocal_sub( x, data )
    % Target function
    % Least-squares formulation given measurements at points
    Q = exp(x(1));
    F = 0;
    n_e_model = cell( data.nr_of_sets, 1 );
    for k = 1:data.nr_of_sets
        n_e_model{k} = ...
            ModelCometIonosphere( data.r{k}, data.theta{k}, data.r_comet, data.v_gas, data.a, data.alpha, data.gamma, Q );
        n_t = length(data.r{k});
        F = F + ...
            sum( ...
                ( ( data.n_e{k} - n_e_model{k} ) ./ data.dn_e{k} ).^2 ...
            )/n_t;
    end
    F = F/data.nr_of_sets;
    if data.debug
        DebugPlot( data.nr_of_sets, data.t, data.n_e, n_e_model, sprintf( 'gamma %.4f Q %.4e', data.gamma, Q*data.m_ave ) );
    end   
end

function [ gamma, dgamma, Q, dQ ] = FindSolutionGammaGlobalQLocalLOS( debug, gamma0, Q0, t_ls, dtheta_homogeneity, r_comet, v_gas, a, alpha, m_ave, t, r, theta, n_e, dn_e, tol )
    % Optimize for gamma using optimization for Q as a subproblem
    nr_of_sets = length(t);
    n = zeros( nr_of_sets, 1 );
    for k = 1:nr_of_sets
        n(k) = length(t{k});
    end
    data = ...
        struct( ...
            'debug', debug, ...
            'nr_of_sets', nr_of_sets, ...
            'n', n, ...
            'm_ave', m_ave, ...
            'Q0', Q0, ...
            't_ls', t_ls, ...
            'dtheta_homogeneity', dtheta_homogeneity, ...
            'tolerance', tol, ...
            'a', a, ...
            'alpha', alpha, ...
            'r_comet', r_comet, ...
            'v_gas', v_gas, ...
            't', { t }, ...
            'r', { r }, ...
            'theta', { theta }, ...
            'n_e', { n_e }, ...
            'dn_e', { dn_e } ...
        );
    x0 = gamma0;
    sx = 0.1;
    [ x, dx ] = ...
        expl_optimization( [], 'FindSolutionGammaGlobalQLocalLOS', x0, sx, tol, @TargetFunctionGammaGlobalQLocalLOS, data, 'enhanced', 'silent' );
    gamma = x(1);
    dgamma = dx(1);
    % Then rerun for the subproblem with the optimal gamma
    [ Q, dQ, ~ ] = SolutionGammaGlobalQLocalLOS( x, data );
end

function F = TargetFunctionGammaGlobalQLocalLOS( x, data )
    % Target function
    % Least-squares formulation given measurements at points
    gamma = x(1);
    if ( gamma <= 0.5 ) 
        F = 1.e30 * ( 1 + 0.5 - gamma );
    elseif ( gamma > 1.5 )
        F = 1.e30 * ( 1 + gamma - 0.5 );
    else
        [ ~, ~, F ] = SolutionGammaGlobalQLocalLOS( x, data );
        % Compute the target function value
        F = sum(F)/length(F);
    end
    % fprintf( 'Trying gamma %10.8f => F = %20.8f\n', gamma, F );
end

function [ Q, dQ, F ] = SolutionGammaGlobalQLocalLOS( x, data )
    % Set up optimization of Q for the given gamma
    % set up an optimization for all local Q
    gamma = x(1);
    t_ls = data.t_ls;
    n_t_ls = length(t_ls);
    Q = zeros( n_t_ls, 1 );
    dQ = zeros( n_t_ls, 1 );
    F = zeros( n_t_ls, 1 );
    for i = 1:n_t_ls
        t_ref = t_ls(i);
        theta_ref = interp1( data.t{1}, data.theta{1}, t_ref );
        set_t = cell( data.nr_of_sets, 1 );
        set_r = cell( data.nr_of_sets, 1 );
        set_theta = cell( data.nr_of_sets, 1 );
        set_n_e = cell( data.nr_of_sets, 1 );
        set_dn_e = cell( data.nr_of_sets, 1 );
        [ factor, selection ] = HomogeneitySelection( data.theta{1}, theta_ref, data.dtheta_homogeneity ); 
        for k = 1:data.nr_of_sets
            set_t{k} = data.t{k}(selection);
            set_r{k} = data.r{k}(selection);
            set_theta{k} = data.theta{k}(selection);
            set_n_e{k} = data.n_e{k}(selection);
            set_dn_e{k} = data.dn_e{k}(selection) .* factor(selection);
        end
        [ Q(i), dQ(i), F(i) ] = ...
            FindSolutionGammaGlobalQLocalLOS_sub( ...
                data.debug, ...
                gamma, data.Q0, ...
                data.r_comet, data.v_gas, ...
                data.a, data.alpha, data.m_ave, ...
                set_t, set_r, set_theta, set_n_e, set_dn_e, ...
                data.tolerance ...            
            );
    end
end

function [ Q, dQ, F ] = FindSolutionGammaGlobalQLocalLOS_sub( debug, gamma, Q0, r_comet, v_gas, a, alpha, m_ave, t, r, theta, n_e, dn_e, tol )
    nr_of_sets = length(t);
    n = zeros( nr_of_sets, 1 );
    for k = 1:nr_of_sets
        n(k) = length(t{k});
    end
    data = ...
        struct( ...
            'debug', debug, ...
            'nr_of_sets', nr_of_sets, ...
            'n', n, ...
            'm_ave', m_ave, ...
            'gamma', gamma, ...
            'a', a, ...
            'alpha', alpha, ...
            'r_comet', r_comet, ...
            'v_gas', v_gas, ...
            't', { t }, ...
            'r', { r }, ...
            'theta', { theta }, ...
            'n_e', { n_e }, ...
            'dn_e', { dn_e } ...
        );
    x0 = log(Q0);
    sx = 1;
    [ x, dx ] = ...
        expl_optimization( [], 'FindSolutionGammaGlobalQLocalLOS_sub', x0, sx, tol, @TargetFunctionGammaGlobalQLocalLOS_sub, data, 'enhanced', 'silent' );
    Q = exp(x(1));
    dQ = Q * dx(1);
    % Get target function value at the optimum
    F = TargetFunctionGammaGlobalQLocalLOS_sub( x, data );
end
    
function F = TargetFunctionGammaGlobalQLocalLOS_sub( x, data )
    % Target function
    % Least-squares formulation given measurements at points
    N_sample = 10;
    Q = exp(x(1));
    F = 0;
    scale = linspace( 0, 1, N_sample )';
    n_e_model = cell( size(data.n_e) );
    for k = 1:data.nr_of_sets
        if k == 1
            n_e_model{k} = ...
                ModelCometIonosphere( data.r{k}, data.theta{k}, data.r_comet, data.v_gas, data.a, data.alpha, data.gamma, Q );
        else
            n_e_LOS_model = zeros( data.n(k), 1 );
            the_r_A = data.r{1};
            the_theta_A = data.theta{1};
            the_r_B = data.r{k};
            the_theta_B = data.theta{k};
            for i = 1:data.n(k)
                x_A = the_r_A(i) .* cosd(the_theta_A(i));
                y_A = the_r_A(i) .* sind(the_theta_A(i));
                x_B = the_r_B(i) .* cosd(the_theta_B(i));
                y_B = the_r_B(i) .* sind(the_theta_B(i));
                x_i = x_A + (x_B-x_A).*scale;
                y_i = y_A + (y_B-y_A).*scale;
                r_i = sqrt( x_i.^2 + y_i.^2 );
                theta_i = atan2d( y_i, x_i );
                d_AB = sqrt( ( x_A - x_B ).^2 + ( y_A - y_B ).^2 );
                n_e_model_i = ModelCometIonosphere( r_i, theta_i, data.r_comet, data.v_gas, data.a, data.alpha, data.gamma, Q );
                n_e_LOS_model(i) = d_AB * sum(n_e_model_i)/N_sample;
            end
            n_e_model{k} = n_e_LOS_model;
        end
        F = F + ...
            sum( ...
                ( ( data.n_e{k} - n_e_model{k} ) ./ data.dn_e{k} ).^2 ...
            )/data.n(k);
    end
    F = F / data.nr_of_sets;
    if data.debug
        DebugPlotLOS( data.nr_of_sets, data.t, data.n_e, n_e_model, sprintf( 'gamma %.4f Q %.4e', data.gamma, Q*data.m_ave ) );
    end   
end

function [ gamma, dgamma, Qi, dQi ] = FindSolutionGammaGlobalQiGlobal( debug, gamma0, Q0, theta_i, r_comet, v_gas, a, alpha, m_ave, t, r, theta, n_e, dn_e, tol )
    nr_of_sets = length(t);
    n = zeros( nr_of_sets, 1 );
    for k = 1:nr_of_sets
        n(k) = length(t{k});
    end
    data = ...
        struct( ...
            'debug', debug, ...
            'nr_of_sets', nr_of_sets, ...
            'n', n, ...
            'theta_i', theta_i, ...
            'a', a, ...
            'alpha', alpha, ...
            'm_ave', m_ave, ...
            'r_comet', r_comet, ...
            'v_gas', v_gas, ...
            't', { t }, ...
            'r', { r }, ...
            'theta', { theta }, ...
            'n_e', { n_e }, ...
            'dn_e', { dn_e } ...
        );
    x0 = [ gamma0; log(Q0) ];
    sx = [ 0.1; ones(size(theta_i)) ];
    [ x, dx ] = ...
        expl_optimization( [], 'FindSolutionGammaGlobalQiGlobal', x0, sx, tol, @TargetFunctionGammaGlobalQiGlobal, data, 'random_light', 'silent' );
    gamma = x(1);
    Qi = exp(x(2:end));
    dgamma = dx(1);
    dQi = Qi .* dx(2:end);
end

function F = TargetFunctionGammaGlobalQiGlobal( x, data )
    % Target function
    % Least-squares formulation given measurements at points
    gamma = x(1);
    if ( gamma <= 0.5 ) 
        F = 1.e30 * ( 1 + 0.5 - gamma );
        F_smooth = -1;
    elseif ( gamma > 1.5 )
        F = 1.e30 * ( 1 + gamma - 0.5 );
        F_smooth = -1;
    else
        Qi = exp(x(2:end));
        F = 0;
        n_e_model = cell( data.nr_of_sets, 1 );
        for k = 1:data.nr_of_sets
            ipl_Qi = InterpolateQiFromThetaToTime( data.theta_i, Qi, data.theta{k} );
            n_e_model{k} = ...
                ModelCometIonosphere( data.r{k}, data.theta{k}, data.r_comet, data.v_gas, data.a, data.alpha, gamma, ipl_Qi );
            n_t = length(data.r{k});
            F = F + ...
                sum( ...
                    ( ( data.n_e{k} - n_e_model{k} ) ./ data.dn_e{k} ).^2 ...
                )/n_t;
        end
        F = F/data.nr_of_sets;
        full_Qi = log([ Qi; Qi(1) ]);
        diff_Qi = diff(full_Qi);
        F_smooth = sum(diff_Qi.^2)/length(Qi);
        F = F + 0.01*F_smooth;
        if data.debug
            DebugPlot( data.nr_of_sets, data.t, data.n_e, n_e_model, sprintf( 'gamma %.4f <Q> %.4e', gamma, mean(Qi)*data.m_ave ) );
        end
    end
    % fprintf( 'Trying gamma %10.8f and Qi => F = %20.8f %20.8f\n', gamma, F, F_smooth );
end

function [ gamma, dgamma, Qi, dQi ] = FindSolutionGammaGlobalQiGlobalLOS( debug, gamma0, Q0, theta_i, r_comet, v_gas, a, alpha, m_ave, t, r, theta, n_e, dn_e, tol )
    nr_of_sets = length(t);
    n = zeros( nr_of_sets, 1 );
    for k = 1:nr_of_sets
        n(k) = length(t{k});
    end
    data = ...
        struct( ...
            'debug', debug, ...
            'nr_of_sets', nr_of_sets, ...
            'n', n, ...
            'theta_i', theta_i, ...
            'a', a, ...
            'alpha', alpha, ...
            'm_ave', m_ave, ...
            'r_comet', r_comet, ...
            'v_gas', v_gas, ...
            't', { t }, ...
            'r', { r }, ...
            'theta', { theta }, ...
            'n_e', { n_e }, ...
            'dn_e', { dn_e } ...
        );
    x0 = [ gamma0; log(Q0) ];
    sx = [ 0.1; ones(size(theta_i)) ];
    [ x, dx ] = ...
        expl_optimization( [], 'FindSolutionGammaGlobalQiGlobalLOS', x0, sx, tol, @TargetFunctionGammaGlobalQiGlobalLOS, data, 'random_light', 'silent' );
    gamma = x(1);
    Qi = exp(x(2:end));
    dgamma = dx(1);
    dQi = Qi .* dx(2:end);
end

function F = TargetFunctionGammaGlobalQiGlobalLOS( x, data )
    % Target function
    % Least-squares formulation given measurements at points
    gamma = x(1);
    if ( gamma <= 0.5 ) 
        F = 1.e30 * ( 1 + 0.5 - gamma );
        F_smooth = -1;
    elseif ( gamma > 1.5 )
        F = 1.e30 * ( 1 + gamma - 0.5 );
        F_smooth = -1;
    else
        N_sample = 10;
        Qi = exp(x(2:end));
        F = 0;
        scale = linspace( 0, 1, N_sample )';
        n_e_model = cell( size(data.n_e) );
        for k = 1:data.nr_of_sets
            if k == 1
                ipl_Qi_A = InterpolateQiFromThetaToTime( data.theta_i, Qi, data.theta{k} );
                n_e_model{k} = ...
                    ModelCometIonosphere( data.r{k}, data.theta{k}, data.r_comet, data.v_gas, data.a, data.alpha, gamma, ipl_Qi_A );
            else
                ipl_Qi_B = InterpolateQiFromThetaToTime( data.theta_i, Qi, data.theta{k} );
                n_e_LOS_model = zeros( data.n(k), 1 );
                the_r_A = data.r{1};
                the_theta_A = data.theta{1};
                the_r_B = data.r{k};
                the_theta_B = data.theta{k};
                for i = 1:data.n(k)
                    x_A = the_r_A(i) .* cosd(the_theta_A(i));
                    y_A = the_r_A(i) .* sind(the_theta_A(i));
                    x_B = the_r_B(i) .* cosd(the_theta_B(i));
                    y_B = the_r_B(i) .* sind(the_theta_B(i));
                    x_i = x_A + (x_B-x_A).*scale;
                    y_i = y_A + (y_B-y_A).*scale;
                    r_i = sqrt( x_i.^2 + y_i.^2 );
                    theta_i = atan2d( y_i, x_i );
                    d_AB = sqrt( ( x_A - x_B ).^2 + ( y_A - y_B ).^2 );
                    the_ipl_Qi = ipl_Qi_A(i) + (ipl_Qi_B(i) - ipl_Qi_A(i)).*scale;
                    n_e_model_i = ModelCometIonosphere( r_i, theta_i, data.r_comet, data.v_gas, data.a, data.alpha, gamma, the_ipl_Qi );
                    n_e_LOS_model(i) = d_AB * sum(n_e_model_i)/N_sample;
                end
                n_e_model{k} = n_e_LOS_model;
            end
            F = F + ...
                sum( ...
                    ( ( data.n_e{k} - n_e_model{k} ) ./ data.dn_e{k} ).^2 ...
                )/data.n(k);
        end
        F = F/data.nr_of_sets;
        full_Qi = log([ Qi; Qi(1) ]);
        diff_Qi = diff(full_Qi);
        F_smooth = sum(diff_Qi.^2)/length(Qi);
        F = F + 0.01*F_smooth;
        if data.debug
            DebugPlotLOS( data.nr_of_sets, data.t, data.n_e, n_e_model, sprintf( 'gamma %.4f <Q> %.4e', gamma, mean(Qi)*data.m_ave ) );
        end
    end
    % fprintf( 'Trying gamma %10.8f and Qi => F = %20.8f %20.8f\n', gamma, F, F_smooth );
end

function ipl_Qi = InterpolateQiFromThetaToTime( theta_i, Qi, theta )
    theta_min = theta_i(1);
    theta_max = theta_min + 360;
    ipl_theta = theta;
    ipl_theta( theta < theta_min ) = ipl_theta( theta < theta_min ) + 360;
    ipl_theta( theta > theta_max ) = ipl_theta( theta > theta_max ) - 360;
    ipl_Qi = interp1( [ theta_i; theta_i(1)+360 ], [ Qi; Qi(1) ], ipl_theta );
end

function q_factor = JetEnhancement( theta, jet )
    % Unpack jet-related data
    f_jet = jet.f;
    theta_jet = -jet.theta;
    dtheta_jet = jet.dtheta;
    
    q_factor = ones(size(theta));
    
    n_jet = length(jet.f);
    for i=1:n_jet
        % Compute angular distance, but beware of the periodicity jump
        angular_separation = theta+180-theta_jet(i);
        angular_separation_alt = angular_separation-360;
        test = ( abs(angular_separation_alt) < abs(angular_separation) );
        angular_separation(test) = angular_separation_alt(test);

        % Compute jet density enhancement factor
        q_factor = q_factor + f_jet(i) * exp( -(angular_separation/dtheta_jet(i)).^2 );
    end
end

function q_gas = GasProduction( theta, r_comet, a, Q_gas, jet )
    
    % Compute the flux on the nucleus surface at the subsolar point
    q_factor = JetEnhancement( theta, jet );
    q_gas0 = ( 1 + a ) * q_factor .* Q_gas ./ ( 4 * pi * r_comet.^2 );
    
    % Compute the flux on the nucleus surface at the base of the streamline
    % as a function of theta
    q_gas = q_gas0 .* ( 1 + a*cosd(theta+180) ) / ( 1 + a );
end

function n_gas = GasComa( r, theta, r_comet, Q_gas, v_gas, a, jet )
    % Returns values of the gas coma neutral density
    % assuming radial expansion but a SZA modulated sublimation
    
    % Dimensionless distance
    rs = r / r_comet;
    
    % Compute the flux on the nucleus surface at the base of the streamline
    % as a function of theta
    q_gas = GasProduction( theta, r_comet, a, Q_gas, jet );
    
    % Obtain the neutral gas coma density at the desired location, i.e., as a
    % function of theta and of r
    n_gas = ( q_gas / v_gas ) ./ rs.^2;
end

function n_e = CometIonosphere( r, theta, r_comet, Q_gas, v_gas, a, alpha, gamma, jet ) 
    % Routine that computes total electron density as a function of theta and r

    % Gas production at the surface
    q_gas = GasProduction( theta, r_comet, a, Q_gas, jet );
    
    % Ionospheric charge density
    n_e = ...
        ( r_comet^(2*gamma) * alpha / ( (3-2*gamma)* v_gas^(gamma+1) ) ) ...
        * q_gas.^gamma .* ( r - r_comet ).^(3-2*gamma) ./ r.^2;
end

function q_gas = ModelGasProduction( theta, r_comet, a, Q )
    
    % Compute the flux on the nucleus surface at the subsolar point
    q_gas0 = ( 1 + a ) * Q ./ ( 4 * pi * r_comet.^2 );
    
    % Compute the flux on the nucleus surface at the base of the streamline
    % as a function of theta
    q_gas = q_gas0 .* ( 1 + a*cosd(theta+180) ) / ( 1 + a );
end

function n_gas = ModelGasComa( r, theta, r_comet, v_gas, a, Q )
    % Returns values of the gas coma neutral density
    % assuming radial expansion but a SZA modulated sublimation
    
    % Dimensionless distance
    rs = r / r_comet;
    
    % Compute the flux on the nucleus surface at the base of the streamline
    % as a function of theta
    q_gas = ModelGasProduction( theta, r_comet, a, Q );
    
    % Obtain the neutral gas coma density at the desired location, i.e., as a
    % function of theta and of r
    n_gas = ( q_gas / v_gas ) ./ rs.^2;
end

function n_e = ModelCometIonosphere( r, theta, r_comet, v_gas, a, alpha, gamma, Q ) 
    % Routine that computes total electron density as a function of theta and r
    % based on the model

    % Gas production at the surface
    q_gas = ModelGasProduction( theta, r_comet, a, Q );
    
    % Ionospheric charge density
    n_e = ...
        ( r_comet^(2*gamma) * alpha / ( (3-2*gamma) .* v_gas^(gamma+1) ) ) ...
        .* q_gas.^gamma .* ( r - r_comet ).^(3-2*gamma) ./ r.^2;
end

function [ r, theta ] = FlybyTrajectory( Vflyby, D_ca, theta_ca, t_ca, t )
    % Trajectory as a function of time
    
    % Cartesian coordinates
    x_sc = Vflyby * sind(theta_ca) * ( t - t_ca ) - D_ca * cosd( theta_ca );
    y_sc = Vflyby * cosd(theta_ca) * ( t - t_ca ) + D_ca * sind( theta_ca );
    
    % Polar coordinates
    r = sqrt( x_sc.^2 + y_sc.^2 );
    theta = atan2d( y_sc, x_sc );
    
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

function [ factor, selection ] = HomogeneitySelection( t, t_ref, dt_homogeneity )
    % Homogeneity factor
    factor = exp( ((t-t_ref)/dt_homogeneity).^2 );
    selection = ( factor < 1000 );
end
    
function [] = DebugPlot( nr_of_sets, t, n_e, n_e_model, txt )
    hf = findobj(0,'Tag','Debug');
    if isempty(hf)
        hf = figure( 'Tag', 'Debug' );
    else
        clf(hf);
    end
    ha = axes(hf);
    for k = 1:nr_of_sets
        semilogy( ha, t{k}, n_e{k}, 'b.', t{k}, n_e_model{k}, 'k-' );
        hold( ha, 'on' );
    end
    title( txt );
    pause( 0.001 );
end

function [] = DebugPlotLOS( nr_of_sets, t, n_e, n_e_model, txt )
    hf = findobj(0,'Tag','Debug');
    if isempty(hf)
        hf = figure( 'Tag', 'Debug' );
    else
        clf(hf);
    end
    ha1 = axes(hf, 'Units','normalized', 'Position', [ 0.15, 0.55, 0.80, 0.35 ] );
    ha2 = axes(hf, 'Units','normalized', 'Position', [ 0.15, 0.15, 0.80, 0.35 ] );
    for k = 1:nr_of_sets
        if k == 1
            semilogy( ha1, t{k}, n_e{k}, 'b.', t{k}, n_e_model{k}, 'k-' );
            hold( ha1, 'on' );
        else
            semilogy( ha2, t{k}, n_e{k}, 'b.', t{k}, n_e_model{k}, 'k-' );
            hold( ha1, 'on' );
        end
    end
    title( ha1, txt );
    pause( 0.001 );
end


    
    