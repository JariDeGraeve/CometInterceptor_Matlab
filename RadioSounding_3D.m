function [] = RadioSounding_3D()

    % Algorithm selection
%     algo = ...
%         struct( ...
%             'do_OneOverRgamma', true, ...
%             'do_OneOverRgamma_time', true, ...
%             'do_OneOverRgamma_sza', true ...
%         );
    
    algo = ...
        struct( ...
              'do_OneOverRgamma', true, ...
              'do_OneOverRgamma_time', true, ...
              'do_OneOverRgamma_sza', true ...
        );
    RadioSounding_Core( algo );
    
end

function [] = RadioSounding_Core( algo )
    
    %% Parameters
    % Debug script parameters
    asymmety_present = 1;
    jets_present = 1;

    % Physical parameters
    r_comet = 2000;            % comet radius, in m
    Q_gas_kg = 1000;           % gas production rate, in kg/s
    m_ave = 18 * 1.6e-27;      % average neutral gas particle mass (relative mass water * atomic mass unit), in kg/particle
    v_gas = 1000;              % neutral gas speed, in m/s
    D_ca_A = 1000000;          % distance at closest approach of A, in m
    t_ca_A = 0;                % time at closest approach of A, in s

%     delta_B1 = [-183000, 150000, -766000];   % coordinates spacecraft B1 relative to A, in m
%     delta_B2 = [-233000, -350000, -466000];  % coordinates spacecraft B2 relative to A, in m
    delta_B1 = [-170000, 270000, -300000];   % coordinates spacecraft B1 relative to A, in m
    delta_B2 = [50000, -390000, -200000];   % coordinates spacecraft B1 relative to A, in m
%     delta_B2 = [-627000, -577000, -127000];  % coordinates spacecraft B2 relative to A, in m

    V_flyby = [42000, 42000, 0];  % flyby speedvector, in m/s
    phi_ca = 30;               % azimuth angle at closest approach, in degrees
    theta_ca = 30;             % polar angle at closest approach, in degrees
    gamma = 0.8;               % ionization rate exponent, dimensionless
    alpha = 0.00002;           % ionization rate, in 1/s
    a = 0;                   % gas coma asymmetry factor, dimensionless
    if(asymmety_present)
        a = 0.8;
    end
    theta_jet = [];       % polar angle of center of gas jet, in degrees
    phi_jet = [];        % azimuth angle of center of gas jet, in degrees
    dangle_jet = [];     % angular half-width of gas jet, in degrees
    f_jet = [];             % density contrast in jet, dimensionless
    if(jets_present)
        theta_jet = [90,60, 80];       
        phi_jet = [50, 0, 50];        
        dangle_jet = [5, 15, 10];
        f_jet = [3,5,6];
    end
%    sc_inclination = 0;        % inclination of orbit
    rel_error = 0.05;          % relative measurement error, dimensionless
    t0 = 150;                  % limits of time interval, in s
    dt = 0.2;                  % measurement time resolution, in s
    
    % Numerical parameters
    t0 = 150;                  % limits of time interval, in s
    dt = 0.2;                  % measurement time resolution, in s
    dt_homogeneity = 5;        % time homogeneity scale, in s
    dsza_homogeneity = 1.0;  % theta homogeneity scale, in degrees
    tolerance = 1e-6;          % tolerance for solving optimization problem
    dt_ls = 5;                 % time interval for which to compute least-squares solutions
               
    % Plot parameters
    font_size = 14;
    marker_size = 8;
    line_width = 2;
    
    % Derive Q_gas in particles per second
    Q_gas = (Q_gas_kg/m_ave);
    n_gas_surface = Q_gas/(4*pi*r_comet^2*v_gas);
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

    rng(1);
    % Add lots of random little jets to randomize the environment
%     theta_rjets = (180).*rand(1,50);
%     phi_rjets = -100 + (100 + 100).*rand(1,50);
%     dangle_rjets = 3 + (6-3).*rand(1,50);
%     f_rjets = 0.1 + (1-0.1).*rand(1,50);

%     jet = struct( 'theta', [theta_jet,theta_rjets], 'phi', [phi_jet,phi_rjets], 'dangle', [dangle_jet,dangle_rjets] ,'f', [f_jet,f_rjets] );
    jet = struct( 'theta', theta_jet, 'phi', phi_jet, 'dangle', dangle_jet ,'f', f_jet );
%     % Conditions in ionosphere along trajectories
    f_A = CometIonosphere( r_A, phi_A, theta_A, r_comet, Q_gas, v_gas, a, jet );
    f_B1 = CometIonosphere( r_B1, phi_B1, theta_B1, r_comet, Q_gas, v_gas, a, jet );
    f_B2 = CometIonosphere( r_B2, phi_B2, theta_B2, r_comet, Q_gas, v_gas, a, jet );
    
    n_t = length(r_A);
    
    %% Generate measurements with noise
%     rng(1);
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
    plot_nn = CometIonosphere( plot_rr, plot_phi, plot_theta, r_comet, Q_gas, v_gas, a, jet );
%      plot_nn = CometReal_OneOverRgamma( plot_rr, plot_phi, plot_theta, Q_gas, gamma_gas );
    xslice = [];   
    yslice = 0;
    zslice = 0;
    slices = slice(ha2, -xx/1000, yy/1000, zz/1000,log10(plot_nn/1e6),xslice,yslice,zslice);
    set(slices,'EdgeColor','none', 'FaceColor', 'interp');

    
    sphere = xx.^2 + yy.^2 + zz.^2;
    p = patch(ha2,isosurface(-xx/1000, yy/1000, zz/1000, sphere, 1e11));
    isonormals(-xx/1000,yy/1000,zz/1000,sphere,p)
    isocolors(-xx/1000,yy/1000,zz/1000,log10(plot_nn/1e6),p)
    p.FaceColor = 'interp';
    p.EdgeColor = 'none';

    view(3)
    hold( ha2, 'on' );
    idx = round(n_t/2);
    plot3( ha2, ...
        r_A.*sind(theta_A).*cosd(phi_A)/1000, r_A.*sind(theta_A).*sind(phi_A)/1000, r_A.*cosd(theta_A)/1000, 'b-', ...
        r_B1.*sind(theta_B1).*cosd(phi_B1)/1000, r_B1.*sind(theta_B1).*sind(phi_B1)/1000, r_B1.*cosd(theta_B1)/1000,'g-', ...
        r_B2.*sind(theta_B2).*cosd(phi_B2)/1000, r_B2.*sind(theta_B2).*sind(phi_B2)/1000, r_B2.*cosd(theta_B2)/1000,'r-', ...
        0, 0, 0, 'ko', ...
        r_A(idx).*sind(theta_A(idx)).*cosd(phi_A(idx))/1000, r_A(idx).*sind(theta_A(idx)).*sind(phi_A(idx))/1000, r_A(idx).*cosd(theta_A(idx))/1000, 'bo', ...
        r_B1(idx).*sind(theta_B1(idx)).*cosd(phi_B1(idx))/1000, r_B1(idx).*sind(theta_B1(idx)).*sind(phi_B1(idx))/1000, r_B1(idx).*cosd(theta_B1(idx))/1000, 'go', ...
        r_B2(idx).*sind(theta_B2(idx)).*cosd(phi_B2(idx))/1000, r_B2(idx).*sind(theta_B2(idx)).*sind(phi_B2(idx))/1000, r_B2(idx).*cosd(theta_B2(idx))/1000, 'ro', ...
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
%                     1e-15, 2e-15, 3e-15, 4e-15, 5e-15, 6e-15, 7e-15, 8e-15, 9e-15 ...
%                     1e-14, 2e-14, 3e-14, 4e-14, 5e-14, 6e-14, 7e-14, 8e-14, 9e-14 ...
%                     1e-13, 2e-13, 3e-13, 4e-13, 5e-13, 6e-13, 7e-13, 8e-13, 9e-13 ...
%                     1e-12, 2e-12, 3e-12, 4e-12, 5e-12, 6e-12, 7e-12, 8e-12, 9e-12 ...
%                     1e-11, 2e-11, 3e-11, 4e-11, 5e-11, 6e-11, 7e-11, 8e-11, 9e-11 ...
%                     1e-10, 2e-10, 3e-10, 4e-10, 5e-10, 6e-10, 7e-10, 8e-10, 9e-10 ...
%                     1e-9, 2e-9, 3e-9, 4e-9, 5e-9, 6e-9, 7e-9, 8e-9, 9e-9 ...
%                     1e-8, 2e-8, 3e-8, 4e-8, 5e-8, 6e-8, 7e-8, 8e-8, 9e-8 ...
                    0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, ...
                    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, ...
                    1, 2, 3, 4, 5, 6, 7, 8, 9, ...
                    10, 20, 30, 40, 50, 60, 70, 80, 90, ...
                    100, 200, 300, 400, 500, 600, 700, 800, 900, ...
                    1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, ...
                    10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, ...
                    1e5, 2e5, 3e5, 4e5, 5e5, 6e5, 7e5, 8e5, 9e5, ...
                    1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6, 9e6, ...
                    1e7, 2e7, 3e7, 4e7, 5e7, 6e7, 7e7, 8e7, 9e7, ...
                    1e8, 2e8, 3e8, 4e8, 5e8, 6e8, 7e8, 8e8, 9e8, ...
                    1e9 ...
                ] ...
            ), ...
        'YTickLabel', ...
            { ...
%                 '10^{-15}', '', '', '', '', '', '', '', '', ...
%                 '10^{-14}', '', '', '', '', '', '', '', '', ...
%                 '10^{-13}', '', '', '', '', '', '', '', '', ...
%                 '10^{-12}', '', '', '', '', '', '', '', '', ...
%                 '10^{-11}', '', '', '', '', '', '', '', '', ...
%                 '10^{-10}', '', '', '', '', '', '', '', '', ...
%                 '10^{-9}', '', '', '', '', '', '', '', '', ...
%                 '10^{-8}', '', '', '', '', '', '', '', '', ...
                '', '', '', '', '', '', '', '', '', ...
                '10^{-1}', '', '', '', '', '', '', '', '', ...
                '10^0', '', '', '', '', '', '', '', '', ...
                '10^1', '', '', '', '', '', '', '', '', ...
                '10^2', '', '', '', '', '', '', '', '', ...
                '10^3', '', '', '', '', '', '', '', '', ...
                '10^4', '', '', '', '', '', '', '', '', ...
                '10^5', '', '', '', '', '', '', '', '', ...
                '10^6', '', '', '', '', '', '', '', '', ...
                '10^7', '', '', '', '', '', '', '', '', ...
                '10^8', '', '', '', '', '', '', '', '', ...
                '' ...
            }, ...
        'FontSize', font_size ...
    );

    drawnow;
         
    Q_start = Q_gas*1.33;
    gamma_start = gamma_gas*1.25;
    
    %% Problem 0:
    % Derive Q, gamma from measurements by A, B1 and B2 global
    if algo.do_OneOverRgamma
        set_t = { t; t; t};
        set_r = { r_A; r_B1; r_B2 };
        set_phi = { phi_A; phi_B1; phi_B2 };
        set_theta = { theta_A; theta_B1; theta_B2 };
        set_f = { f_A_obs; f_B1_obs; f_B2_obs };
        set_df = { f_A_obs * rel_error; f_B1_obs * rel_error;  f_B2_obs * rel_error };    
        [ gamma_0, dgamma_0, Q_0, dQ_0 ] = ...
            FindSolution_OneOverRgamma( ...
                false, Q_start, gamma_start, ...
                set_t, set_r, set_phi, set_theta, set_f, set_df, ...
                tolerance, v_gas ...
            );
        fprintf( 'Problem 0\n' );
        fprintf( '\n' );
        fprintf( '      Q                             gamma\n' );
        fprintf( '    : %12.6e +/- %12.6e  %10.6f +/- %10.6f\n', Q_0*m_ave, dQ_0*m_ave, gamma_0, dgamma_0 );
        fprintf( '\n' );
    end

    %% Problem 1:
    % Derive Q, gamma from measurements by A, B1 and B2 at a given time
    if algo.do_OneOverRgamma_time
        t_ls = ( -t0 : dt_ls : t0 )';
        n_t_ls = length(t_ls);
        gamma_1 = zeros( n_t_ls, 1 );
        dgamma_1 = zeros( n_t_ls, 1 );
        Q_1 = zeros( n_t_ls, 1 );
        dQ_1 = zeros( n_t_ls, 1 );
        fprintf( 'Problem 1\n' );
        fprintf( '\n' );
        fprintf( '      Q                             gamma\n' );
        for i = 1:n_t_ls
            t_ref = t_ls(i);
            [ factor, selection ] = HomogeneitySelection( t, t_ref, dt_homogeneity );           
            set_t = { t(selection); t(selection); t(selection) };
            set_r = { r_A(selection); r_B1(selection); r_B2(selection) };
            set_theta = { theta_A(selection); theta_B1(selection); theta_B2(selection) };
            set_phi = { phi_A(selection); phi_B1(selection); phi_B2(selection) };
            set_f = { f_A_obs(selection); f_B1_obs(selection); f_B2_obs(selection) };
            set_df = ...
                { ...
                    rel_error * f_A_obs(selection) .* factor(selection); ...
                    rel_error * f_B1_obs(selection) .* factor(selection);  ...
                    rel_error * f_B2_obs(selection) .* factor(selection)  ...
                };

            [ gamma_1(i), dgamma_1(i), Q_1(i), dQ_1(i) ] = ...
                FindSolution_OneOverRgamma( ...
                    false, Q_start, gamma_start, ...
                    set_t, set_r, set_phi, set_theta, set_f, set_df, ...
                    tolerance, v_gas ...
                );
            fprintf( '    : %12.6e +/- %12.6e  %10.6f +/- %10.6f\n', Q_1(i)*m_ave, dQ_1(i)*m_ave, gamma_1(i), dgamma_1(i) );
        end
        fprintf( '\n' );
    end

    %% Problem 2:
    % Derive Q, gamma from measurements by A, B1 and B2 at a given solar zenith angle
    if algo.do_OneOverRgamma_sza
        t_ls = ( -t0 : dt_ls : t0 )';
        n_t_ls = length(t_ls);
        gamma_2 = zeros( n_t_ls, 1 );
        dgamma_2 = zeros( n_t_ls, 1 );
        Q_2 = zeros( n_t_ls, 1 );
        dQ_2 = zeros( n_t_ls, 1 );
        fprintf( 'Problem 2\n' );
        fprintf( '\n' );
        fprintf( '      Q                             gamma\n' );
        sza_A = SolarZenithAngle(phi_A, theta_A);
        sza_B1 = SolarZenithAngle(phi_B1, theta_B1);
        sza_B2 = SolarZenithAngle(phi_B2, theta_B2);
        for i = 1:n_t_ls
            t_ref = t_ls(i);
            sza_ref = interp1( t, sza_A, t_ref );
            [ factor_A, selection_A ] = HomogeneitySelection( sza_A, sza_ref, dsza_homogeneity );           
            [ factor_B1, selection_B1 ] = HomogeneitySelection( sza_B1, sza_ref, dsza_homogeneity );           
            [ factor_B2, selection_B2 ] = HomogeneitySelection( sza_B2, sza_ref, dsza_homogeneity );           
            set_t = { t(selection_A); t(selection_B1); t(selection_B2) };
            set_r = { r_A(selection_A); r_B1(selection_B1); r_B2(selection_B2) };
            set_theta = { theta_A(selection_A); theta_B1(selection_B1); theta_B2(selection_B2) };
            set_phi = { phi_A(selection_A); phi_B1(selection_B1); phi_B2(selection_B2) };
            set_f = { f_A_obs(selection_A); f_B1_obs(selection_B1); f_B2_obs(selection_B2) };
            set_df = ...
                { ...
                    rel_error * f_A_obs(selection_A) .* factor_A(selection_A); ...
                    rel_error * f_B1_obs(selection_B1) .* factor_B1(selection_B1);  ...
                    rel_error * f_B2_obs(selection_B2) .* factor_B2(selection_B2)  ...
                };
            [ gamma_2(i), dgamma_2(i), Q_2(i), dQ_2(i) ] = ...
                FindSolution_OneOverRgamma( ...
                    false, Q_start, gamma_start, ...
                    set_t, set_r, set_phi, set_theta, set_f, set_df, ...
                    tolerance, v_gas ...
                );
            fprintf( '    : %12.6e +/- %12.6e  %10.6f +/- %10.6f\n', Q_2(i)*m_ave, dQ_2(i)*m_ave, gamma_2(i), dgamma_2(i) );
        end
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
    if algo.do_OneOverRgamma
        plot( ha3A, [ -1, 1 ]*t0, [ 1, 1]*gamma_0, '-', 'Color', [ 0.6 0.6 0.6 ], 'LineWidth', 0.5*line_width );
    end
    if algo.do_OneOverRgamma_time
        plot( ha3A, t_ls, gamma_1, 'd-', 'Color', [ 0 0 1 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    if algo.do_OneOverRgamma_sza
        plot( ha3A, t_ls, gamma_2, 'd-', 'Color', [ 0.5 0 0.5 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
%     hxl = get( ha3A, 'XLabel' );
%     set( hxl, 'String', '$t$ [s]', 'Interpreter', 'latex', 'FontSize', font_size );
    hyl = get( ha3A, 'YLabel' );
    set( hyl, 'String', '$\gamma$', 'Interpreter', 'latex', 'FontSize', font_size );
    set( ha3A, 'YLim', [ 0.5, 1.5 ]*gamma_gas );

    set( ha3B, ...
        'YLim', [ 1e-6, 1e6 ]*Q_gas, ...
        'YScale', 'log', ...
        'YTickMode', 'manual', ...
        'YTick', [  0.000000001, 0.00000001, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, ...
                     10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000, ...
%                    1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, ...
%                    1e20, 1e21, 1e22, 1e23, 1e24, 1e25, 1e26, 1e27, 1e28, 1e29 ...
                   ], ...
        'YTickLabel',{ '10^{-9}', '10^{-8}', '10^{-7}', '10^{-6}', '10^{-5}', '10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}', '10^{0}', ...
                         '10^{1}', '10^{2}', '10^{3}', '10^{4}', '10^{5}', '10^{6}', '10^{7}', '10^{8}', '10^{9}', ...
%                        '10^{10}', '10^{11}', '10^{12}', '10^{13}', '10^{14}', '10^{15}', '10^{16}', '10^{17}', '10^{18}', '10^{19}', ...
%                        '10^{20}', '10^{21}', '10^{22}', '10^{23}', '10^{24}', '10^{25}', '10^{26}', '10^{27}', '10^{28}', '10^{29}' ...
                       });
    hold( ha3B, 'on' );
    plot( ha3B, t, Q_gas*m_ave*ones(size(t)), '-', 'Color', [ 0 0 0 ], 'LineWidth', line_width );
    if algo.do_OneOverRgamma
        plot( ha3B, [ -1, 1 ]*t0, [ 1, 1]*Q_0*m_ave, '-', 'Color', [ 0.6 0.6 0.6 ], 'LineWidth', 0.5*line_width );
    end
    if algo.do_OneOverRgamma_time
        plot( ha3B, t_ls, Q_1*m_ave, 'd-', 'Color', [ 0 0 1 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    if algo.do_OneOverRgamma_sza
        plot( ha3B, t_ls, Q_2*m_ave, 'd-', 'Color', [ 0.5 0 0.5 ], 'MarkerSize', marker_size, 'LineWidth', 0.5*line_width );
    end
    hxl = get( ha3B, 'XLabel' );
    set( hxl, 'String', '$t$ [s]', 'Interpreter', 'latex', 'FontSize', font_size );
    hyl = get( ha3B, 'YLabel' );
    set( hyl, 'String', '$Q$ [kg/s$^{-1}$]', 'Interpreter', 'latex', 'FontSize', font_size );
%     set( ha3B, 'YLim', [ 1e-12, 10 ]*Q_gas );
end

function [ gamma, dgamma, Q, dQ ] = FindSolution_OneOverRgamma( debug, Q0, gamma0, t, r, phi, theta, f, df, tol, v_gas )
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
            f_model{k} = CometModel_OneOverRgamma( data.r{k}, data.phi{k}, data.theta{k}, Q, gamma, data.v_gas );
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
    q_gas = q_gas0 .* ( 1 + a*cosd(sza+180) ) / ( 1 + a ); 
end

function q_gas = GasProductionNormalized( phi, theta, r_comet, a, Q_gas, jet )
    
    % Unnormalized gas production
    q_gas0 = GasProduction(phi,theta,r_comet,a,Q_gas,jet);
    % Total gas production in model
    integr = integral2(@(phi,theta) GasProduction(phi,theta,r_comet,a,Q_gas,jet) .* sin(theta) .* r_comet.^2, 0, 2*pi, 0, pi);
    % Correct gas production so total gas production is Q_gas;
    q_gas = q_gas0 * (Q_gas/integr);

end

function n_n = CometIonosphere( r, phi, theta, r_comet, Q_gas, v_gas, a, jet ) 

    % Gas production at the surface
    q_gas = GasProductionNormalized( phi, theta, r_comet, a, Q_gas, jet );
    
    % Gas density at the surface
    n_n = q_gas ./ (v_gas * (r ./ r_comet) .^2);

end

function f = CometModel_OneOverRgamma( r, phi, theta, Q, gamma, v_gas ) 
    % Routine that computes the power law model
    % f = Q / r^gamma
    %f = (Q * 4 * pi *r_comet^2 *1000) ./ (r/r_comet).^(gamma); % Q is dichtheid aan opp op 1m straal (4*pi * v_gas)
    %f = (Q ./ r.^gamma);
    f = Q ./ (r.^gamma * 4 * pi * v_gas); 
    %4*pi*r_comet^2*v_gas)
    
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

function [ factor, selection ] = HomogeneitySelection( t, t_ref, dt_homogeneity )
    % Homogeneity factor
    factor = exp( ((t-t_ref)/dt_homogeneity).^2 );
    selection = ( factor < 1000 );
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


    
    