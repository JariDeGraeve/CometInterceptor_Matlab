function [ varargout ] = expl_optimization( s, title, x0, sx, tol, fh, data, exploration_type, mode, varargin )
% -----------------------------------------------------------------------------
%
% NAME
%
%	semantics
%
% PURPOSE
%
%	This program performs exploratory multidimensional optimization
%   of a well-behaved problem.
%   The routine works in several phases. In each phase, a number of 
%   attempts is made to find a (local) minimum using BFGS optimization
%   and a subsequent search of the solution environment
%   to assess whether this is a local or a global minimum. If a lower 
%   function value is obtained during a phase, a new phase is started
%   with a smaller domain step size. This is repeated until no better 
%   solution is found.
%   Optionally, the optimization can be constrained to a linear subspace.
%
% CALLING SEQUENCE
%
%	[ x ] = expl_optimization( s, title, x0, sx, tol, fh, data, exploration_type, mode )
%	[ x, dx ] = expl_optimization( s, title, x0, sx, tol, fh, data, exploration_type, mode )
%	[ ... ] = expl_optimization( s, title, x0, sx, tol, fh, data, exploration_type, mode, A, b )
%
% INPUT PARAMETERS
%   
%   s
%     a mim_session object, or [] if you want to run this
%     routine outside MIM
%
%   title
%     string saying what this optimization is all about, 
%     shown in diagnostic window title
%
%   x0
%     starting solution (n-dim column vector)
%
%   sx
%     characteristic scales
%
%   tol
%     tolerance on the solution, relative to sx
%
%   fh
%     a function handle
%
%   data
%     data structure transmitted to the fh function
%
%   exploration_type
%     string, indicates type of exploration used
%       'none' : no exploration
%       'standard' : exploration in each dimension
%       'enhanced' : exploration along all combinations of dimensions
%       'random_light' : random exploration in all dimensions, small nr of tries
%       'random_medium', 'random' : random exploration in all dimensions
%       'random_heavy' : random exploration in all dimensions, more tries
%
%   mode
%     string, indicates mode of operation
%       'silent' : no output produced
%       'monitor' : an optimization_dialog window shows the progress of the
%           optimization process graphically
%       'print' : an optimization_dialog_text window shows the progress of the
%           optimization process as text
%       'monitor_print' : both of the above
%       'debug' : both of the above + additional info
%     note that graphical monitoring of the exploratory step is 
%     currently supported only for the 1-dimensional case
%     
%   A, b
%     constraint matrix and vector; the subspace is defined by A x = b
%     if given, A must be non-empty
%
% OUTPUT PARAMETERS
%
%   x, dx
%     solution and error (dx = tol * sx)
%     note that dx may lose part of its meaning in the constrained case, as
%     the constraint will introduce certain correlations between different
%     error components
%
% -----------------------------------------------------------------------------

    % Create random number generator
    rn = mim_random( 'normal' );
    
    % Determine exploration type
    switch exploration_type
        case 'random_light'
            do_exploration = true;
            stochastic_exploration = true;
            stochastic_light = true;
            stochastic_heavy = false;
        case { 'random_medium', 'random' }
            do_exploration = true;
            stochastic_exploration = true;
            stochastic_light = false;
            stochastic_heavy = false;
        case 'random_heavy'
            do_exploration = true;
            stochastic_exploration = true;
            stochastic_light = false;
            stochastic_heavy = true;
        case 'standard'
            do_exploration = true;
            stochastic_exploration = false;
            deterministic_best = false;
        case 'enhanced'
            do_exploration = true;
            stochastic_exploration = false;
            deterministic_best = true;
        otherwise 
            do_exploration = false;
    end

    % Determine mode of operation
    switch mode
        case 'monitor'
            show = true;
            print = false;
            debug = false;
        case 'print'
            show = false;
            print = true;
            debug = false;
        case 'monitor_print'
            show = true;
            print = true;
            debug = false;
        case 'debug'
            show = true;
            print = true;
            debug = true;
        otherwise 
            show = false;
            print = false;
            debug = false;
    end

    if ~isa( s, 'mim_session' )
        % if the routine is called with a non-specified
        % session, ignore showing or printing
        show = false;
        print = false;
    end

    if print || show
        tt = [ 'mim - ', get( s, 'session' ), ' - ', title ];
    end

    % Algorithm parameters
    phase_tol = 0.001;
            % tolerance with which to solve the BFGS problem in
            % each phase; must be sufficiently small to allow BFGS to
            % really do its work
    step_fraction_max = 0.9;
    step_fraction_min = 0.5;
            % fraction by which the scale size is reduced when going to
            % a new phase
    max_phases = 1000;
            % maximum number of phases, i.e. repeated BFGS optimizations at
            % different places so as to avoid getting stuck in local
            % minima. Therefore, this is a limit on the number of local
            % minima that you want to consider. Note that when you are
            % hitting the precision limit because of the numerical
            % calculation of gradients and hessians, a high number of
            % apparent local minima exists - while this is typically not
            % the case for the actual target function that one wants to
            % consider.
    
    % Set up
    is_constrained = ~isempty(varargin);
    if is_constrained
        A = varargin{1};
        b = varargin{2};
        [ p, m ] = size(A);
        n = m - p;          % dimensionality of search space
        U = null(A);        % basis for search space
        UU = U * U';
        x0 = nearest_in_subspace( x0, A, b, UU );
                            % starting solution must be in the subspace
    else
        n = length(x0);     % dimensionality of search space
        % U = eye(n);       % basis for search space
        % UU = eye(n);
        % x0                % starting solution
    end
    d = length(x0);         % dimensionality of the unconstrained solution space
    
    % Initial Cholesky factor for the Hessian
    R = eye(n);
    
    % Perform phases until global optimum is believed to be found
    step_scale = 1;
    x_current = x0;
    if print
        the_text = ...
            mim_acculist( ...
                { [ sprintf( '%04d ', 0 ),  sprintf( '%04d ', 0 ), sprintf( '%25.16e ', x_current ) ] } ...
            );
    end
    i_phase = 0;    
    fcounter = 0;
    history_step = zeros(max_phases,1);
    history_f    = zeros(max_phases,1);
    log_tol = log(tol);

    if debug
        fprintf('Start expl_optimization\n');    
    end
    
    while ( step_scale > tol )
        
        % Begin new phase with reduced step
        i_phase = i_phase + 1;
        ssx = step_scale * sx;        
        
        % Do BFGS
        % at the present scale, and with a precision that goes only
        % a few orders of magnitude beyond the step size in the current
        % phase
        % reuse the hessian estimate obtained previously
        attempt_relative_precision = max( [ phase_tol, 0.1*tol/step_scale ] );
        attempt_absolute_precision = attempt_relative_precision * step_scale;
        
        if debug
            fprintf( '%04d ssx %.16f rel prec %.16f abs prec %.16f\n', i_phase, mean(ssx), attempt_relative_precision, attempt_absolute_precision );
        end

        [ x, dx, f_evaluations, F, R ] = ...
            BFGS_optimization( s, title, x_current, 0.1*ssx, attempt_relative_precision, fh, data, mode, R, varargin{:} );
        fcounter = fcounter + f_evaluations;
        
        if debug
            fprintf( '%25.16e ', x );
            fprintf( ': %25.16e\n', F );
        end
                    
        % prepare print info
        if print
            the_text = ...
                add( ...
                    the_text, ...
                    [ sprintf( '%04d ', i_phase ), ...
                        sprintf( '%25.16e ', x ), sprintf( ': %25.16e', F ) ] ...
                );
        end
        
        % Explore environment
        % determine step size, which is at least the
        % size of the change of solution in this phase
        delta_x = 1.5 * (x - x_current);
        delta = abs( delta_x );
        delta_too_small = ( delta < ssx );
        delta( delta_too_small ) = ssx( delta_too_small );
        
        found = false;
        test_scale = step_scale;
        while ~found && ( test_scale > attempt_absolute_precision )
            if do_exploration
                if stochastic_exploration
                    % perform stochastic sampling of the surroundings of
                    % the present point, looking for a better solution. The
                    % sampling is done in all directions, and with roughly
                    % the current step size.
                    if stochastic_light
                        n_random_tries = 2*d;
                    elseif stochastic_heavy
                        n_random_tries = 32*d;
                    else
                        n_random_tries = 8*d;
                    end
                    for i_random_tries = 1:n_random_tries
                        % generate random point
                        tau = generate( rn, d );
                        sigma = generate( rn, 1 );
                        tau = delta.*tau * ( abs(sigma).^0.25 /(norm(tau)+1.e-99) );
                        % try
                        xtry = x + tau;
                        if is_constrained
                            xtry = nearest_in_subspace( xtry, A, b, UU );
                        end
                        Ftry = fh( xtry, data );
                        fcounter = fcounter + 1;
                        if Ftry < F
                            found = true;
                            x_current = xtry;
                            break;
                        end
                    end
                else
                    % deterministic exploration
                    % explore at scale test_scale, corresponding to vector delta
                    if deterministic_best
                        % exploration along all combinations of dimensions
                        if is_constrained
                            [ x_current, F, fcounter, found ] = ...
                                try_up_to_dim( fh, data, x_current, x_current, F, delta, n, fcounter, is_constrained, A, b, UU );
                        else
                            [ x_current, F, fcounter, found ] = ...
                                try_up_to_dim( fh, data, x_current, x_current, F, delta, n, fcounter, is_constrained, [], [], [] );
                        end
                    else
                        % exploration only along the coordinate axes
                        for i = 1:d
                            % try left and right
                            xleft = x;
                            xleft(i) = x(i) - delta(i);
                            if is_constrained
                                xleft = nearest_in_subspace( xleft, A, b, UU );
                            end
                            Fleft = fh( xleft, data );
                            fcounter = fcounter + 1;
                            if Fleft < F
                                found = true;
                                x_current = xleft;
                                break;
                            end
                            xright = x;
                            xright(i) = x(i) + delta(i);
                            if is_constrained
                                xright = nearest_in_subspace( xright, A, b, UU );
                            end
                            Fright = fh( xright, data );
                            fcounter = fcounter + 1;
                            if Fright < F
                                found = true;
                                x_current = xright;
                                break;
                            end
                        end
                    end
                end
            else % no exploration
                % found = false;
            end

            % determine step_fraction to use
            step_fraction = ...
                step_fraction_max - ...
                ( step_fraction_max - step_fraction_min ) * log(test_scale) / log_tol;
            if found
                % search with larger scale, but make sure that test_scale
                % remains between tol and 1
                sqrt_step_fraction = sqrt(step_fraction);
                if test_scale <= sqrt_step_fraction
                    test_scale = test_scale / sqrt_step_fraction;
                    delta = delta / sqrt_step_fraction;
                % else
                    % continue with the present values
                end
            else
                % search with finer scale
                test_scale = test_scale * step_fraction;
                delta = delta * step_fraction;
            end
        end
        
        % At this point, either a better value was found, or you stopped
        % searching because you are approaching the precision of the
        % present BFGS solution.
            
        % Check whether this was a local or a global optimum
        if found
            % an improvement has been found
            % corresponding changes have already been made
            if debug
                fprintf( 'step %g\n', test_scale );
            end
        else
            % no improvement found
            x_current = x;
        end

        % modify the step scale
        % build in a safety feature that will force the step to go down
        % ultimately
        step_scale = min( [ test_scale, 100*step_fraction^(i_phase/5) ] );
        history_step(i_phase) = step_scale;
        history_f(i_phase) = F;        
        
        if i_phase > 1
            i_phase_compare = ceil(i_phase/3);
            if abs( history_step(i_phase) - history_step(i_phase_compare ) )/ ( i_phase - i_phase_compare ) < 0.01 * history_step(i_phase)
                % abort
                step_scale = tol/2;
            end
            if i_phase == max_phases - 1
                step_scale = tol/2;
            end
        end

        if debug
            clf
            subplot(2,1,1);
            semilogy( 1:i_phase, history_step(1:i_phase), 'ro' );
            subplot(2,1,2);
            semilogy( 1:i_phase, history_f(1:i_phase), 'ro' );
            waitforbuttonpress
        end

    end
                       
    % Display results if desired
    if print
        uiwait( text_dialog( s, tt, the_text ) );
    end
    
    % debugging info
    if debug
        disp( 'Final optimum' );
        disp( x_current );
        disp( 'Function evaluations' );
        disp( fcounter );
    end
    
    % Return results
    switch nargout
        case 1
            varargout = { x_current };
        case 2
            varargout = { x_current, ssx };
    end
    
end


function x = nearest_in_subspace( x, A, b, UU )
    % We assume that ~isempty(A)
    
    % The basic procedure goes as follows:
    %   find arbitrary point x0 in the subspace
    %   c = b - A * x;
    %   dx = A \ c;
    %   x0 = x + dx;
    %   project to subspace
    %   a = U' * ( x - x0 );
    %   project from subspace
    %   x = x0 + U * a;
    
    % Implementation
    % where UU = U * U'
    dx = A \ ( b - A * x );
    x0 = x + dx;
    x = x0 - UU * dx;

end


function [ x, F, f_counter, found ] = ...
                try_up_to_dim( fh, data, xtry, x, F, delta, m, f_counter, is_constrained, A, b, UU )
    % Recursive procedure to explore alternatives in the n-dimensional
    % space up to dimension m.
    % x is the solution vector (length n).
    % As soon as a better solution is found, stop the exploration and
    % restart BFGS.
    if m > 0
        % try center, left, and right
        xtry_here = xtry;
        [ x, F, f_counter, found ] = ...
            try_up_to_dim( fh, data, xtry_here, x, F, delta, m-1, f_counter, is_constrained, A, b, UU );
        if ~found
            xtry_here(m) = xtry(m) - delta(m);
            [ x, F, f_counter, found ] = ...
                try_up_to_dim( fh, data, xtry_here, x, F, delta, m-1, f_counter, is_constrained, A, b, UU );
            if ~found
                xtry_here(m) = xtry(m) + delta(m);
                [ x, F, f_counter ] = ...
                    try_up_to_dim( fh, data, xtry_here, x, F, delta, m-1, f_counter, is_constrained, A, b, UU );
            end
        end
    else
        % evaluate for the given xtry
        if is_constrained
            xtry = nearest_in_subspace( xtry, A, b, UU );
        end
        if any( xtry-x ~= 0. )
            Ftry = fh( xtry, data );
            f_counter = f_counter + 1;        
            if Ftry < F
                found = true;
                x = xtry;
                F = Ftry;
            else
                found = false;
                % x, F remain the optimum
            end
        else
            found = false;
            % x, F remain the optimum
        end            
    end
end


% =============================================================================