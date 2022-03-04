function [ varargout ] = BFGS_optimization( s, title, x0, sx, tol, fh, data, mode, varargin )
% -----------------------------------------------------------------------------
%
% NAME
%
%	semantics
%
% PURPOSE
%
%	This program performs BFGS multidimensional optimization
%   of a well-behaved problem.
%   This is a Quasi-Newton method based on an update of the Hessian
%   approximation using the BFGS approach. In fact, not the hessian is
%   stored, but its Cholesky decomposition upper triangular matrix.
%   This implementation computes the gradients and the hessian
%   numerically. To improve the conditioning, the problem is solved using a
%   scaling of the search space.
%   Optionally, the optimization can be constrained to a linear subspace.
%
% CALLING SEQUENCE
%
%	[ x ] = BFGS_optimization( s, title, x0, sx, tol, fh, data, mode )
%	[ x, dx ] = BFGS_optimization( s, title, x0, sx, tol, fh, data, mode )
%	[ x, dx, fcounter ] = BFGS_optimization( s, title, x0, sx, tol, fh, data, mode )
%	[ x, dx, fcounter, v ] = BFGS_optimization( s, title, x0, sx, tol, fh, data, mode )
%	[ x, dx, fcounter, v, R ] = BFGS_optimization( s, title, x0, sx, tol, fh, data, mode )
%	[ ... ] = BFGS_optimization( s, title, x0, sx, tol, fh, data, mode, R )
%	[ ... ] = BFGS_optimization( s, title, x0, sx, tol, fh, data, mode, A, b )
%	[ ... ] = BFGS_optimization( s, title, x0, sx, tol, fh, data, mode, R, A, b )
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
%   mode
%     string, indicates mode of operation
%       'silent' : no output produced
%       'monitor' : an optimization_dialog window shows the progress of the
%           optimization process graphically
%       'print' : an optimization_dialog_text window shows the progress of the
%           optimization process as text
%       'monitor_print' : both of the above
%       'debug' : both of the above + additional info
%     
%   R
%     initial estimate for the cholesky factor of the hessian
%     note that H = R'*R
%
%   A, b
%     constraint matrix and vector; the subspace is defined by A x = b
%
% OUTPUT PARAMETERS
%
%   x, dx
%     solution and error (dx = tol * sx)
%     note that dx may lose part of its meaning in the constrained case, as
%     the constraint will introduce certain correlations between different
%     error components
%
%   fcounter
%     counts the number of function evaluations that have been performed
%
%   v
%     target function value at x
%
%   R 
%     cholesky factor of the hessian of the target function at x
%     note that H = R'*R
%
% -----------------------------------------------------------------------------

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

    if show
        [ str_parameters, str_target_function, str_target_function_convergence ] = ...
            get_resource( ...
                s, 'string', 'mim', ...
                'parameters', 'target_function', 'target_function_convergence' ...
            );
        hs = ...
            optimization_dialog( ...
                s, 'create', tt, ...
                str_parameters, str_target_function, str_target_function_convergence, ...
                {}, ...
                { 'F' }, ...
                tol, 'background', 3 ...
            );
        hs = optimization_dialog( s, 'restart', hs, 1 );    
    end
    if print
        ht = ...
            optimization_dialog_text( ...
                s, 'create', tt, ...
                get_resource( s, 'string', 'mim', 'optimization_table_title' ), ...
                'background', 1 ...
            );
    end

    % Process input arguments and set up the computation
    % m denotes the dimensionality of x-space
    % n denotes the dimensionality of the search space
    % (for unconstrained problems, n == m; otherwise n < m)
    % For unconstrained problems, R refers to the non-scaled space
    % For constrained problems, R refers to the non-scaled subspace and
    % A and b refer to the original non-scaled space.
    switch length(varargin)
        case 0
            is_constrained = false;
            R = [];
        case 1
            is_constrained = false;
            R = varargin{1};
        case 2
            is_constrained = true;
            A = varargin{1};
            b = varargin{2};
            R = [];    
        case 3
            is_constrained = true;
            R = varargin{1};
            A = varargin{2};
            b = varargin{3};
    end
    % deal with the constraint matrices, scales, and the starting solution
    if is_constrained
        [ p, m ] = size(A);
        n = m - p;
        % scaling of the constraint matrix
        for i = 1:m
            A(:,i) = A(:,i) * sx(i);
        end
        U = null(A);        % basis for search space
        UU = U * U';
        x0 = nearest_in_subspace( x0./sx, A, b, UU );
                            % starting solution in the scaled search subspace
        sx_orig = sx;
        sx = sqrt( eig( U'*diag(sx_orig.^2)*U ) );
                            % equivalent scales in the non-scaled subspace
    else
        n = length(x0);
        % there is no constraint matrix
        % U = eye(n);       % basis for search space
        % UU = eye(n);
        x0 = x0 ./ sx;      % starting solution in the scaled search space
    end
    % scaling of the Hessian or its Cholesky decomposition
    if isempty(R)
        R = eye(n);
    else
        for i = 1:n    
            R(:,i) = R(:,i) * sx(i);
        end
    end
    
    % The problem will be solved in the scaled space or subspace for
    % an unconstrained or a constrained problem, respectively.

    % Auxiliary quantities
    zn = zeros(n,1);

    % Fixed parameters
    numerical_gradient_relative_step = 0.01;
    max_line_search_steps = 1000;
    tolerance_reduction_factor = exp(log(tol)/n);
        % ideally, you could arrive at the final result in n line search steps
    line_search_relative_step_threshold = exp(log(tol)/n-2);
        % solve line search steps to a precision that is a bit higher than
        % the tolerance reduction that you must achieve in each step
    interval_factor = 5;
    interval_factor_plus_1 = ...
        interval_factor + 1;
    interval_factor_plus_1_2 = ...
        interval_factor_plus_1.^2;
    interval_factor_over_interval_factor_plus_1 = ...
        interval_factor / interval_factor_plus_1;
        % interval_factor is the maximum ratio between the half-intervals
        % to be used in the parabola fit; it is also used to ensure a
        % consistent decrease in the length of the total interval that is
        % bracketing the minimum in a line search
    descent_direction_limiting_factor = 5;
    
    % initial point and function value
    iteration = 0;
    rel_tol = 1.0;
    if is_constrained
        x = zn; 
        v = fh( ( x0 + U * x ).*sx_orig, data );
    else
        x = x0; 
        v = fh( x.*sx, data );              
    end
    if show
        hs = optimization_dialog( s, 'set', hs, 1, v );
    end
    if print
        ht = ...
            optimization_dialog_text( ...
                s, 'set', ht, ...
                sprintf( '%06d  %9.3e  %20.15e', iteration, 1, v ) ...
            );
    end

    % initial gradient
    % can be computed as a one-sided gradient
    G = zn;
    step = numerical_gradient_relative_step * rel_tol;
    for ii = 1:n
        x_try = x;
        x_try(ii) = x_try(ii) + step;
        if is_constrained
            G(ii) = ( fh( ( x0 + U * x_try ).*sx_orig, data ) - v ) / step;
        else
            G(ii) = ( fh( x_try.*sx, data ) - v ) / step;
        end
    end

    % performance instrumentation
    % a counter keeps track of the number of function evaluations
    fcounter = 1 + n;

    % initial estimate of the hessian
    % rather than storing the hessian, we use its
    % cholesky decomposition upper triangular matrix
    %   R = chol(H), such that R'*R = H.
    if debug
        H = R'*R; 
    end
    HessianOK = false;
    
    % initial displacement
    % find descent step dx for line search
    % by solving H dx = -G
    % i.e.       R' R dx = -G
    dx = solve_for_dx( G, R, zn, n );

    % optimization loop
    while ( rel_tol > tol )
        
        % At this point you have (in the scaled space or subspace):
        % - a point x, the current best solution, with the corresponding value
        %   of the target function being v
        % - the gradient at that point G
        % - an estimate of the hessian at that point, either known explicitly
        %   as H (only in debugging mode) or as its cholesky factor R
        % - a step dx to be applied to find a second point along the
        %   descent direction, possibly closer to the minimum
        %   (the more so if the hessian is a good approximation)
        % - the norm of the latest gradient Gnorm
        % - the relative tolerance rel_tol, which starts at 1 and descends
        %   down to tol, or 0.5 * tol, as its smallest value

        % store previous iteration results
        iteration = iteration + 1;
        dxold = dx;
        xold = x;
        Gold = G;
        if debug
            Hold = H; % previous Hessian approximation
        end
        Rold = R;
        
        if debug
            fprintf( 'Iteration %6d\n', iteration );
            fprintf('x: ');
            fprintf( '%25.16f ', x );
            fprintf('\n');
            fprintf( 'v: %25.16f\n', v );
            fprintf('G: ');
            fprintf( '%25.16f ', G );
            fprintf('\n');
            fprintf('dx: ');
            fprintf( '%25.16f ', dx );
            fprintf('\n');
            
            if false 
                % change if you want to include the following plotting code
                % snippet, for exploring the environment
                ssx = rel_tol;
                ntt = 18;
                x_all = cell((2*ntt+1)*ones(1,2));
                v_all = cell((2*ntt+1)*ones(1,2));
                xtest = x;
                for i1 = -ntt:ntt
                    xtest(1) = x(1) + i1*ssx/9;
                    for i2 = -ntt:ntt
                        xtest(2) = x(2) + i2*ssx/9;
                            x_all{i1+ntt+1,i2+ntt+1} = xtest;
                            v_all{i1+ntt+1,i2+ntt+1} = fh( xtest.*sx, data );
                    end
                end
                cut1 = zeros(2*ntt+1,2*ntt+1);
                for i1 = -ntt:ntt
                    for i2 = -ntt:ntt
                        cut1(i1+ntt+1,i2+ntt+1) = v_all{i1+ntt+1,i2+ntt+1};
                    end
                end
                [ xx, yy ] = meshgrid( x(1)+(-ntt:ntt)*ssx/9, x(2)+(-ntt:ntt)*ssx/9 );
                clf;
                set( gca, 'XLim', [ x(1) - 10*ssx, x(1) + 10*ssx ], 'YLim',  [ x(2) - 10*ssx, x(2) + 10*ssx ] );
                hold on
                contour(xx, yy, cut1',50); 
                surface(xx,yy,cut1'-max(max(abs(cut1)))-100,'EdgeColor','none');
                plot( x(1), x(2), 'ko' );
                d_normalized = dx/norm(dx);
                plot( x(1)+ (0:0.05:1)*ssx*d_normalized(1), x(2)+(0:0.05:1)*ssx*d_normalized(2), 'b--');
                plot( x(1)+ (0:0.05:1)*ssx*G(1)/norm(G), x(2)+(0:0.05:1)*ssx*G(2)/norm(G), 'r--');
                plot( x(1)+ ssx, x(2), 'kd');
                plot( x(1), x(2)+ ssx, 'kd');
                grid on
                axis equal; axis tight;
                hold off
                % end exploring for test purposes
                draw_now
                waitforbuttonpress
            end

        end

        % Descent direction step (in the scaled space)
        norm_dx = norm(dx);
        if norm_dx < 1.e-16
            % step anomalously small
            % handle problematic case
            dx = tol * ones(n,1);
        elseif norm_dx < tol
            % step smaller than the specified tolerance;
            % make sure dx is large enough to avoid problems due to
            % roundoff; assume that the tol level is reasonable
            dx = ( tol / norm_dx ) * dx;
        elseif norm_dx > rel_tol * descent_direction_limiting_factor
            % step too large at this point;
            % for the sake of robustness, we limit the step size
            dx = ( rel_tol * descent_direction_limiting_factor / norm_dx ) * dx;
        % else
            % step dx is acceptable as it is
        end
        
        % Normalized descent direction vector
        d = dx / norm(dx);

        % Perform line search along that direction
        
        % First point
        % x, v are given

        % Second point, along the descent direction
        % Because of the computation of dx, and the modifications to it,
        % this point is neither too close, nor too far away.
        x_p = x + dx;
        if is_constrained
            v_p = fh( ( x0 + U * x_p ).*sx_orig, data );
        else
            v_p = fh( x_p.*sx, data );
        end
        fcounter = fcounter + 1;
        if debug
            fprintf( 'v_p %25.16f\n', v_p );
        end

        % Third point, along the descent direction
        % Note that, since we are not doing steepest descent, we do not know
        % the gradient along that direction.
        if v_p >= v
            % In principle, the direction should be a descent direction, so
            % the minimum should be situated between x and x_p.
            % In practice, the direction might not be really a descent direction
            % (e.g. due to roundoff problems).
            % A robust solution is to take a third point to the left (on the
            % opposite side of x_p) to try to bracket the minimum.
            % This is done while respecting the maximum size ratio between 
            % the left and right halves of the bracketing interval.
            ds = dx/interval_factor;
            x_try = x - ds;
            if is_constrained
                v_try = fh( ( x0 + U * x_try ).*sx_orig, data );
            else
                v_try = fh( x_try.*sx, data );
            end
            fcounter = fcounter + 1;
            % Test the third point: does it indeed give a higher value?
            while v_try <= v
                % The given direction was not a descent direction: in the
                % opposite sense you have found a better value. Shift the
                % interval until you bracket the minimum. This loop is
                % finite since ds is nonzero and interval_factor > 1. Note
                % that the maximum size ratio of the half-intervals is
                % respected while doing so.
                ds = ds * interval_factor;
                x_p = x;
                v_p = v;
                x = x_try;
                v = v_try;
                x_try = x - ds;
                if is_constrained
                    v_try = fh( ( x0 + U * x_try ).*sx_orig, data );
                else
                    v_try = fh( x_try.*sx, data );
                end
                fcounter = fcounter + 1;
            end
            % At this point v_try is higher, brackets the minimum.
            x_pp = x_try;
            v_pp = v_try;
        else % v_p < v
            % The given direction was clearly a descent direction.
            ds = dx * interval_factor;
            x_pp = x; x = x_p;
            v_pp = v; v = v_p;
            x_try = x + ds;
            if is_constrained
                v_try = fh( ( x0 + U * x_try ).*sx_orig, data );
            else
                v_try = fh( x_try.*sx, data );
            end
            fcounter = fcounter + 1;
            % Test the third point: does it still give a lower value? Or
            % has a higher value been found to bracket the minimum from the
            % right side?
            while v_try <= v
                %  Go as far as possible by shifting the interval until you 
                % bracket the minimum from both sides. This loop is
                % finite since ds is nonzero and interval_factor > 1. Note
                % that the maximum size ratio of the half-intervals is
                % respected while doing so.
                ds = ds * interval_factor;
                x_pp = x;
                v_pp = v;
                x = x_try;
                v = v_try;
                x_try = x + ds;
                if is_constrained
                    v_try = fh( ( x0 + U * x_try ).*sx_orig, data );
                else
                    v_try = fh( x_try.*sx, data );
                end
                fcounter = fcounter + 1;
            end
            % At this point v_try is higher, brackets the minimum.
            x_p = x_try;
            v_p = v_try;
        end
        
        % At this point, and viewing along d, you have the three points:
        %
        %    |
        %    |
        %    |        o v_pp
        %    |
        %    |                                           o v_p
        %    |
        %    |                        o v
        %    |
        %    |
        %    +---------+--------------+------------------+----------
        %            x_pp             x                 x_p
        %
        % All you know is that the middle point has the lowest v; nothing
        % is known in particular on whether v_pp or v_p is the highest.
        %
        % It is also guaranteed that the lengths x - x_pp and x_p - x do
        % not differ more than the given factor interval_factor.
        
        % Start line search by repeated improvement of the bracketing triple.
        % Improvements are made using three rules:
        % - we require that the length of the bracketing interval is decreasing
        %   in each step
        % - we try to locate the minimum by fitting a parabola and adopting the
        %   minimum of the parabola as the next best guess
        % The loop is stopped at the required precision or if a given maximum
        % number of steps is exceeded.

        line_search_initial_scale = ( x_p - x_pp )' * d;
        line_search_current_scale = line_search_initial_scale;
        icount = 0;
        while ...
            ( line_search_current_scale > ...
                max( [ line_search_initial_scale * line_search_relative_step_threshold, tol ] ) ...
            ) && ( icount < max_line_search_steps )

            icount = icount + 1;

            % Compute the abscissae for x_p and x_pp.
            % The abscissa for x is 0.
            abscis_p = ( x_p - x )' * d + 1.e-99;
            abscis_pp = ( x_pp - x )' * d - 1.e-99;
%             if abscis_p == 0
%                 if abscis_pp < 0
%                     abscis_p = 1.e-99;
%                 else
%                     abscis_p = -1.e-99;
%                 end
%             elseif abscis_pp == 0
%                 if abscis_p < 0
%                     abscis_pp = 1.e-99;
%                 else
%                     abscis_pp = -1.e-99;                
%                 end
%             end

            % As stated before, x is the central point with the current minimum v,
            % x_p is to the right, and x_pp to the left of it, so that
            % abscis_p > 0 and abscis_pp < 0 have opposite sign.
            %
            % An absolute requirement to be able to guarantee convergence 
            % is that the interval size should tend to zero. That requires 
            % that max( abscis_p, -abscis_pp ) is reduced. Moreover,
            % it must be reduced as sharply as possible.
            
            % There are three points available, from which we will attempt
            % to construct an interpolating parabola to find a new estimate
            % of the minimum. 
            %
            % To have a useful fit, the quantities abscis_p and abscis_pp 
            % should not be too different: We enforce a maximum ratio
            % between both
            %   1/interval_factor <= -abscis_p/abscis_pp <= interval_factor
            while -abscis_pp > interval_factor * abscis_p
                % The left half of the interval is too large.
                % We make the left half of the interval shorter,
                % while keeping the right half of the same length.
                % We break up the left half of the interval according to
                % the interval_factor itself (for efficiency reasons: 
                % It is then likely that only one step of the loop will be
                % needed).
                delta = 1.0001 * abscis_pp / interval_factor_plus_1; 
                xnew = x + delta * d;        
                if is_constrained
                    vnew = fh( ( x0 + U * xnew ).*sx_orig, data );
                else
                    vnew = fh( xnew.*sx, data );
                end
                fcounter = fcounter + 1;
                if vnew >= v
                    % We can simply bring x_pp closer. It is likely that
                    % the half interval ratio bounds are satisfied (when 
                    % the original violation of the bounds was not
                    % dramatic), but in general that is not absolutely 
                    % guaranteed. That is why this operation has to be 
                    % possibly repeated.
                    x_pp = xnew;
                    v_pp = vnew;
                    abscis_pp = ( x_pp - x )' * d - 1.e-99;
                else % vnew < v
                    % As vnew is the smallest of the three, you can now
                    % shift the bracketing interval. Because of our choice
                    % of xnew, this new bracketing interval satisfies the
                    % half interval ratio bounds, so that you will leave
                    % the loop after this step.
                    x_p = x;
                    v_p = v;
                    x = xnew;
                    v = vnew;
                    % x_pp, v_pp ok
                    abscis_p = ( x_p - x )' * d + 1.e-99;
                    abscis_pp = ( x_pp - x )' * d - 1.e-99;
                end
            end
            while abscis_p > - interval_factor * abscis_pp
                % The right half of the interval is too large.
                % Avoid getting stuck when abscis_pp is zero.
                % We make the right half of the interval shorter,
                % while keeping the left half of the same length.
                % We break up the right half of the interval according to
                % the interval_factor itself (for efficiency reasons: 
                % It is then likely that only one step of the loop will be
                % needed).
                delta = 1.0001 * abscis_p / interval_factor_plus_1; 
                xnew = x + delta * d;        
                if is_constrained
                    vnew = fh( ( x0 + U * xnew ).*sx_orig, data );
                else
                    vnew = fh( xnew.*sx, data );
                end
                fcounter = fcounter + 1;
                if vnew >= v
                    % We can simply bring x_p closer. It is likely that
                    % the half interval ratio bounds are satisfied (when 
                    % the original violation of the bounds was not
                    % dramatic), but in general that is not absolutely 
                    % guaranteed. That is why this operation has to be 
                    % possibly repeated.
                    x_p = xnew;
                    v_p = vnew;
                    abscis_p = ( x_p - x )' * d + 1.e-99;
                else % vnew < v
                    % As vnew is the smallest of the three, you can now
                    % shift the bracketing interval. Because of our choice
                    % of xnew, this new bracketing interval satisfies the
                    % half interval ratio bounds, so that you will leave
                    % the loop after this step.
                    x_pp = x;
                    v_pp = v;
                    x = xnew;
                    v = vnew;
                    % x_p, v_p ok
                    abscis_p = ( x_p - x )' *d + 1.e-99;
                    abscis_pp = ( x_pp - x )' * d - 1.e-99;
                end
            end

            % At this point, we have the three points ordered properly as
            % x_pp - x - x_p,
            % where the lengths of both half intervals are not too
            % disparate.

            % Fit a parabola through the three points
            %   p(t) = a0 + a1 ( t - x )+ a2 ( t - x )^2 / 2
            % where the fit is computed from interpolation.

            % Since the three given points surround the minimum, you
            % are sure that the parabola will give you a minimum that
            % is situated inside the interval.

            % The sign of a2 gives the curvature sense,
            % the parabola has a minimum if a2 > 0
            % and a maximum when a2 < 0;
            % given the nature of the set of the three points, we know
            % that a2 > 0, or, when the minimum is fairly flat, a2 == 0
            % the top of the parabola is located at
            %   d p(t) / dt = a1 + a2 ( t - x ) = 0
            % or
            %   t = x - a1 / a2
            alpha_p = abscis_pp * ( v_p - v );
            alpha_pp = abscis_p * ( v - v_pp );
            line_step = 0.5 * ( abscis_pp * alpha_p + abscis_p * alpha_pp ) ...
                                / ( alpha_p + alpha_pp - 1.e-99 );

            % Determine the point by redefining the position of
            % the top of the parabola by taking into account robustness
            % constraints and the requirement for
            % a contraction of the original interval.
            if line_step >= 0
                if line_step > abscis_p * interval_factor_over_interval_factor_plus_1
                    % This is only possible for fairly degenerate situations,
                    % where the top of the parabola more or less coincides with
                    % the right point; in which case we better do not rely on
                    % the parabola.
                    line_step = abscis_p * interval_factor_over_interval_factor_plus_1;
                elseif line_step < abscis_p / interval_factor_plus_1_2
                    % Avoid being too close to the middle point, but allow
                    % being close enough - thereby deliberately breaking the
                    % half interval requirement if needed, because of the
                    % confidence we have in the parabola fit.
                    line_step = abscis_p / interval_factor_plus_1_2;
                % else
                    % Keep the top of the parabola
                end
            else  
                if line_step < abscis_pp * interval_factor_over_interval_factor_plus_1
                    % This is only possible for fairly degenerate situations,
                    % where the top of the parabola more or less coincides with
                    % the left point; in which case we better do not rely on 
                    % the parabola.
                    line_step = abscis_pp * interval_factor_over_interval_factor_plus_1;
                elseif line_step < abscis_pp / interval_factor_plus_1_2
                    % Avoid being too close to the middle point, but allow
                    % being close enough - thereby deliberately breaking the
                    % half interval requirement if needed, because of the
                    % confidence we have in the parabola fit.
                    line_step = abscis_pp / interval_factor_plus_1_2;
                % else
                    % Keep the top of the parabola
                end
            end
            
            % The new point now certainly lies within the original
            % bracketing interval. In addition, it does not lie too close
            % to any of the three prior points, although it can be quite
            % close to the center (that is the goal of using the fit).

            % Apply the change to determine the new point, and evaluate the
            % target function there.
            xnew = x + line_step * d;                
            if is_constrained
                vnew = fh( ( x0 + U * xnew ).*sx_orig, data );
            else
                vnew = fh( xnew.*sx, data );
            end
            fcounter = fcounter + 1;

            % Construct the new bracketing set with three points.
            if vnew <= v
                % Make the interval strictly narrower by choosing the left
                % or the right sub-interval.
                if line_step < 0
                    % x_pp, v_pp : remain
                    x_p = x;
                    v_p = v;
                else % line_step > 0
                    % x_p, v_p : remain
                    x_pp = x;
                    v_pp = v;
                end
                x = xnew;
                v = vnew;
            else % vnew > v
                % Make the interval strictly narrower by replacing the
                % corresponding endpoint.
                if line_step < 0
                    % x_p, v_p : remain
                    x_pp = xnew;
                    v_pp = vnew;
                else % line_step > 0
                    % x_pp, v_pp : remain
                    x_p = xnew;
                    v_p = vnew;
                end
                % x, v : remain
            end

            % Line search scale that is achieved now
            line_search_current_scale = ( x_p - x_pp )' * d;
            
        end % while ...
            %       ( line_search_current_scale > ...
            %         max( [ line_search_initial_scale * line_search_relative_step_threshold, tol ] ) ...
            %       ) && ( icount < max_line_search_steps )

        % At this point, the line search has ended: It has produced a 
        % better x and v embedded in a narrower interval.

        % Determine the gradient at this point along this direction. It is
        % not really zero since we do not proceed until full convergence.
        % Use a central (though asymmetric) difference.
        abscis = line_search_current_scale; % that is, ( x_p - x_pp )' * d;        
        if abs(abscis) > 1.e-99
            G_along_line = ( v_p - v_pp ) / abscis;
        else
            % safe value to handle having a zero abscis value
            G_along_line = 0;
        end
        
        % Diagnostics
        if show
            hs = optimization_dialog( s, 'set', hs, rel_tol, v );
        end
        if print
            ht = ...
                optimization_dialog_text( ...
                    s, 'set', ht, ...
                    sprintf( '%06d  %9.3e  %20.15e', iteration, rel_tol, v ) ...
                );
        end
        
        change_in_solution = norm( x - xold );
        
        % You can accept the change in this direction
        % as the current value, since this was the steepest descent
        % direction (at the previous point).
        
        % Update the tolerance to use.
        % Note that the sum of all steps taken in the line search
        % might be bigger than rel_tol; to avoid infinite iteration,
        % we force rel_tol to decrease anyhow.
        rel_tol = min( [ max( [ rel_tol * tolerance_reduction_factor, 1.5*change_in_solution ] ), 0.95 * rel_tol ] );
                
        if rel_tol > tol
            
            % Further improvement will be necessary.
            
            % Gradient estimate.
            % In general, the gradient must be computed second-order
            % precise, so that you can expect to use it for constructing the
            % Hessian. This, however, is expensive, so do not do it
            % - in the first n steps when you probably are not yet in the
            %   quasi-newton regime.
            % - in the final steps as soon as the Hessian has been 
            %   determined sufficiently well.

            % You can exploit the fact that you just finished a line search,
            % so you know that the gradient along direction d
            % is the remaining gradient from the line search.
            % Do this by rotating to an appropriate frame, where the
            % first axis is along d, compute the gradient
            % in that frame (first component is zero), and rotate back.
            % For the case of n == 1, this doesn't work as G == 0 will be the
            % result, which is not quite true because the line search has not
            % been done exactly; we use a different approach then.
            OC = find_orthonormal_complement( d, 'full_set' );
            G = zn;
            step = numerical_gradient_relative_step * rel_tol;        
            if n == 1
                G(1) = G_along_line;
            % else
                % it is much more useful to consider G(1) == 0 rather than that
                % very crude approximation
            end
            if ( iteration > n-1 ) && ~HessianOK
                for ii = 2:n
                    x_try = x + step * OC(:,ii);
                    x_try2 = x - step * OC(:,ii);
                    if is_constrained
                        G(ii) = ( fh( ( x0 + U * x_try ).*sx_orig, data ) - fh( ( x0 + U * x_try2 ).*sx_orig, data ) ) / (2*step);
                    else
                        G(ii) = ( fh( x_try.*sx, data ) - fh( x_try2.*sx, data ) ) / (2*step);
                    end
                end
                fcounter = fcounter + 2*(n-1);
            else
                for ii = 2:n
                    x_try = x + step * OC(:,ii);
                    if is_constrained
                        G(ii) = ( fh( ( x0 + U * x_try ).*sx_orig, data ) - v ) / step;
                    else
                        G(ii) = ( fh( x_try.*sx, data ) - v ) / step;
                    end
                end
                fcounter = fcounter + (n-1);
            end
            G = OC * G;

            % Update the Hessian or its Cholesky decomposition using the 
            % BFGS approach; avoid division by zero.
            xx = x - xold;
            if norm(xx) > 0.1*tol
                % Update the hessian; if xx is too small, such an update 
                % will be badly conditioned; moreover, it will likely not 
                % add much.
                GG = G - Gold;
                
                % Standard BFGS update
                update = GG ./ sqrt(GG'*xx + 1.e-99);
                
%                 % Do not use the stuff below: not robust enough, although
%                 % it seems to work for perfect quadratic functions;
%                 % probably rather badly conditioned as one approaches the
%                 % optimum
%                 % alternative update
%                 % exploits the more detailed information along the search
%                 % line, in points v_p and v_pp
%                 GG_par = ( GG'*d ) * d;
%                 GG_prp = GG - GG_par;
%                 a_p = ( x_p - x )'*d;
%                 a_pp = ( x_pp - x )'*d;
%                 a_ppp = a_p - a_pp;
%                 a = xx'*d;
%                 dvp = v_p - v;
%                 dvpp = v_pp - v;
%                 GGppp = 2*( dvp/a_p - dvpp/a_pp)*d + ( a_ppp/a )*GG_prp;
%                 %update = GGppp ./ sqrt( GGppp'*(x_p-x_pp) +1.e-99 );
                
                % Standard downdate
                downdate = Gold ./ ( sqrt(-dxold'*Gold + 1.e-99) );
                
                % modify Hessian and/or its Cholesky factor
                R = cholupdate( Rold, update, '+' );
                [ R, failure ] = cholupdate( R, downdate, '-' );
                    % returning second argument avoids error signal;
                    % ignores downdate if problems occur
                    % avoid getting stuck by continuing the iteration
                    % because rel_tol has been decreased anyhow
                if debug
                    if failure > 0
                        disp('downdate failed');
                    end
                    H = Hold + update*update' - downdate*downdate';       
                end
                HessianOK = ( update'*update + downdate'*downdate < rel_tol );
            else
                HessianOK = true;
            end
            if debug
                disp('H using Cholesky:');
                disp(H); 
                disp('H straightforward:');
                disp(R'*R);
                disp('Hessian OK');
                disp(HessianOK);
            end

            % Find descent step dx for the next line search
            % by solving H dx = -G
            % i.e.       R' R dx = -G
            dx = solve_for_dx( G, R, zn, n );
            if debug
                disp('d using Cholesky:');
                disp(dx); 
                disp('d straightforward:');
                disp(-inv(H)*G);
            end
            
        end

    end % while optimizing

    % Return results and close open windows
    if show
        if ishandle(hs)
            optimization_dialog( s, 'delete', hs );
        end
    end
    if print
        if ishandle(ht)
            optimization_dialog_text( s, 'delete', ht );
        end
    end
    if debug
        disp( [ '# f evaluations : ', int2str(fcounter) ] );
    end
    
    % return results
    % go back to the non-scaled space
    if is_constrained
        x = ( x0 + U * x ) .* sx_orig;
        dx = tol * sx_orig;
    else
        x = x .* sx;
        dx = tol * sx;
    end

    switch nargout
        case 1
            varargout = { x };
        case 2
            varargout = { x, dx };
        case 3
            varargout = { x, dx, fcounter };
        case 4
            varargout = { x, dx, fcounter, v };
        case 5
            for i = 1:n    
                R(:,i) = R(:,i) / sx(i);
            end
            varargout = { x, dx, fcounter, v, R };
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


function dx = solve_for_dx( G, R, zn, n )
    % Find descent step dx for line search
    % by solving H dx = -G
    % i.e.       R' R dx = -G
    % Formally:
    %   dx = - inv(H) * G;
    % In practice:
    %   (a) forward step
    %       let rd = R * dx and gg = -G
    %       the system then is R' rd = gg
    rd = zn;
    rd(1) = -G(1)/R(1,1);
    for i = 1:n-1
        ip1 = i+1;
        G(ip1:n) = G(ip1:n) + R(i,ip1:n)' * rd(i);
        rd(ip1) = -G(ip1)/R(ip1,ip1);
    end
    %   (b) backward step
    %       the system now is R d = rd
    dx = zn;
    dx(n) = rd(n)/R(n,n);
    for i = n:-1:2
        im1 = i-1;
        rd(1:im1) = rd(1:im1) - R(1:im1,i) * dx(i);
        dx(im1) = rd(im1)/R(im1,im1);
    end
end


% =============================================================================
