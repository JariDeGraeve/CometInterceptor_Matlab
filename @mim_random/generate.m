function a = generate( p, n )
% -----------------------------------------------------------------------------
%
% NAME
%
%	generate
%
% PURPOSE
%
%	Generate random numbers
%
% CALLING SEQUENCE
%
%	a = generate( p, n )
%
% INPUT PARAMETERS
%
%   p
%       a mim_random object
%
%   n
%       number of random numbers to generate
%
% OUTPUT PARAMETERS
%
%   a
%     a column vector with n random numbers
%
% -----------------------------------------------------------------------------
    
    % generate random number sequence
    switch p.type
        case 'uniform'
            x = rand( n, 1 );
        case 'normal'
            x = randn( n, 1 );
    end
    % transform
    a = p.alpha * x + p.beta;
    % counter
    p.counter = p.counter + n;

end

% =============================================================================
