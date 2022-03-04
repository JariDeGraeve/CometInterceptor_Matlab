function p = mim_random( varargin )
% -----------------------------------------------------------------------------
%
% NAME
%
%	mim_random
%
% PURPOSE
%
%	Constructor for the mim_random class.
%
% CALLING SEQUENCE
%
%	p = mim_random( q )
%	p = mim_random( type, parameters ... )
%
% INPUT PARAMETERS
%
%   1st form (default constructor)
%
%     q
%       a mim_random object
%
%   2nd form
%
%     type
%       string, 'uniform', 'normal' ...
%
%     parameters
%       additional parameters
%       uniform : low, high (default 0, 1)
%       normal : mean, stdev (default 0, 1)
%
% OUTPUT PARAMETERS
%
%   p
%     a mim_random object
%
% -----------------------------------------------------------------------------

    persistent mim_random_id
    
    if isempty( mim_random_id )
        % remember initialization
        mim_random_id = 0;
        % seed for random number generator
        rng( 'shuffle' );
    end
    mim_random_id = mim_random_id + 1;
    
    if ( nargin == 1 ) && isa( varargin{1}, 'mim_random' )
        p = varargin{1};
    else
        switch varargin{1}
            case 'uniform'
                if nargin > 2
                    the_low = varargin{2};
                    the_high = varargin{3};
                else
                    the_low = 0;
                    the_high = 1;
                end
                alpha = the_high - the_low;
                beta = the_low;
            case 'normal'
                if nargin > 2
                    the_mean = varargin{2};
                    the_stdev = varargin{3};
                else
                    the_mean = 0;
                    the_stdev = 1;
                end
                alpha = the_stdev;
                beta = the_mean;
        end
        p = struct( ...
                'type', varargin{1}, ...
                'alpha', alpha, ...
                'beta', beta, ...
                'counter', 0, ...
                'id', mim_random_id ...
            );
        p = class( p, 'mim_random' );
    end

end

% =============================================================================
