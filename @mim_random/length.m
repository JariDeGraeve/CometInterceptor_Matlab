function n = length( p )
% -----------------------------------------------------------------------------
%
% NAME
%
%	history
%
% PURPOSE
%
%	Get the length of the random generator sequence
%
% CALLING SEQUENCE
%
%	n = length( p )
%
% INPUT PARAMETERS
%
%   p
%     a mim_random object
%
% OUTPUT PARAMETERS
%
%   n
%     number of random numbers that have been produced
%
% -----------------------------------------------------------------------------

    n = p.counter;

end

% =============================================================================
