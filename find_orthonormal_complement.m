function [ B ] = find_orthonormal_complement( A, varargin )
% -----------------------------------------------------------------------------
%
% NAME
%
%	find_orthonormal_complement
%
% PURPOSE
%
%   This routine finds the orthonormal complement of a given subspace.
%   The subspace that is given is supposed to be in the form of a set
%   of orthonormal vectors (columns of A). The complement B will be 
%   returned as a set of orthonormal vectors (columns of B) as well.
%   Therefore [ A B ]' * [ A B ] = I
%
%   Note that this operation offers
%     B = find_orthonormal_complement( A, 'complement' )
%   as an alternative for
%     B = null( [ A, zeros(n,n-m) ] );
%   Its conditioning is somewhat worse, but it is faster.
%
% CALLING SEQUENCE
%
%	[ B ] = find_orthonormal_complement( A )
%	[ B ] = find_orthonormal_complement( A, 'complement' )
%	[ AB ] = find_orthonormal_complement( A, 'full_set' )
%	
% INPUT PARAMETERS
%   
%   A
%     system matrix
%     n x m, n >= m, m > 0
%
%   mode
%     string
%       'complement' : the routine returns the complement
%       'full_set' : the routine returns the full set of
%                    vectors forming a complete basis
%
% OUTPUT PARAMETERS
%
%   B or AB 
%     the set of column vectors defining the complement, n x (n-m), or 
%     the full set of columns vectors defining a complete basis, n x n
% 
% -----------------------------------------------------------------------------

    % Analyze argument list to check mode
    if isempty( varargin )
        mode = 'complement';
    else
        mode = varargin{1};
    end

    % Size of the matrix
    [ n, m ] = size( A );
                                % n is the dimension of the vector space
                                % m is the number of given column vectors 
                                % (assumed to form an orthonormal set)

    % Allocate storage for the full basis
    B = [ A, zeros(n,n-m) ];

    % Count how many orthonormal basis vectors we already have;
    % progressively add new orthonormal basis vectors
    candidates = eye(n,n-m);
    for k = m+1:n

        % Add the k-th basis vector

        % Find a vector that is "sufficiently perpendicular" to
        % the already selected set of k-1 basis vectors

        % First attempt:
        % Use e_{k-m} as a candidate vector,
        candidate = candidates(:,k-m);

        % Project onto the k-1 basis vectors you already have
        % Two alternative implementations of projection
        % 1. Loop implementation
    %     for j = 1:k-1
    %         B_j = B(:,j);
    %         candidate = candidate - B_j(k-m) * B_j;
    %     end
        % 2. Vector implementation
        % Note that the columns > k-1 of B contain all zeros,
        % so make full matrix-vector operations to avoid data
        % repacking
        projection = B(k-m,:);
        candidate = candidate - B * projection';
        % End of 2 alternative implementations

        % Check if there is something significant remaining
        the_norm2 = sum(candidate.*candidate);
        if the_norm2 < 0.1

            % Quality not good enough: make a number of new attempts
            % The choice of candidates is made in a deterministic manner
            % You are sure to have a good result in at most n steps
            candidate = [ 1:n ]'*(2/n);
            i_candidate = 1;
            while true

                % Project onto the k-1 basis vectors you already have
                candidate_prime = candidate';
                % Two alternative implementations of projection
                % 1. Loop implementation
                for j = 1:k-1
                    B_j = B(:,j);
                    candidate = candidate - ( candidate_prime * B_j ) * B_j;
                end
                % 2. Vector implementation
    %             projection = candidate_prime * B(:,1:k-1);
    %             candidate = candidate - B(:,1:k-1) * projection';
                % End of 2 alternative implementations

                % Check if there is something significant remaining
                the_norm2 = sum(candidate.*candidate);
                if the_norm2 > 0.1
                    break;
                else
                    % try another candidate
                    candidate = [ candidate(n); candidate(1:n-1) ];
                    i_candidate = i_candidate + 1;
                    if i_candidate > n
                        break; 
                        % avoid infinite loop, but you should normally not
                        % encounter this situation
                    end
                end

            end

        end

        % Use that processed candidate vector as the k-th basis vector.
        % That vector is indeed perpendicular to all previous
        % vectors since any projections along their directions have been removed
        % before. Renormalize it to unit length (its length is definitely > 0).
        B(:,k) = candidate/sqrt(the_norm2);

    end % adding basis vectors

    switch mode
        case 'complement'
            B = B(:,m+1:n);
        % otherwise
            % B = B;
    end

end

% =============================================================================