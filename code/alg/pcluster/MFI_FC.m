function [exemplars, Z, P, silohouette_index] = MFI_FC(S, K, R, restarts);
% MFI_FC.m: A "multistart fast interchange" procedure for the
%        p-median problem.  The procedure is the Teitz and 
%        Bart (1968) VERTEX SUBSTITUTION HEURISTIC using the
%        using fast updates from Hansen and Mladenovic (1997), 
%        which are based on Whitaker's (1983) work.

% NOTE: This program differs from MFI.m in two respects:
%    1) This program reads in an n x n DISSIMILARITY matrix directly. 
%       Accordingly, this program is not restricted to Euclidean distances, 
%       and even accommodates ASYMMETRIC matrices.  To avoid "unusual" 
%       solutions, however, elements should be nonnegative.
%    2) This program reads a vector of 'costs' for objects to serve as
%       exemplars.  Thus, the 'FC' stands for the 'fixed charge' for
%       selecting an object as an exemplar. Objects with lower costs have
%       greater "preference" for serving as an exemplar.  If there are no 
%       differences among exemplars with respect to preference, we recommend
%       setting the n x 1 cost vector to all zeros.

% INPUTS:   S = an n x n matrix of DISSIMILARITIES for n objects
%           K = the desired number of medians/clusters/exemplars
%           R = an n x 1 vector of exemplar preferences
%           restarts = the desired number of restarts
% OUTPUTS:  exemplars = the best found set of medians/exemplars
%           Z = the best-found objective function value (including the
%               fixed costs of the exemplars).
%           silohoette_index = index from Kaufman & Rousseeuw (1990)
%           P an n by 1 vector of cluster assignments

state = 1;                       % fix state for testing
rand('state', state);
tic; 
[n,n1] = size(S);

nreps = restarts; gbest = 9.9e+12; fstore = zeros(nreps,1); 
c1 = zeros(n,1); c2 = zeros(n,1); p = K;

for klk = 1:nreps
    s = randperm(n);                        % Randomly order the exemplar
    a1 = S(:,s(1:p));                       % candidates and let the first
    fbest = sum(R(s(1:p)));                 % p be the selected subset
    for i = 1:n;                            % Note: In this program
        dmin = 9.9e+12;                     % i indexes the rows, which
        for j = 1:p                         % assigned to exemplars or
            if a1(i,j) < dmin               % selected columns indexed by j
                dmin = a1(i,j); jsel = j;
            end
        end
        fbest = fbest + dmin; c1(i) = s(jsel); % fbest = initial objective; 
    end                                        % c1(i) = exemplar to which
    for i = 1:n                                % object i is most similar
        dmin = 9.9e+12;                        % c2(i) = exemplar to which
        for j = 1:p                            % object is 2nd most similar
            jj = s(j);                         % is computed in this block
            if c1(i) == jj
                continue
            end
            if S(i,jj) < dmin
                dmin = S(i,jj); jsel = jj;
            end
        end
        c2(i) = jsel;
    end
    trig = 0;
    while trig == 0   % The exchange process begins here
        wstar = 0;
% The block below finds, for each unseleted point (goin), the best exemplar
% to remove (goout) if adding goin.  The best (goin, goout) pair is
% identified and (wstar) is the resulting change in the objective function.
% If wstar >= 0, then the algorithm terminates because there is no viable
% exchange that will improve the objective function further.
        for goin = p+1:n 
            ii = s(goin); w = R(ii); u = zeros(n,1);
            for i = 1:n
                if S(i,ii) < S(i,c1(i))
                    w = w + S(i,ii) - S(i,c1(i));
                else
                    u(c1(i)) = u(c1(i)) + min(S(i,ii),S(i,c2(i))) - S(i,c1(i));
                end
            end
            dmin = 9.9e+12;
            for j = 1:p
                jj = s(j);
                if u(jj)-R(jj) < dmin
                    dmin = u(jj)-R(jj);
                    goout = j;
                end
            end
            w = w + dmin; 
            if w < wstar
                wstar = w; goinb = goin; gooutb = goout;
            end
        end
        if wstar >= -.00001
            trig = 1;
            continue
        end
% The block below updates the c1 and c2 vectors if an (goin, goout) swap
% results in an improvement.
        fbest = fbest + wstar; goinc = s(goinb); gooutc = s(gooutb);
        idum = s(goinb); s(goinb) = s(gooutb); s(gooutb) = idum;
        for i = 1:n
            if c1(i) == gooutc
                if S(i,goinc) <= S(i,c2(i))
                    c1(i) = goinc;
                else
                    c1(i) = c2(i);
                    dmin = 9.9e+12;
                    for j = 1:p
                        jj = s(j);
                        if c1(i) == jj
                            continue
                        end
                        if S(i,jj) < dmin
                            dmin = S(i,jj); jsel = jj;
                        end
                    end
                    c2(i) = jsel;
                end
            else
                if S(i,c1(i)) > S(i,goinc)
                    c2(i) = c1(i); c1(i) = goinc;
                elseif S(i,goinc) < S(i,c2(i))
                    c2(i) = goinc;
                elseif c2(i) == gooutc
                    dmin = 9.9e+12;
                    for j = 1:p
                        jj = s(j);
                        if c1(i) == jj
                            continue
                        end
                        if S(i,jj) < dmin
                            dmin = S(i,jj); jsel = jj;
                        end
                    end
                    c2(i) = jsel;
                end
            end
        end    
    end
    if fbest < gbest
        gbest = fbest; sbest = s; cbest = c1;
    end
    fstore(klk) = fbest;
end
% DOUBLE-CHECK THE BEST-FOUND SOLUTION
P = zeros(n,1);  exemplars = sbest(1:p); Z = sum(R(exemplars));
for i = 1:n
    dmin = 9.9e+12;
    for j = 1:p
        g = exemplars(j);
        if S(i,g) < dmin
            P(i) = j; dmin = S(i,g);
        end
    end
    Z = Z + dmin;
end
X = [];
TRI = zeros(n.*(n-1)./2,1);
ict = 0;
for i = 1:n-1
    for j = i+1:n
        ict = ict + 1; TRI(ict) = (S(i,j)+S(j,i))./2;
    end
end
u = silhouette(X, P, TRI');
W = 1:n;                           % This block is important for singletons
for k = 1:K                        % MATLAB's default is that objects in
    y = W(P==k); ylen = length(y); % singletons have silhoettes = 1,
    if ylen == 1                   % but Kaufman and Rousseeuw (1990, p. 85)
        u(y) = 0;                  % assume silhouettes = 0.
    end
end
silohouette_index = mean(u);
toc
