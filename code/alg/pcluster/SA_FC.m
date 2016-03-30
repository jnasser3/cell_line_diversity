function [exemplars, Z, P, silohouette_index] = SA_FC(S, K, R, cool);
%SA_FC.m  -- Simulated Annealing Heuristic for the P-median Problem
%       We use the method published by Brusco & Koehn (2009),
%       except that this program operates on DISSIMILARITY matrix.

% INPUTS:   S = an n x n matrix of DISSIMILARITIES for n objects
%           K = the desired number of medians/clusters/exemplars
%           R = an n x 1 vector of exemplar preferences
%           cool = the simulated annealing parameter.  We recommend
%                  0.9, but 0.8 often suffices and is faster.
% OUTPUTS:  exemplars = the best found set of medians/exemplars
%           Z = the best-found objective function value
%           silohoette_index = index from Kaufman & Rousseeuw (1990)
%           P an n by 1 vector of cluster assignments

% NOTE: This program differs from SA.m in two respects:
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

% NOTE: There is no 'tmin' in this program -- stops when no inferior 
%       solution is accepted

state = 1;                 % fix state for testing
rand('state', state);tic;
[m,n] = size(S); p = K;
           
t = randperm(n);           % Initial random permutation
s = t(1:p);                % s = the current set of exemplars
x = t(p+1:n);              % u = unselected objects that are not exemplars
z = sum(R(s(1:p)));        % fixed cost of selected exemplars
z = z + sum(max(S(:,s),[],2)); % z = the current objective function value
zbest = z;
sbest = s;
pbest = p;                 % The starting solution
tempmax = 0;               % Find an initial temperature by sampling
for ijk = 1:200
  i3 = int32(rand.*(n-p)+1); 
  i = double(i3); 
  if i > n-p
    i = n-p;
  end
  goin = x(i);             % goin is the entering object
  j3 = int32(rand*p+1);
  goout = double(j3);         
  if goout > p                
    goout = p;             % goout is the exemplar position where goin enters
  end
  strial = s;
  strial(goout) = goin;
  ztrial = sum(max(S(:,strial),[],2));
  ztrial = ztrial + sum(R(strial(1:p)));
  if abs(z-ztrial) > tempmax
      tempmax = abs(z-ztrial);
  end
end

iloop = 10.*n;             % SA temperature length
temp = tempmax;            % SA initial temperature

c1 = zeros(m,1); c2 = zeros(m,1);
a1 = S(:,s);                       
    z = sum(R(s(1:p))); 
    for i = 1:m;                            % Note: In this program
        dmin = 9.9e+15;                    % i indexes the rows, which
        for j = 1:p                         % assigned to exemplars or
            if a1(i,j) < dmin               % selected columns indexed by j
                dmin = a1(i,j); jsel = j;
            end
        end
        z = z + dmin; c1(i) = s(jsel);      % z = initial objective; 
    end                                        % c1(i) = exemplar to which
    for i = 1:m                                % object i is most similar
        dmin = 9.9e+15;                       % c2(i) = exemplar to which
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
% THE SA ALGORITHM BEGINS ON THE NEXT LINE
terminate = 0;
while terminate == 0
    terminate = 1;
    for kkk = 1:iloop  
        goinidx = ceil(rand*(n-p));  
        goin = x(goinidx);
              
        w = R(goin); u = zeros(n,1);
        for i = 1:m
            if S(i,goin) < S(i,c1(i))
                    w = w + S(i,goin) - S(i,c1(i));
            else
                    u(c1(i)) = u(c1(i)) + min(S(i,goin),S(i,c2(i))) - S(i,c1(i));
            end
        end
        dmin = 9.9e+15;
        for j = 1:p
            jj = s(j);
            if u(jj)-R(jj) < dmin
                dmin = u(jj)-R(jj); gooutidx = j;
            end
        end
        w = w + dmin; goout = s(gooutidx); 
                     
        if w <= 0
            x(goinidx) = goout;
            s(gooutidx) = goin;
            z = z + w;
            for i = 1:m
                if c1(i) == goout
                    if S(i,goin) <= S(i,c2(i))
                        c1(i) = goin;
                    else
                        c1(i) = c2(i);
                        dmin = 9.9e+15;
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
                    if S(i,c1(i)) > S(i,goin)
                        c2(i) = c1(i); c1(i) = goin;
                    elseif S(i,goin) < S(i,c2(i))
                        c2(i) = goin;
                    elseif c2(i) == goout
                        dmin = 9.9e+15;
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
            if z < zbest
                zbest = z;
                sbest = s;
            end
        else
            s1 = rand;                             
            rcrit = exp(-w./temp);             
            if s1 <= rcrit
                terminate = 0;
                x(goinidx) = goout;
                s(gooutidx) = goin;
                z = z + w;
                for i = 1:m
                    if c1(i) == goout
                        if S(i,goin) <= S(i,c2(i))
                            c1(i) = goin;
                        else
                            c1(i) = c2(i);
                            dmin = 9.9e+15;
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
                        if S(i,c1(i)) > S(i,goin)
                            c2(i) = c1(i); c1(i) = goin;
                        elseif S(i,goin) < S(i,c2(i))
                            c2(i) = goin;
                        elseif c2(i) == goout
                            dmin = 9.9e+15;
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
        end
    end
    temp = temp .* cool;                 % Temperature reduction
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
