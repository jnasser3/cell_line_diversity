function elements = compute_maximally_diverse(ds,p,varargin)
% elements = compute_maximally_diverse(ds,p,varargin)
%
% Given an input matrix of pairwise distances (aka dissimilarities),
% returns a set of elements that are maximally diverse.
%
% Computes the maximum elements by solving a quadratic program. If D is a
% matrix of pairwise distances, we seek a 0-1 vector x of support p such that 
%
% x^tDx + f^tx
%
% is maximized. The convex relaxation is to expand the search space to x
% being entrywise in [0,1] and ||x||_1 <= p. 
%
% The known elements are modeled through the parameter vector f. 
%
% Input
%       ds: (Required) A square struct of pairwise distances
%       p: (Required) The number of elements to return
%       exclude: Elements not to consider. Default ''
%       known: Id's that must be included in the set of maximally
%              diverse elements
%
% Output
%       elements: the id's of the most diverse elements

params = {'exclude',...
         'known'};
dflts = {'',...
         ''};
args = parse_args(params,dflts,varargin{:});

%Subset ds if needed
exclude = get_array_input(args.exclude,'');
if ~isempty(exclude)
    objects = setdiff(ds.rid,exclude);
    ds = ds_slice(ds,'cid',objects,'rid',objects);
end

n = numel(ds.rid);

%% Run the quadratic program.
%Define Objective function. quadprog works on minimization problems. Since
%maximizing diversity is equivalent to minimizing similarity, convert out
%distance matrix into a similarity matrix by flipping the sign.
H = double(-1*ds.mat);

% Take into account preference for the lines through the f parameter.
known = get_array_input(args.known,'');
known_idx = ismember(ds.rid,known);
fudge_factor = 2*n*max(max(ds.mat));
f = double(-1*fudge_factor * known_idx);

%Bounds on x
lb = zeros(n,1);
ub = ones(n,1);

%Bounds on the L1 norm of x
k = p + numel(known);
A = ones(n,n);
b = k*ones(n,1);

x = quadprog(H,f,A,b,[],[],lb,ub);
[~,sort_idx] = sort(x,'descend');
elements = ds.rid(sort_idx(1:k));
end

