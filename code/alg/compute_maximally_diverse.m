function [elements, objective,exitflag,output,lambda] = compute_maximally_diverse(ds,p,varargin)
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
%       least_diverse: (NOT YET SUPPORTED). Boolean. If set to true, computes the least diverse
%              set of elements. Default false
%
% Output
%       elements: the id's of the most diverse elements. Inlcude's id's of
%                 known elements
%       objective: the value of the objective function

params = {'exclude',...
         'known',...
         'least_diverse'};
dflts = {'',...
         '',...
         false};
args = parse_args(params,dflts,varargin{:});

%Subset ds if needed
exclude = get_array_input(args.exclude,'');
if ~isempty(exclude)
    objects = setdiff(ds.rid,exclude);
    ds = ds_slice(ds,'cid',objects,'rid',objects);
end

known = get_array_input(args.known,'');
n = numel(ds.rid);

%% Run the quadratic program.
%Define Objective function. quadprog works on minimization problems. Since
%maximizing diversity is equivalent to minimizing similarity, convert our.
%distance matrix into a similarity matrix by flipping the sign. If we are
%indeed solving a minimization problem (finding least diverse elements)
%then we do not need to flip.

%at the moment only support most diverse.
sign_flipper = -1;

H = double(sign_flipper*ds.mat);

%Bounds on x
lb = zeros(n,1);
ub = ones(n,1);

%Bounds on the L1 norm of x
k = p + numel(known);
Aeq = ones(1,n);
b = k;

% Take into account preference for the elements through the f parameter.
f = set_preference_vector(ds,known);

options = optimoptions('quadprog','Display','iter');
[x, objective, exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,b,lb,ub,[],options);
[~,sort_idx] = sort(x,'descend');
elements = ds.rid(sort_idx(1:k));
end

function f = set_preference_vector(ds,known)
    %Since quadprog always solves a minimization problem, the sign of the
    %preference vector does not depend on whether we are trying to find the
    %most or least diverse lines
    
    n = numel(ds.rid);
    known_idx = ismember(ds.rid,known);
    fudge_factor = 2*n*max(max(ds.mat));
    f = double(-1*fudge_factor * known_idx);
end
