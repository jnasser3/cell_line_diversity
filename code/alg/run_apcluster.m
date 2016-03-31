function [exemplars,membership] = run_apcluster(ds,varargin)
%exemplars = run_apcluster(pdi,p,varargin)
%
%Returns approximate solution to the p-medians problem using affinity propogation
%Given pairwise
%distances between n objects, returns a set of p ( < n) objects that
%minimize the sum of the minimum distances of each object to it's nearest
%landmark point.
%
%
%Inputs:
%       ds: an n x n struct of pairwise distances
%       exclude: a cell array or .grp of objects to exclude. Default ''
%       known: a cell array or .grp of known exemplars. Default ''
%       pref_quantile: A number in [0,1]. Controls the preferences. Lower
%               means more exemplars. Eg: 0 means all points are exemplars. 
%               1 means only a single exemplar is chosen. Default .01
%       verbose: Print logging to screen. Default true
%
%Outputs:
%       exemplars: cell array of the estimated exemplars. 
%       membership: Exemplar cluster membership

params = {'exclude',...
          'known',...
          'pref_quantile',...
          'verbose'};
dflts = {'',...
    '',...
    .01,...
    true};
args = parse_args(params,dflts,varargin{:});

addpath('/cmap/projects/cell_line_diversity/code/pcluster');

%Subset ds if needed
exclude = get_array_input(args.exclude,'');
if ~isempty(exclude)
    objects = setdiff(ds.rid,exclude);
    ds = ds_slice(ds,'cid',objects,'rid',objects);
end

%Set preferences vector
% w = -1*mean([median(tri2vec(ds.mat)) max(tri2vec(ds.mat))]);
known = get_array_input(args.known,'');
preferences = mk_known_pref(-1*ds.mat,ds.rid,known,args.pref_quantile);

%Run the clustering algo. Note that apcluster requires a similarity matrix,
%so change the sign on the distance matrix
[idx, netsim, dpsim, expref] = apcluster(-1*ds.mat,preferences);

exemplar_idx = unique(idx);
exemplars = ds.rid(exemplar_idx);
membership = ds.rid(idx);

end


function pref = mk_known_pref(sim,objects,known,t)
%Given a subset of the objects which are known landmarks, returns a weight
%vector promoting the inclusion of these landmarks. %
%
%Known objects objects get a 'preference' of 1, unkown objects get a
%preference of w (which should be < 0).
n = numel(objects);

[pref_min, pref_max] = preferenceRange(sim)
const = t*pref_min + (1-t)*pref_max

pref = const*ones(n,1);
lm_idx = ismember(objects,known);
pref(lm_idx) = 1;
end