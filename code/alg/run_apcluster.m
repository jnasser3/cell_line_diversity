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
%       verbose: Print logging to screen. Default true
%
%Outputs:
%       exemplars: cell array of the estimated exemplars. 
%       membership: Exemplar cluster membership

params = {'exclude',...
          'known',...
          'verbose'};
dflts = {'',...
    '',...
    true};
args = parse_args(params,dflts,varargin{:});

addpath('/cmap/projects/cell_line_diversity/code/pcluster');

%Subset ds if needed
exclude = get_array_input(args.exclude,'');
if ~isempty(exclude)
    objects = setdiff(ds.rid,exclude);
    ds = ds_slice(ds,'cid',objects,'rid',objects);
end

ds

%Set preferences vector
w = -1*mean([median(tri2vec(ds.mat)) max(tri2vec(ds.mat))])
known = get_array_input(args.known,'');
preferences = mk_known_pref(ds.rid,known,w)

%Run the clustering algo. Note that apcluster requires a similarity matrix,
%so change the sign on the distance matrix
[idx, netsim, dpsim, expref] = apcluster(-1*ds.mat,preferences);

exemplar_idx = unique(idx);
exemplars = ds.rid(exemplar_idx);
membership = ds.rid(idx);

end


function wt = mk_known_pref(objects,known,w)
%Given a subset of the objects which are known landmarks, returns a weight
%vector promoting the inclusion of these landmarks. %
%
%Known objects objects get a 'preference' of 0, unkown objects get a
%preference of w (which should be < 0).

n = numel(objects);
wt = w*ones(n,1);
lm_idx = ismember(objects,known);
wt(lm_idx) = 1;
end