function [exemplars,membership] = run_apcluster(ds,varargin)
%exemplars = run_apcluster(pdi,p,varargin)
%
%Returns approximate solution to the p-medians problem using affinity propogation
%Given pairwise
%distances between n objects, returns a set of p ( < n) objects that
%minimize the sum of the minimum distances of each object to it's nearest
%landmark point.
%
%Inputs:
%       ds: an n x n struct of pairwise distances
%       include: a cell array or .grp of objects to include in the
%                consideration. Default all in ds.
%       exclude: a cell array or .grp of objects to exclude from consideration. Default ''
%       known: a cell array or .grp of known exemplars to cluster around. Default ''
%       preference: a cell array or .grp of exemplars we want to give preference to.
%       pref_quantile: A number in [0,1]. Controls the preferences. Lower
%               means more exemplars. Eg: 0 means all points are exemplars. 
%               1 means only a single exemplar is chosen. Default .01
%       pref_factor: Preference factor weight. A higher number means
%               preferences are given more weight. Default 1.5
%       verbose: Print logging to screen. Default true
%
%Outputs:
%       exemplars: cell array of the estimated exemplars. 
%       membership: Exemplar cluster membership

params = {'exclude',...
          'include',...
          'known',...
          'preference',...
          'pref_quantile',...
          'pref_factor',...
          'verbose'};
dflts = {'',...
    '',...
    '',...
    '',...
    .01,...
    1.5,...
    true};
args = parse_args(params,dflts,varargin{:});

addpath('/cmap/projects/cell_line_diversity/code/pcluster');

%Subset ds if needed
exclude = get_array_input(args.exclude,'');
if ~isempty(exclude)
    objects = setdiff(ds.rid,exclude);
    ds = ds_slice(ds,'cid',objects,'rid',objects);
end

%Subset ds if needed
include = get_array_input(args.include,ds.cid);
include = intersect(include,ds.cid);
ds = ds_slice(ds,'cid',include,'rid',include,'ingore_missing',true);

%Set preferences vector
known = get_array_input(args.known,'');
preference = get_array_input(args.preference,'');
initial_preference_weight = mk_known_pref(-1*ds.mat,...
    ds.rid,...
    known,...
    preference,...
    args.pref_quantile,...
    args.pref_factor);

%Run the clustering algo. Note that apcluster requires a similarity matrix,
%so change the sign on the distance matrix
[idx, ~, ~, ~] = apcluster(-1*ds.mat,initial_preference_weight);

exemplar_idx = unique(idx);
exemplars = ds.rid(exemplar_idx);
membership = ds.rid(idx);

end

function pref = mk_known_pref(sim,objects,known,candidate,t,pref_factor)
%Given a subset of the objects which are known landmarks, returns a weight
%vector promoting the inclusion of these landmarks. %
%
%Known objects objects get a 'preference' of 1, unkown objects get a
%preference controlling the total number of exemplars apcluster chooses.
n = numel(objects);

[pref_min, pref_max] = preferenceRange(sim);
const = t*pref_min + (1-t)*pref_max;

pref = pref_factor*const*ones(n,1);

candidate_idx = ismember(objects,candidate);
pref(candidate_idx) = const;

known_idx = ismember(objects,known);
pref(known_idx) = -1*const*100;


end