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
%       num_exemplar: A positive integer specifying the exact number of
%           exemplars to return (including known exemplars). Set equal to 0
%           to just use preferences and let affinity propogation choose
%           optimal number. Default 0. If not 0, should always be greater
%           than number of elements in known.
%       preference: a cell array or .grp of exemplars we want to give preference to.
%       pref_quantile: A number typically in [0,1]. Controls the preferences. Lower
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
          'num_exemplar',...
          'pref_quantile',...
          'pref_factor',...
          'verbose'};
dflts = {'',...
    '',...
    '',...
    '',...
    0,...
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

%Get array variables
known = get_array_input(args.known,'');
preference = get_array_input(args.preference,'');

%Make sure num_exemplar is valid.
if args.num_exemplar
    assert(isint(args.num_exemplar),'Number of exemplars must be a positive integer')
    assert(args.num_exemplar > 0,'Number of exemplars must be a positive integer')
    assert(args.num_exemplar >= numel(unique(known)), 'Number of exemplars should be greater than known!')
end

%Run ap using a one-shot approach or bisection
pref_quantile = args.pref_quantile;
if ~args.num_exemplar
    %%One Shot
    %Set preferences vector
    preference_vector = mk_preference_vector(-1*ds.mat,...
        ds.rid,...
        known,...
        preference,...
        pref_quantile,...
        args.pref_factor);

    %Run the clustering algo. Note that apcluster requires a similarity matrix,
    %so change the sign on the distance matrix
    [idx, ~, ~, ~] = apcluster(-1*ds.mat,preference_vector);
else
    %%Bisection
    %Set preferences vector
    idx = 0;
    low_quantile = 0;
    high_quantile = .1;
    K = args.num_exemplar;
        
    if args.verbose
        fprintf('Running Bisection to target %d exemplars \n',K);
    end
    
    %Run a bisection on the preference quantile to search for the exact
    %number of exemplars
    while numel(unique(idx)) ~= K
        
        temp_quantile = (high_quantile + low_quantile)/2;
        
        if args.verbose
            fprintf('Trying a quantile of %.5f \n',temp_quantile);
        end
        
        preference_vector = mk_preference_vector(-1*ds.mat,...
            ds.rid,...
            known,...
            preference,...
            temp_quantile,...
            args.pref_factor);
        
        [idx, ~, ~, ~] = apcluster(-1*ds.mat,preference_vector);
        
        if args.verbose
            fprintf('Got %d exemplars\n',numel(unique(idx)));
        end
        
        if numel(unique(idx)) > K
            low_quantile = temp_quantile;
        elseif numel(unique(idx)) < K
            high_quantile = temp_quantile;
        end
    end
        
end

exemplar_idx = unique(idx);
exemplars = ds.rid(exemplar_idx);
membership = ds.rid(idx);

end

function pref_vector = mk_preference_vector(sim,objects,known,preferred,t,pref_factor)
%Create a vector specifying the initial preferance each object has for
%being an exemplar.
%
%The higher the preference, the more likely it is to be chosen as an
%exemplar. We currently support three levels of hierarchy: known, preferred
%and everything else. Known objects get an (essentially) infinite
%preference. Preferred objects get a preference that is pref_factor better
%than everything else
n = numel(objects);

%Get preference for 'everything else'
[pref_min, pref_max] = preferenceRange(sim);
const = t*pref_min + (1-t)*pref_max;

pref_vector = pref_factor*const*ones(n,1);

candidate_idx = ismember(objects,preferred);
pref_vector(candidate_idx) = const;

known_idx = ismember(objects,known);
pref_vector(known_idx) = -1*const*100;


end