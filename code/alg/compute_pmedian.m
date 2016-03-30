function [exemplars,loss,si] = compute_pmedian(ds,p,varargin)
%exemplars = compute_pmedian(pdi,p,varargin)
%
%Returns approximate solution to the p-medians problem. Given pairwise
%distances between n objects, returns a set of p ( < n) objects that
%minimize the sum of the minimum distances of each object to it's nearest
%landmark point.
%
%A solution to the conditional p-medians problem is also supported. In this
%case we have a set of known landmark points and wish to choose an
%additional p landmarks.
%
%Inputs:
%       ds: an n x n struct of pairwise distances
%       p: the number of landmark objects to return
%       exclude: a cell array or .grp of objects to exclude. Default ''
%       known_ex: a cell array or .grp of known exemplars. Default ''
%       exemplar_method: How to incorporate exemplars into conditional
%                   pmedians. Options: {'hack','xform_pdist'}
%       algorithm: Specifies which algorithm to use. Options:
%                  'mfi': vertex substitution type algorithm
%                  'sa': simulated annealing
%       mfi_restarts: Number of times to run the mfi algorithm. Default 1.
%       verbose: Print logging to screen. Default true
%
%Outputs:
%       exemplars: cell array of the estimated exemplars. This *DOES NOT*
%               include the known_ex.

params = {'exclude',...
          'known_ex',...
          'exemplar_method',...
          'algorithm',...
          'mfi_restarts',...
          'verbose'};
dflts = {'',...
    '',...
    'xform_pdist',...
    'mfi',...
    1,...
    true};
args = parse_args(params,dflts,varargin{:});

addpath('/cmap/projects/cell_line_diversity/code/pcluster');

%Subset ds if needed
exclude = get_array_input(args.exclude,'');
if ~isempty(exclude)
    objects = setdiff(ds.rid,exclude);
    ds = ds_slice(ds,'cid',objects,'rid',objects);
end

%Transform into conditional pmedians problem
wt = zeros(1,numel(ds.rid));
if numel(args.known_ex) > 0
    switch args.exemplar_method
        case 'xform_pdist'
            ds = mk_conditional_pdist(ds,args.known_ex);
        case 'hack'
            wt = mk_known_wts(ds.rid,args.known_ex);
    end
end

%Display information to user
if args.verbose
    fprintf('Computing exemplars using the vertex substitution method');
end

%Run pmedians
switch args.algorithm
    case 'mfi'
        [ex_idx, loss, ~, si] = MFI_FC(ds.mat,p,wt,args.mfi_restarts);
    case 'sa'
        %wt = mk_known_wts(objects,args.known_ex);
        [ex_idx, loss, ~, si] = SA_FC(ds.mat,p,wt,.8);
end

exemplars = ds.rid(ex_idx);

end

function wt = mk_known_wts(objects,lm_objects)
%Given a subset of the objects which are known landmarks, returns a weight
%vector promoting the inclusion of these landmarks.

n = numel(objects);
wt = zeros(n,1);
lm_idx = ismember(objects,lm_objects);
wt(lm_idx) = -100000;
end