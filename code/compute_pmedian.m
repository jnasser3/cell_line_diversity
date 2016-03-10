function [lm_objects,loss,si] = compute_pmedian(ds,p,varargin)
%locs = compute_pmedian(pdi,p,varargin)
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
%       objects: a cell array or .grp of objects to consider. Default to
%                all objects in the ds
%       known_lm: a cell array or .grp of known lm objects. Default to none
%       algorithm: Specifies which algorithm to use. Default mfi
%       mfi_restarts: Number of times to run the mfi algorithm. Default 1.
%
%Outputs:
%       lm_objects: cell array of landmark objects

params = {'objects',...
          'known_lm',...
          'algorithm',...
          'mfi_restarts'};
dflts = {'',...
    '',...
    'mfi',...
    1};
args = parse_args(params,dflts,varargin{:});

addpath('/cmap/projects/cell_line_diversity/code/pcluster');

%get objects to use
objects = get_array_input(args.objects,unique(ds.rid));

%Run pmedians
switch args.algorithm
    case 'mfi'
        wt = mk_known_wts(objects,args.known_lm);
        [lm_idx, loss, ~, si] = MFI_FC(ds.mat,p,wt,args.mfi_restarts);
end

lm_objects = ds.rid(lm_idx);

end

function wt = mk_known_wts(objects,lm_objects)
%Given a subset of the objects which are known landmarks, returns a weight
%vector promoting the inclusion of these landmarks.

n = numel(objects);
wt = zeros(n,1);
lm_idx = ismember(objects,lm_objects);
wt(lm_idx) = -100000;
end