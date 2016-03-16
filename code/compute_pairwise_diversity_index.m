function [pdi, all_corrs] = compute_pairwise_diversity_index(ds, varargin)
%pdi = compute_pairwise_diversity_index(ds, varargin)
%
%Given a ds of signatures, computes the pairwise diversity index for all
%cell lines
%
%Input: ds: a genes x perturbagens struct containing the gex data
%       prune_ds: Boolean. If true we only consider signatures that have a
%                   prune_field appearing across ALL cell lines in the ds. Setting
%                   prune_ds to true ensures that all pairwise distances between lines
%                   are computed using the same set of perturbagens. Default true.
%       prune_field: The column field to match on. Default is 'pert_id_dose_time'
%       metric: The metric used to quantify diversity.
%       corr_type: Type of correlation to use. Default is spearman.
%       kde_method: Function used to compute kernel density estimate.
%                   Options are {'matlab','fast_kde'}. Default is 'fast_kde'
%       cell_lines: cell array or .grp of lines for which to compute
%                   pairwise diversities. Default: all in ds.
%       pert_ids: cell array or .grp of pert_ids on which to compute
%                   pairwise diversities.
%       make_plot: Boolean, whether to make plots. Default false
%       show_plot: Boolean, whether to make plot visible. Default false
%       save_plot: Boolean, whether to save plots to disk. Default false
%       plot_dir:  Directory to save plots to. 
%       save_gctx: Boolean, whether to save the pairwise distance gctx.
%                  Defualt false
%       gctx_dir: Directory to save gctx to
%                
%Output: pdi: a structure of Pairwise Diversity Index. pdi.mat is a
%        cell_lines x cell_lines matrix with pdi.mat(i,j) being the diversity index
%        between cell_line i and cell_line j

params = {'prune_ds',...
          'prune_field',...
          'cell_lines',...
          'pert_ids',...
          'metric',...
          'corr_type',...
          'kde_method',...
          'make_plot',...
          'show_plot',...
          'save_plot',...
          'plot_dir',...
          'save_gctx',...
          'gctx_dir'};
dflts = {true,...
         'pert_id_dose_time',...
         '',...
         '',...
         'rel_bioa_corr',...
         'spearman',...
         'fast_kde',...
         false,...
         false,...
         false,...
         '.',...
         false,...
         '.'};
args = parse_args(params,dflts,varargin{:});

%% Preprocessing
%get cell lines and pert ids
cell_lines = get_array_input(args.cell_lines,unique(ds.cdesc(:,ds.cdict('cell_id'))));
pert_ids = get_array_input(args.pert_ids,unique(ds.cdesc(:,ds.cdict('pert_id'))));

%Make sure we have at least two cell lines
assert(numel(cell_lines) >= 2, 'Must have at least two cell lines')
assert(numel(pert_ids) >= 2, 'Must have at least two pert ids')

%subset ds to pert_ids
if ~isempty(pert_ids)
    ds = subset_ds(ds,'pert_id',pert_ids);
end

%Remove signatures that have -666 for ccq75 (ones that have 1 replicate)
bad_idx = cell2mat((ds.cdesc(:,ds.cdict('distil_cc_q75')))) == -666;
ds = ds_slice(ds, 'cid', ds.cid(bad_idx), 'exclude_cid', true);

%If requested subset ds to only signatures that appear once in each line.
if args.prune_ds
    ds = select_training_data_function(ds,'cell_lines',cell_lines,'match_field',args.prune_field);  
end

%initialize pdi
mat = zeros(numel(cell_lines));
pdi = mkgctstruct(mat,'rid',cell_lines,'cid',cell_lines);

%% Compute the diversity index for each pair of lines. 
for ii = 1:numel(cell_lines)
    ds1 = subset_ds(ds,'cell_id',cell_lines(ii));

    for jj = (ii+1):numel(cell_lines)
        ds2 = subset_ds(ds,'cell_id',cell_lines(jj));
        pdi.mat(ii,jj) = compute_diversity_index_two(ds1,ds2,args);
    end
end

%Fill in the lower triangle
pdi.mat = pdi.mat + pdi.mat';

%% Save the gctx and make config file
if args.save_gctx
    mkgctx(fullfile(args.gctx_dir,'pdi.gctx'),pdi)
    mkconfigyaml(fullfile(args.gctx_dir,'config.yaml'),args)
end


end

function ds = subset_ds(ds,field,to_include)
    good_idx = find(ismember(ds.cdesc(:,ds.cdict(field)),to_include));
    ds = ds_slice(ds, 'cidx', good_idx);
end