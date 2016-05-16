function [pdi, data_struct] = compute_pairwise_diversity_index(ds, varargin)
%pdi = compute_pairwise_diversity_index(ds, varargin)
%
%Given a ds of signatures, computes the pairwise diversity index for all
%cell lines. ds should be a genex x samples GCT struct containing signatures across > 1
%cell lines. Each cell line must be annotated in the cell_id field of
%cdesc.
%
%Input: ds: a genes x samples struct containing the gex data
%       prune_ds: Boolean. If true we only consider signatures that have a
%                   prune_field appearing across ALL cell lines in the ds. Setting
%                   prune_ds to true ensures that all pairwise diversities between lines
%                   are computed using the same set of perturbagens. Default true.
%       prune_field: The column field to prune on. Default is 'pert_id'
%       metric: The metric used to quantify diversity. Default
%           rel_bio_corr. Other options are not recommended.
%       corr_type: Type of correlation to use. Default is spearman.
%       corr_bkg_type: Background of correlations against which to compute 
%                   diversities. Options are 'global' or 'local'. Default
%                   is 'global'.
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
%       save_gctx: Boolean, whether to save the pairwise diversity gctx.
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
          'corr_bkg_type',...
          'kde_method',...
          'make_plot',...
          'show_plot',...
          'save_plot',...
          'plot_dir',...
          'save_gctx',...
          'gctx_dir'};
dflts = {true,...
         'pert_id',...
         '',...
         '',...
         'rel_bioa_corr',...
         'spearman',...
         'global',...
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

%Subset ds to pert_ids
ds = ds_subset(ds,'pert_id',pert_ids);

%Remove signatures that have less than two replicates
bad_idx = cell2mat((ds.cdesc(:,ds.cdict('distil_nsample')))) < 2;
ds = ds_slice(ds, 'cid', ds.cid(bad_idx), 'exclude_cid', true);

%If requested subset ds to only signatures that appear once in each line.
if args.prune_ds
    ds = select_training_data_function(ds,'cell_lines',cell_lines,'match_field',args.prune_field);
end

%Compute global correlation background
if strcmp(args.corr_bkg_type,'global')
    %ds_corr_bkg = select_training_data_function(ds,'cell_lines',cell_lines,'match_field',args.prune_field);
    ds_corr_bkg = ds;
    
%     subsample = 100;
%     idx = randperm(numel(ds.cid),subsample); 
    
    corrs = fastcorr(ds_corr_bkg.mat(:,:),'type',args.corr_type);

    [~, corr_bkg] = mat2vec_diag(corrs);
    
    [corr_xform,corr_mesh] = compute_xform(corr_bkg,'cdf','fast_kde');
    args.corr_xform = corr_xform;
    args.corr_mesh = corr_mesh;
end

%initialize pdi
mat = zeros(numel(cell_lines));
pdi = mkgctstruct(mat,'rid',cell_lines,'cid',cell_lines);

%% Compute the diversity index for each pair of lines. 
for ii = 1:numel(cell_lines)
    ds1 = ds_subset(ds,'column','cell_id',cell_lines(ii));

    for jj = (ii+1):numel(cell_lines)
        ds2 = ds_subset(ds,'column','cell_id',cell_lines(jj));
        
        [this_pdi,data_struct] = compute_diversity_index_two(ds1,ds2,args);
        
        assert(~isnan(this_pdi),'Pairwise diversity evaluated to NaN')
        pdi.mat(ii,jj) = this_pdi;
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

% function ds = subset_ds(ds,field,to_include)
%     good_idx = find(ismember(ds.cdesc(:,ds.cdict(field)),to_include));
%     ds = ds_slice(ds, 'cidx', good_idx);
% end