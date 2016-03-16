function [pdi, all_corrs] = compute_pairwise_diversity_index(ds, varargin)
%pdi = compute_pairwise_diversity_index(ds, varargin)
%
%Given a ds of signatures, computes the pairwise diversity index for all
%cell lines
%
%Input: ds: a genes x perturbagens struct containing the gex data
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

params = {'cell_lines',...
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
dflts = {'',...
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

%get cell lines and pert ids
lines = get_array_input(args.cell_lines,unique(ds.cdesc(:,ds.cdict('cell_id'))));
pert_ids = get_array_input(args.pert_ids,unique(ds.cdesc(:,ds.cdict('pert_id'))));

%Make sure we have at least two cell lines
assert(numel(lines) >= 2, 'Must have at least two cell lines')
assert(numel(pert_ids) >= 2, 'Must have at least two pert ids')

%initialize pdi
mat = zeros(numel(lines));
pdi = mkgctstruct(mat,'rid',lines,'cid',lines);

%Compute the diversity index for each pair of lines. 
for ii = 1:numel(lines)
    idx1 = strcmp(ds.cdesc(:,ds.cdict('cell_id')),lines(ii));
    ds1 = ds_slice(ds,'cidx',find(idx1));
    ds1 = subset_ds(ds1,...
        'pert_id',pert_ids);

    for jj = (ii+1):numel(lines)
        idx2 = strcmp(ds.cdesc(:,ds.cdict('cell_id')),lines(jj));
        ds2 = ds_slice(ds,'cidx',find(idx2));
        ds2 = subset_ds(ds2,...
            'pert_id',pert_ids);

        [DI, all_corrs] = compute_diversity_index_two(ds1,ds2,args);
        
        pdi.mat(ii,jj) = DI;
    end
end

%Fill in the lower triangle
pdi.mat = pdi.mat + pdi.mat';

%save the gctx and make config file
if args.save_gctx
    mkgctx(fullfile(args.gctx_dir,'pdi.gctx'),pdi)
    mkconfigyaml(fullfile(args.gctx_dir,'config.yaml'),args)
end


end

function ds = subset_ds(ds,field,to_include)
%remove samples that have -666 for distil_cc_q75
bad_idx = cell2mat((ds.cdesc(:,ds.cdict('distil_cc_q75')))) == -666;
ds = ds_slice(ds, 'cid', ds.cid(bad_idx), 'exclude_cid', true);

%subset ds to field to include
good_idx = find(ismember(ds.cdesc(:,ds.cdict(field)),to_include));
ds = ds_slice(ds, 'cidx', good_idx);
end