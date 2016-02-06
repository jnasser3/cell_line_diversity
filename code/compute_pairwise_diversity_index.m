function pdi = compute_pairwise_diversity_index(ds, varargin)
%Given a ds of signatures, computes the pairwise diversity index for all
%cell lines
%
%Input: ds: a genes x perturbagens struct containing the gex data
%       cell_lines: cell array or .grp of lines for which to compute
%                   pairwise diversities
%       make_plot: Boolean, whether to make plots. Default 0
%       save_plot: Boolean, whether to save plots to disk. Default 0
%       plot_dir: Directory to save plots to. 
%                 Default
%                 /cmap/projects/cell_line_diversity/analysis/pdi/fig
%
%Output: pdi: a structure of Pairwise Diversity Index. pdi.mat is a
%cell_lines x cell_lines matrix with pdi.mat(i,j) being the diversity index
%between cell_line i and cell_line j

params = {'cell_lines',...
          'make_plot',...
          'save_plot',...
          'plot_dir'};
dflts = {'/cmap/projects/cell_line_diversity/data/phase1_core_lines.grp',...
         false,...
         false,...
         '/cmap/projects/cell_line_diversity/analysis/pairwise_diversity_index/fig'};
args = parse_args(params,dflts,varargin{:});

if iscell(args.cell_lines)
    lines = args.cell_lines;
else
    lines = parse_grp(args.cell_lines);
end

%initialize pdi
mat = zeros(numel(lines));
pdi = mkgctstruct(mat,'rid',lines,'cid',lines);

%Compute the diversity index for each pair of lines. 
for ii = 1:numel(lines)
    idx1 = strcmp(ds.cdesc(:,ds.cdict('cell_id')),lines(ii));
    ds1 = ds_slice(ds,'cidx',find(idx1));
    
    for jj = (ii+1):numel(lines)
        idx2 = strcmp(ds.cdesc(:,ds.cdict('cell_id')),lines(jj));
        ds2 = ds_slice(ds,'cidx',find(idx2));
        
        pdi.mat(ii,jj) = compute_diversity_index_two(ds1,ds2,...
            'make_plot',args.make_plot,...
            'save_plot',args.save_plot,...
            'plot_dir',args.plot_dir);
    end
end

%Fill in the lower triangle
pdi.mat = pdi.mat + pdi.mat';
end

