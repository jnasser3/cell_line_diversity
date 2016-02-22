function pdi = compute_pairwise_diversity_index(ds, varargin)
%pdi = compute_pairwise_diversity_index(ds, varargin)
%
%Given a ds of signatures, computes the pairwise diversity index for all
%cell lines
%
%Input: ds: a genes x perturbagens struct containing the gex data
%       metric: The metric used to quantify diversity.
%       cell_lines: cell array or .grp of lines for which to compute
%                   pairwise diversities. Default: all in ds.
%       make_plot: Boolean, whether to make plots. Default 0
%       show_plot: Boolean, whether to make plot visible. Default 0
%       save_plot: Boolean, whether to save plots to disk. Default 0
%       plot_dir:  Directory to save plots to. 
%                
%
%Output: pdi: a structure of Pairwise Diversity Index. pdi.mat is a
%cell_lines x cell_lines matrix with pdi.mat(i,j) being the diversity index
%between cell_line i and cell_line j

params = {'cell_lines',...
          'metric',...
          'make_plot',...
          'show_plot',...
          'save_plot',...
          'plot_dir'};
dflts = {'',...
         'basic_wtd_corr',...
         false,...
         false,...
         false,...
         ''};
args = parse_args(params,dflts,varargin{:});

if isempty(args.cell_lines)
    lines = unique(ds.cdesc(:,ds.cdict('cell_id')));
elseif iscell(args.cell_lines)
    lines = args.cell_lines;
else
    lines = parse_grp(args.cell_lines);
end

%Make sure we have at least two cell lines
assert(numel(lines) >= 2, 'Must have at least two cell lines')

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
            'metric',args.metric,...
            'make_plot',args.make_plot,...
            'show_plot',args.show_plot,...
            'save_plot',args.save_plot,...
            'plot_dir',args.plot_dir);
    end
end

%Fill in the lower triangle
pdi.mat = pdi.mat + pdi.mat';
end

