function pdi = compute_pairwise_diversity_index(ds, varargin)
%Given a ds of signatures, computes the pairwise diversity index for all
%cell lines
%
%Input: cell_lines: cell array or .grp of lines for which to compute
%                   pairwise diversities
%
%Output: pdi: a structure of Pairwise Diversity Index

params = {'cell_lines',...
          'make_plot'};
dflts = {'/cmap/projects/cell_line_diversity/data/phase1_core_lines.grp',...
         'false'};
args = parse_args(params,dflts,varargin{:});

if iscell(args.cell_lines)
    lines = args.cell_lines;
else
    lines = parse_grp(args.cell_lines);
end

%initialize pdi
mat = zeros(numel(lines));
pdi = mkgctstruct(mat,'rid',lines,'cid',lines);

for ii = 1:numel(lines)
    idx1 = strcmp(ds.cdesc(:,ds.cdict('cell_id')),lines(ii));
    ds1 = ds_slice(ds,'cidx',find(idx1));
    
    for jj = (ii+1):numel(lines)
        idx2 = strcmp(ds.cdesc(:,ds.cdict('cell_id')),lines(jj));
        ds2 = ds_slice(ds,'cidx',find(idx2));
        
        pdi.mat(ii,jj) = compute_diversity_index_two(ds1,ds2,...
            'make_plot',args.make_plot);
    end
end

%Fill in the lower triangle
pdi.mat = pdi.mat + pdi.mat';
end

