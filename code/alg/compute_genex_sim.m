function [sim, excl] = compute_genex_sim(varargin)
% compute_genex_sim - given a set of cell lines, compute the baseline gene expression
% using a particular metric on CCLE data
% 
% [sim, excl] = compute_genex_sim(varargin)
% Returns:
%     sim       - a gctx with cid = rid = set of cell lines with element ij in mat given by the 
%                 distance between cell line i and j using the given metric in the given CCLE ds
%     excl      - any cell lines input that were excluded
%
% Parameters:
%     cells     - A cell array of cell lines.  If empty (default), computes on all available cell lines
%     metric    - string; the distance function.  Supported: 'euclidean' (default), 'cosine' (distance: 1-cosine)
%     mahal_num_pcs - Dimension of PCA to perform when using mahalanobis
%           distance. Default 30
%     ccleds    - string, which determines which ccleds to use.  Supported: rnaseq (default)
%     gene_space - .grp file or cell array. Space of genes in which to
%                  compute distances. Default is all genes in the ds.
%     plot      - Boolean, make plots summarizing distances; default 0
%     outdir    - directory to save plots and data; default ../analysis/baseline_gex_distance
%     savefiles - boolean, whether to write files or just return them; default 0
%     norm      - normalize the stdev of all genes to 1. 

pnames = {'cells', ...
    'metric', ...
    'mahal_num_pcs',...
    'ccleds', ...
    'gene_space',...
    'plot', ...
    'outdir', ...
    'savefiles', ...
    'norm'};
dflts = {{}, ...
    'euclidean', ...
    30,...
    'rnaseq2015', ...
    '',...
    0, ...
    '/cmap/projects/cell_line_diversity/analysis/baseline_gex_distance', ...
    0, ...
    1};
args = parse_args(pnames, dflts, varargin{:});


[ds, annot, excl] = get_gex_data(args);
if args.norm
  ds = norm_data(ds);
end

sim = calc_gex_sim(ds, args);

if args.plot
  mk_gex_plots(sim, annot, ds, args);
end

end


function [ds, annot, excl] = get_gex_data(args)
  datapath = '/cmap/data/vdb/cline';
  trainingpath = '/cmap/projects/cell_line_diversity/data/';

  switch args.ccleds
    case 'rnaseq'
      ds = parse_gctx(fullfile(datapath, 'ccle_rnaseq_rpkm_n932x23686.gctx'));
    case 'rnaseq2015'
      ds = parse_gctx(fullfile(trainingpath, 'CCLE_BASELINE_RNASEQ_L1KFULL_RPKM_LOG2_LABELED_n1022x12450.gctx'));
    case 'affy'
      ds = parse_gctx(fullfile(datapath, 'cline_gene_n1515x12716.gctx'));
    case 'affyzs'
      ds = parse_gctx(fullfile(datapath, 'cline_gene_zs_n1515x12716.gctx'));
  end

  %Gene filter
  gene_space = get_array_input(args.gene_space,ds.rid);
  ds = ds_slice(ds,'rid',gene_space,...
      'ignore_missing',true);
  
  cix = ifelse(isempty(args.cells), 1:numel(ds.cid), find(ismember(ds.cid, args.cells)));
  cset = ds.cid(cix);
  ds = gctsubset(ds, 'csubset', cix);
  excl = setdiff(args.cells, cset);
  annot = parse_tbl(fullfile(datapath, 'clinedb.txt'));
end

function ds = norm_data(ds)
    s = std(ds.mat, 0, 2);
    s(s == 0) = 1;
    ds.mat = ds.mat./repmat(s, 1, numel(ds.cid));
end

function ret = calc_gex_sim(ds, args)
  switch lower(args.metric)
    case 'euclidean'
        sim = squareform(pdist(ds.mat', 'euclidean')); 
    
    case 'mahal'
        data = ds.mat';
        [~, pca_data] = pca(data);
        pca_data = pca_data(:,1:args.mahal_num_pcs);
        sim = squareform(pdist(pca_data, 'mahalanobis'));
        
    case 'seuclidean'
        sim = squareform(pdist(ds.mat', 'seuclidean')); 

    case 'cosine'
        sim = squareform(pdist(ds.mat', 'cosine'));
    
    case 'spearman'
        %Note: dont transpose matrix
        sim = 1 - fastcorr(ds.mat, 'type', 'spearman');
        
    case 'pearson'
        %Note: dont transpose matrix
        sim = 1 - fastcorr(ds.mat, 'type', 'pearson'); 
  end

  ret = mkgctstruct(sim, 'rid', ds.cid, 'cid', ds.cid, 'rdesc', ds.cdesc, 'cdesc', ds.cdesc, 'chd', ds.chd, 'rhd', ds.chd);
end


function mk_gex_plots(sim, annot, ds, args)
  mx = tsne_d(sim.mat, [], 2, 30);

  ix = find(ismember(annot.cell_id, sim.cid));
  t = struct_subset(annot, ix);
  jx = match_vectors(sim.cid, t.cell_id);
  t = struct_subset(t, jx);

  [u,c,g] = cellcount(t.cell_lineage);
  [~,cx] = sort(c, 'descend');
  mycmap = varycolor(11);

  figure('Position', [10 10 1200 1000]); hold on;
  for k = 1:10
    scatter(mx(g{cx(k)},1), mx(g{cx(k)},2), 150, mycmap(k,:), '.');
  end
  for k = 1:10
    scatter(mx(g{cx(k+10)},1), mx(g{cx(k+10)}, 2), 10, mycmap(k,:), 's');
  end
  go = vertcat(g{cx(21:end)});
  scatter(mx(go,1), mx(go,2), 10, mycmap(11,:), '>');
  
  legend(vertcat(u(cx(1:20)), 'Other'), 'Location', 'NorthEast');
  xl = get(gca, 'xlim'); grid on;
  set(gca, 'xlim', [xl(1) xl(2)+(0.2*(xl(2) - xl(1)))]);
  title(sprintf('tSNE Plot for %s,%s metric = %s', args.ccleds, ifelse(args.norm, ' normalized',' '), args.metric));
  print(gcf, '-dpng', '-r250', fullfile(args.outdir, sprintf('tsne_%s_%s%s_all.png', args.ccleds, args.metric, ifelse(args.norm, '_norm', ''))));

end
