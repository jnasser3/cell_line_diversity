function [sim, ds, pcspace, excl] = compute_achilles_sim(varargin)
% compute_achilles_sim - given a set of cell lines, compute the distance in Achilles
% 
% [sim, excl] = compute_achilles_sim(varargin)
% Returns:
%     sim       - a gctx with cid = rid = set of cell lines with element ij in mat given by the 
%                 distance between cell line i and j using the given metric in the given CCLE ds
%                 Note that it is a DISTANCE, not a similarity, i.e. smaller values indicate greater
%                 similarity.  Note that correlations like Spearman return (1 - spearman).
%     ds        - gctx struct with the appropriate data from Achilles
%     pcspace   - a linear transformation of the achilles dataset using principal components
%     excl      - any cell lines input that were excluded
%     npcs      - number of principal components to use in the dimensional reduction
%
% Parameters:
%     cells     - A cell array of cell lines.  If empty (default), computes on all available cell lines
%     metric    - string; the distance function.  Supported: 'euclidean' (default), 'cosine' (distance: 1-cosine)
%     achds     - string, which determines which ccleds to use.  Supported: QC (default), 'FC' (fold change; 
%                 raw, unnormalized data)
%     plot      - Boolean, make plots summarizing distances; default 0
%     outdir    - directory to save plots and data; default ../analysis/baseline_gex_distance
%     savefiles - boolean, whether to write files or just return them; default 0
%     norm      - normalize the stdev of all genes to 1. 

pnames = {'cells', ...
    'metric', ...
    'achds', ...
    'plot', ...
    'outdir', ...
    'savefiles', ...
    'norm', ...
    'metricparam'};
dflts = {{}, ...
    'cosine', ...
    'QC', ...
    0, ...
    '/cmap/projects/cell_line_diversity/analysis/baseline_gex_distance', ...
    0, ...
    0, ...
    30};
args = parse_args(pnames, dflts, varargin{:});
args.achdir = '/cmap/projects/rnai_analysis/genedeath/achilles';

[ds, excl, pcspace, pcmap] = get_ach_data(args);
if args.norm
  ds = norm_data(ds);
end

sim = calc_ach_sim(ds, args);

if args.plot
  mk_ach_plots(sim, annot, ds, args);
end

end


function [ds, excl, pcspace, coeff] = get_ach_data(args)
if strcmpi(args.achds, 'QC')
  ach1 = parse_gctx(fullfile(args.achdir, 'Achilles_QC_v2.9_cmap_n102x101669.gctx'));
  ach2 = parse_gctx(fullfile(args.achdir, 'Achilles_QC_v2.4.rnai2.gct'));
  ach2.rhd{1,1} = 'gene'; ach2.rhd{2,1} = 'hp_seq';
  ach2.rdesc(:,2) = strtok(ach2.rid, '_');
  ach2.cdesc = horzcat(ach2.cdesc, ach2.cid);
  ach2.chd = vertcat(ach2.chd, 'full_cell_id');
  ach2.cid = strtok(ach2.cid, '_');

  [~,ia,ib] = intersect(ach1.rid, ach2.rid);
  ach1 = gctsubset(ach1, 'rsubset', ia);
  ach2 = gctsubset(ach2, 'rsubset', ib);
  ds = gctmerge_tool(ach1, ach2);

  ds.rdict = containers.Map(ds.rhd, 1:numel(ds.rhd));
  ds.cdict = containers.Map(ds.chd, 1:numel(ds.chd));

  excl = [];
  if ~isempty(args.cells)
    excl = setdiff(ds.cid, args.cells);
    ds = gctsubset(ds, 'csubset', find(ismember(ds.cid, args.cells)));
  end
    
elseif strcmpi(args.achds, 'FC')
  % Not yet supported
  ach = parse_gctx('Achilles_v2.19.1_FC.rnai.gct');
  ach = parse_gctx(fullfile(args.topdir, 'achilles', 'Achilles_v2.4_FC.rnai.gct'));
  error('Fold change data not yet supported because there exist replicates');
end

nanix = sum(isnan(ds.mat),2);
ds = gctsubset(ds, 'rsubset', find(nanix == 0));

% do PCA
[coeff, pcspace, lt, tsq, expl] = pca(ds.mat');
figure; plot(cumsum(expl));  grid on;
xlabel('Sum of first N Principal Components');
ylabel('Total variance explained by first N components');
title({sprintf('Achilles %s data - cumulative variance explained by principal components', args.achds); ...
    sprintf('N cell lines: %d, N shRNAs: %d', numel(ds.cid), numel(ds.rid))});

% Remove PC1 - it essentially just separates the two datasets
pcspace(:,1) = 0;
ds.mat = (pcspace * coeff')' + repmat(mean(ds.mat,2), 1, size(ds.mat,2));

end

function cows()
  datapath = '/cmap/data/vdb/cline';
  trainingpath = '/cmap/projects/cell_line_diversity/data/';

  switch args.ccleds
    case 'rnaseq'
      ds = parse_gctx(fullfile(datapath, 'ccle_rnaseq_rpkm_n932x23686.gctx'));
    case 'rnaseq2015'
      ds = parse_gctx(fullfile(trainingpath, 'CCLE_BASELINE_RNASEQ_L1KFULL_RPKM_LOG2_n1022x12450.gctx'));
    case 'affy'
      ds = parse_gctx(fullfile(datapath, 'cline_gene_n1515x12716.gctx'));
    case 'affyzs'
      ds = parse_gctx(fullfile(datapath, 'cline_gene_zs_n1515x12716.gctx'));
  end

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

function ret = calc_ach_sim(ds, args)
  switch lower(args.metric)
    case 'euclidean'
      sim = squareform(pdist(ds.mat', 'euclidean'));
    
    case 'cosine'
      amag = diag(diag(sqrt(ds.mat' * ds.mat)));
      sim = 1 - amag \ ds.mat' * ds.mat / amag;
      % sim = squareform(pdist(ds.mat', 'cosine'));

    case 'seuclidean'
      sim = squareform(pdist(ds.mat', 'seuclidean'));

    case 'mahal'
      [~, pca_data] = pca(ds.mat');
      pca_data = pca_data(:, 1:args.metric_param);
      sim = squareform(pdist(pca_data, 'mahalanobis'));

    case 'spearman'
      sim = 1 - fastcorr(ds.mat, 'type', 'Spearman');

    case 'pearson'
      sim = 1 - fastcorr(ds.mat, 'type', 'Pearson');
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
