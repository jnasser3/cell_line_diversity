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
%     ccleds    - string, which determines which ccleds to use.  Supported: rnaseq (default)
%     plot      - Boolean, make plots summarizing distances; default 0
%     outdir    - directory to save plots and data; default ../analysis/baseline_gex_distance
%     savefiles - boolean, whether to write files or just return them; default 0

pnames = {'cells', ...
    'metric', ...
    'ccleds', ...
    'plot', ...
    'outdir', ...
    'savefiles'};
dflts = {{}, ...
    'euclidean', ...
    'rnaseq', ...
    0, ...
    '/cmap/projects/cell_line_diversity/analysis/baseline_gex_distance', ...
    0};
args = parse_args(pnames, dflts, varargin{:});


[ds, excl] = get_gex_data(args);
sim = calc_gex_sim(ds, args);

if args.plot
  mk_gex_plots(sim, ds, args);
end

end


function [ds, excl] = get_gex_data(args)
  datapath = '/cmap/data/vdb/cline';

  switch args.ccleds
    case 'rnaseq'
      ds = parse_gctx(fullfile(datapath, 'ccle_rnaseq_rpkm_n932x23686.gctx'));
    case 'affy'
      ds = parse_gctx(fullfile(datapath, 'cline_gene_n1515x12716.gctx'));
    case 'affyzs'
      ds = parse_gctx(fullfile(datapath, 'cline_gene_zs_n1515x12716.gctx'));
  end

  cix = ifelse(isempty(args.cells), ds.cid, ds.cid(ismember(ds.cid, args.cells)));
  ds = gctsubset(ds, 'csubset', cix);
  excl = setdiff(args.cells, cix);
end


function sim = calc_gex_sim(ds, args)
  switch args.metric
    case 'euclidean'
      % Ghetto
      sim = zeros(numel(ds.cid), numel(ds.cid));
      for k = 1:numel(ds.cid)
        sim(:,k) = sqrt(sum((ds.mat - repmat(ds.mat(:,k), 1, numel(ds.cid))).^2));
      end

    case 'cosine'
      % 1-cosine(x)
      amag = diag(diag(sqrt(ds.mat' * ds.mat)));
      dist = 1 - amag \ ds.mat' * ds.mat / amag;

  end
end


function mk_gex_plots(sim, ds, args)

end
