function evaluate_cpc006(varargin)

params = {'cell_id', ...
    'outdir', ...
    'flatpc1'};
dflts = {'all', ...
    '/cmap/projects/cell_line_diversity/analysis/cpc006/bycell', ...
    0};
args = parse_args(params, dflts, varargin{:});

args.lms = parse_grp('/cmap/data/vdb/spaces/lm_epsilon_n978.grp');


[ds, rds] = get_cpc6_data(args);

[brewpfx, cfq10, dmsoq10] = analyze_cpc6_cell(args, ds, rds);

end


function [u, cellfrac, dmsofrac] = analyze_cpc6_cell(args, ds, rds)
  [u,c,g] = cellcount(ds.cdesc(:, ds.cdict('brew_prefix')));
  cellfrac = [];
  dmsofrac = [];
  figure;
   
  for ii = 1:numel(u)
    rx = find(ismember(rds.cid, vertcat(ds.cdesc{g{ii}, ds.cdict('distil_id')})));
    rsub = gctsubset(rds, 'csubset', rx);

    rcorrs = fastcorr(rsub.mat, 'type', 'Spearman');
    rtcorrs = tri2vec(rcorrs);
    
    [uu,cc,gg] = cellcount(rsub.cdesc(:, rsub.cdict('pert_id')));
    dmsoix = cellstrfind(uu, 'DMSO');
    pertcorrs = [];
    for k = 1:numel(uu)
      pertcorrs{k} = tri2vec(rcorrs(gg{k}, gg{k}));
    end

    pertranks = cellfun(@(x) arrayfun(@(y) mean(rtcorrs > y), x), pertcorrs, 'UniformOutput', 0);
    pertq = fdr_calc(vertcat(pertranks{setdiff(1:numel(uu), dmsoix)}));
    cellfrac = vertcat(cellfrac, mean(pertq < 0.1));

    dmsocorrs = pertcorrs{dmsoix};
    dmsoranks = cellfun(@(x) arrayfun(@(y) mean(dmsocorrs >= y), x), pertcorrs, 'UniformOutput', 0);
    dmsoq = fdr_calc(vertcat(dmsoranks{setdiff(1:numel(uu), dmsoix)}));
    dmsofrac = vertcat(dmsofrac, mean(dmsoq < 0.1));

    % Make figure
    clf; hold on; grid on;
    [b,a] = ksdensity(rtcorrs, 'bandwidth', 0.01);
    plot(a,b,'k--','LineWidth',2.5);
    [b,a] = ksdensity(pertcorrs{dmsoix}, 'bandwidth', 0.01);
    plot(a,b,'b--','LineWidth',2.5);
    [b,a] = ksdensity(vertcat(pertcorrs{setdiff(1:numel(uu), dmsoix)}), 'bandwidth', 0.01);
    plot(a,b,'r','LineWidth',2.5);
    xlabel('Spearman Correlation');
    ylabel('Density');
    legend('All pairs, null', 'DMSO, null', 'Pert_id corrs', 'Location', 'NorthWest');
    xlim([-1 1]);
    title({sprintf('Distribution of instance correlations for %s, N perts = %d', u{ii}, numel(uu)); ...
        sprintf('Fraction of replicate pairs with q < 0.1 vs all: %0.3f; vs dmso: %0.3f', mean(pertq < 0.1), ...
            mean(dmsoq < 0.1))});
    print(gcf, '-dpng', '-r250', fullfile(args.outdir, sprintf('%s_instcorrs_bypert%s.png', u{ii}, ...
        ifelse(args.flatpc1, '_fpc1', ''))));
  end
end


function [ds, rds] = get_cpc6_data(args)
  % get default fields
  z = sig_info('{"cell_id":"cows"}');

  if strcmpi(args.cell_id, 'all')
    s = sig_info('{"sig_id":{"$regex":"^CPC006"}}', 'fields', vertcat(fieldnames(z), 'distil_id'));
    sc = struct2cell(s)';
    sc(:,cellstrfind(fieldnames(s), 'distil_id')) = mongolist2cell(sc(:,cellstrfind(fieldnames(s), 'distil_id')));

    rs = inst_info(vertcat(sc{:,cellstrfind(fieldnames(s), 'distil_id')}));
    rsc = struct2cell(rs)';

    ds = parse_gctx('/cmap/data/build/a2y13q1/modzs.gctx', 'rid', args.lms, 'cid', sc(:,1));
    rds = parse_gctx('/cmap/data/build/a2y13q1/zspc.gctx', 'rid', args.lms, 'cid', rsc(:,1));

    ds.cdesc = sc;
    ds.chd = fieldnames(s);
    ds.cdict = containers.Map(ds.chd, 1:numel(ds.chd));

    rds.cdesc = rsc;
    rds.chd = fieldnames(rs);
    rds.cdict = containers.Map(rds.chd, 1:numel(rds.chd));
  end

  if args.flatpc1
    ds = project_pc1(ds);
    rds = project_pc1(rds);
  end
end
