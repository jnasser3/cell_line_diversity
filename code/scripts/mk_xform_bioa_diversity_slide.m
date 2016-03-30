%ds = parse_gctx('/cmap/projects/cell_line_diversity/data/cmap_training_data/cp_pruned_n1890x978.gctx');

load /cmap/projects/cell_line_diversity/data/ts_cp_ccq75.mat

bioa = cell2mat(ds.cdesc(:,ds.cdict('distil_cc_q75')));
bioa(bioa <= 0) = [];

[~,xform,x] = xform_bioa_diversity(bioa,ts_cp_ccq75);

[ef,ex] = ecdf(ts_cp_ccq75);
plot(ex,ef,'DisplayName','Background Reproducibility')
hold on
plot(x,xform,'DisplayName','Contribution to Diversity')
xlabel('Replicate Reproducibility')
title('Reproducibility Contribution to Diversity')
legend('show')
grid on
