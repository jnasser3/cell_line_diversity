
%ds = parse_gctx('/cmap/projects/cell_line_diversity/data/cmap_training_data/cp_pruned_n1890x978.gctx');

corrs = fastcorr(ds.mat,'type','spearman');

[~,xform,x] = xform_corr_diversity(diag(corrs),tri2vec(corrs));

[fks,xks] = ksdensity(tri2vec(corrs),'npoints',1000);
plot(xks,fks,'-r','DisplayName','Null Correlations')
hold on
plot(xks,fks/max(fks),'-k','DisplayName','Null Correlations - Scaled')
hold on
plot(x,xform,'-b','DisplayName','Contribution to Diversity')
xlabel('Correlation')
title('Correlation Contribution to Diversity')
legend('show')
grid on
