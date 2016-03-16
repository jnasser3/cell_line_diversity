
ds = parse_gctx('/cmap/projects/cell_line_diversity/data/cmap_training_data/cp_pruned_n1890x978.gctx');
idx = ismember(ds.cdesc(:,ds.cdict('cell_id')),{'HA1E','HCC515'});
corrs = fastcorr(ds.mat(:,idx),'type','spearman');

[~,xform,x] = xform_corr_diversity(diag(corrs),tri2vec(corrs),...
    'cdf','fast_kde');

figure;
[~,fks,xks,cks] = kde(tri2vec(corrs),1024,[-1,1]);
plot(xks,fks,'-r','DisplayName','PDF of Null Correlations')
hold on
% plot(xks,fks/max(fks),'-k','DisplayName','Null Correlations - Scaled')
% hold on
% plot(xks,cks,'-k','DisplayName','Emprical CDF')
plot(x,xform,'-b','DisplayName','Contribution to Diversity')
xlabel('Correlation')
title('Correlation Contribution to Diversity')
legend('show')
grid on
