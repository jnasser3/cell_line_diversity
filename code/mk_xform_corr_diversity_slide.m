
ds = parse_gctx('/cmap/projects/cell_line_diversity/data/cmap_training_data/cp_pruned_n1890x978.gctx');
idx = ismember(ds.cdesc(:,ds.cdict('cell_id')),{'HA1E','HCC515'});
corrs = fastcorr(ds.mat(:,idx),'type','spearman');

[matched_corrs, unmatched_corrs] = mat2vec_diag(corrs);
[~,xform,x] = xform_corr_diversity(matched_corrs,unmatched_corrs,...
    'cdf','fast_kde');

%%
figure;
[~,fks,xks,cks] = kde(unmatched_corrs,1024,[-1,1]);
plot(xks,fks,'-r','DisplayName','PDF of Null Correlations')
hold on
% plot(xks,fks/max(fks),'-k','DisplayName','Null Correlations - Scaled')
% hold on
% plot(xks,cks,'-k','DisplayName','Emprical CDF')
plot(x,xform,'-b','DisplayName','Contribution to Diversity')
xlim([-1 1])
ax = gca;
ax.XTick = -1:.2:1;
%ax.YTick = [0:.25:1 1:1:4];
xlabel('Correlation')
title('Correlation Contribution to Diversity')
legend('show')
grid on
