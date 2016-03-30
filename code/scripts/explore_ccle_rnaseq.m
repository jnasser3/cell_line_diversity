ccle = parse_gctx('/cmap/projects/cell_line_diversity/data/CCLE_BASELINE_RNASEQ_L1KFULL_RPKM_LOG2_LABELED_n1022x12450.gctx');

[p,n] = size(ccle.mat);

%% Mean variance scatter
scatter(mean(ccle.mat,2),std(ccle.mat,[],2))
grid on
xlabel('mean')
ylabel('std')
title(sprintf('Mean - Variance scatter of CCLE baseline expression\n %d genes',p))

%% All expression levels
histogram(ccle.mat)
set(gca,'Yscale','log')
grid on
xlabel('Expression Level')
ylabel('Count')
title('All CCLE expression levels')

%% Sparsity pattern
spy(~ccle.mat)