ccle = parse_gctx('/cmap/projects/cell_line_diversity/data/CCLE_BASELINE_RNASEQ_L1KFULL_RPKM_LOG2_LABELED_n1022x12450.gctx');

[p,n] = size(ccle.mat);

%% Mean variance scatter
figure;
scatter(mean(ccle.mat,2),std(ccle.mat,[],2))
grid on
xlabel('mean')
ylabel('std')
title(sprintf('Mean - Variance scatter of CCLE baseline expression\n %d genes',p))

%% All expression levels
figure
histogram(ccle.mat)
set(gca,'Yscale','log')
grid on
xlabel('Expression Level')
ylabel('Count')
title('All CCLE expression levels')

%% Sparsity pattern
figure;
spy(~ccle.mat)

%% Pairwise distances
ccle_pdist = parse_gctx('/cmap/projects/cell_line_diversity/analysis/baseline_gex_distance/CCLE_RNA_SEQ_DISTANCE_COSINE_n1022x1022.gctx');
blacklist = parse_grp('/cmap/projects/cell_line_diversity/data/cell_lines_blacklist.grp');
describe(tri2vec(ccle_pdist.mat))
