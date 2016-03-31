ccle = parse_gctx('/cmap/projects/cell_line_diversity/data/CCLE_BASELINE_RNASEQ_L1KFULL_RPKM_LOG2_LABELED_n1022x12450.gctx');

[p,n] = size(ccle.mat);

%%
% Sanity check to confirm mahalanobis distance calculations
% Since CCLE has more genes than samples, first PCA down to a manageable
% number of dimensions. Confirm that computing mahalanobis distances in
% this space is equivalent to computing euclidean distances in the presence
% of the covaraince matrix. C must be diagonal in this case since the
% results of the PCA ensure we have orthogonal features.
num_pcs = 10;
num_samples = 1000;

test_mat = ccle.mat(:,1:num_samples);

[COEFF, SCORE] = pca(test_mat');

mat = SCORE(:,1:num_pcs);

mahal_dist = squareform(pdist(mat,'mahalanobis'));

C = diag(diag(cov(mat)));
xform_SCORE = inv(sqrt(C))*mat';

mahal_dist2 = squareform(pdist(xform_SCORE','euclidean'));

norm(mahal_dist - mahal_dist2,'fro')