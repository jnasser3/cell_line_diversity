%% Get PCs of CCLE baseline expression
gex_tr = transpose_gct(ccle);
[~, pc_score] = pca(gex_tr.mat);
PCA_labels = gen_labels(size(pc_score, 2), 'prefix', strcat('GEX','_PCA_'));
row_meta = gctmeta(gex_tr, 'row');

gex_pca = mkgctstruct(pc_score, 'rid', gex_tr.rid, 'cid', PCA_labels);
gex_pca = annotate_ds(gex_pca, row_meta, 'dim', 'row');
gex_pca = transpose_gct(gex_pca);

%% Get some pcs
k = 5;
feat = ds_slice(gex_pca,'ridx',1:k);

%% POSCON: Test regression of ccle cosine distance on the PCs
[fstat, pval, ~] = mdmr(ccle_cosine,feat);

%% Negative control: Test regression of ccle cosine distance on random numbers
% rand_struct = mkgctstruct(randn(2,numel(ccle_cosine.rid)),'cid',ccle_cosine.rid,'rid',{'random1','random2'});
% [fstat, pval, null0] = mdmr(ccle_cosine,rand_struct,'nperms',100)

%% Test regression of cpc006 pdi on the PCs
[fstat2, pval2, ~] = mdmr(cpc006_pdi,feat)

%% Test regression of core_ts pdi on the PCs
% [fstat3, pval3, ~] = mdmr(core_ts_pdi,feat)