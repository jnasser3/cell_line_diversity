%Investigate this data set for outliers

ds1 = parse_gctx('/cmap/projects/stacks/STK035_BASE_RNASEQ/STK035_BASE_RNASEQ_L1KFULL_LOG2_n1023x12450.gctx');

%compute euclidean norm for each cell lines gex
stats_table = table(ds1.cid,sqrt(sum(ds1.mat .^ 2))',sum(ds1.mat < 0)',...
    'VariableNames',{'Cell','EuclideanNorm','NumNegative'});

C = fastcorr(ds1.mat,'type','pearson');
bad_idx = find((sum(C) == 0));

good_ds = ds_slice(ds1,'cid',ds1.cid(bad_idx),'exclude_cid',true)



%%THis is all useless, forgot to log2 beforehand
% % % % %% There are two wonky lines
% % % % % 'CH157MN' - has all negative entries
% % % % % 'KMS27' - has some extremely large values
% % % % CH157MN_idx = find(strcmp(ds1.cid,'fh_CH157MN_CENTRAL_NERVOUS_SYSTEM-Tumor'));
% % % % KMS27_idx = find(strcmp(ds1.cid,'fh_KMS27_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE-Tumor'));
% % % 
% % % %% Pick a cell line and compare to see if ds1 is really log2'd
% % % % idx1 = find(strcmp(ds1.cid,'fh_ZR7530_BREAST-Tumor'));
% % % % idx2 = find(strcmp(ds2.cid,'ZR7530'));
% % % % lm_path = mortar.common.Spaces.probe('lm');
% % % % lm_genes = parse_grp(lm_path);
% % % % ds1_test = ds_slice(ds1,'rid',lm_genes,'cidx',idx1,'ignore_missing',true);
% % % % ds2_test = ds_slice(ds2,'rid',ds1_test.rid,'cidx',idx2,'ignore_missing',true);
