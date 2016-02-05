%Load the original ccle rna seq data set
ds = parse_gctx('/cmap/projects/stacks/STK035_BASE_RNASEQ/STK035_BASE_RNASEQ_L1KFULL_LOG2_n1023x12450.gctx');

%Remove the bad sample.
C = fastcorr(ds.mat,'type','pearson');
bad_idx = find((sum(C) == 0));
good_ds = ds_slice(ds,'cid',ds.cid(bad_idx),'exclude_cid',true);

%annotate tissue type
tissue_type = {};
for ii = 1:numel(good_ds.cid)
    try
        str = good_ds.cid{ii};
        underbar_loc = strfind(str,'_');
        t1 = underbar_loc(2);
        dash_loc = strfind(str,'-');
        t2 = dash_loc(end);
        tissue_type(ii) = {str((t1+1):(t2-1))};
    catch
        tissue_type(ii) = {'-666'};
    end
end
tissue_type_struct = struct('cid',good_ds.cid,'tissue_type',tissue_type');
good_ds = annotate_ds(good_ds,tissue_type_struct,...
    'dim','column',...
    'keyfield','cid');

%make the gctx
mkgctx('/cmap/projects/cell_line_diversity/data/CCLE_BASELINE_RNASEQ_L1KFULL_RPKM_LOG2.gctx',good_ds)