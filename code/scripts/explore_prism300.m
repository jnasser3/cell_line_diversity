ds = parse_gctx('/cmap/projects/repurposing/rnwork/PREP_300/PREP_ZS_n4992x300.gctx');

%% Find cells with lots of Nan's
missing = sum(isnan(ds.mat),2);
missing_table = array2table(missing,'VariableNames',{'Num_NAN'},'RowNames',ds.rdesc(:,ds.rdict('cell_id')));
missing_table = sortrows(missing_table);

%Remove cells with lots of Nans
cutoff = 100;
cutoff_idx = (missing > cutoff);
ds = ds_slice(ds,'rid',ds.rid(~cutoff_idx));

%% Compute pairwise spearmans
%Ignore mismatch dimensionality
mat = 1 - corr(ds.mat',...
    'type','spearman',...
    'rows','pairwise');
cell_meta = ds_get_meta(ds, 'row', ds.rhd);
prism_pdist = mkgctstruct(mat,'rid',ds.rid,'cid',ds.rid);
prism_pdist = ds_add_meta(prism_pdist,'row',ds.rhd,cell_meta);
prism_pdist = ds_add_meta(prism_pdist,'column',ds.rhd,cell_meta);

%% tSNE
%run and annotate tsne
perplexity = 10;
outdir = '/cmap/projects/cell_line_diversity/analysis/prism';
outfile = 'prism_tsne';
tsne = tsne_d(prism_pdist.mat,[],2,perplexity);

tsne_ds = mkgctstruct(tsne,'rid',prism_pdist.rid,'cid',{'TSNE1','TSNE2'});
annot = meta2struct(prism_pdist,'dim','row');
tsne_ds = annotate_ds(tsne_ds,annot,'dim','row');

%Tableauify
temp1 = table(tsne_ds.mat(:,1),tsne_ds.mat(:,2),tsne_ds.rid,...
    'VariableNames',{'TSNE_1','TSNE_2','cell_name'});
temp2 = cell2table(tsne_ds.rdesc,'VariableNames',tsne_ds.rhd);
tsne_table = [temp1 temp2];

%Write table
writetable(tsne_table,fullfile(outdir,outfile),...
    'FileType','text',...
    'Delimiter','\t');

%% Get A2 signatures for the prism compounds
% prism_cp = unique(ds.cdesc(:,ds.cdict('pert_id')));
% sigs = mongo_info('sig_info',prism_cp,'query_field','pert_id');
% 
% %remove empty sigs
% empty_idx = cellfun(@isempty,{sigs.sig_id});
% sigs(empty_idx) = [];

%% Correspondence with cpc006
cpc006 = parse_gctx('/cmap/projects/cell_line_diversity/data/cmap_training_data/CPC006_n22098x978.gctx');

%get common lines
cpc006_lines = unique(cpc006.cdesc(:,cpc006.cdict('cell_id')));
prism_lines = unique(ds.rdesc(:,ds.rdict('cell_id')));
common_lines = intersect(cpc006_lines, prism_lines);

%get common perts
cpc006_cp = unique(cpc006.cdesc(:,cpc006.cdict('pert_id')));
prism_cp = unique(ds.cdesc(:,ds.cdict('pert_id')));
common_cp = intersect(cpc006_cp, prism_cp);

cpc006_pdi_for_prism = compute_pairwise_diversity_index(cpc006,...
    'prune_ds',true,...
    'cell_lines',common_lines,...
    'pert_ids',common_cp);

%Subset prism to common cp
prism_for_cpc006 = ds;
cp_idx = ismember(prism_for_cpc006.cdesc(:,prism_for_cpc006.cdict('pert_id')),common_cp);
prism_for_cpc006 = ds_slice(prism_for_cpc006,'cid',prism_for_cpc006.cid(cp_idx));
mat = 1 - corr(prism_for_cpc006.mat',...
    'type','spearman',...
    'rows','pairwise');
cell_meta = ds_get_meta(prism_for_cpc006, 'row', prism_for_cpc006.rhd);
prism_pdist_for_cpc006 = mkgctstruct(mat,'rid',prism_for_cpc006.rid,'cid',prism_for_cpc006.rid);
prism_pdist_for_cpc006 = ds_add_meta(prism_pdist_for_cpc006,'row',prism_for_cpc006.rhd,cell_meta);
prism_pdist_for_cpc006 = ds_add_meta(prism_pdist_for_cpc006,'column',prism_for_cpc006.rhd,cell_meta);

distance_matrix_correspondence(cpc006_pdi_for_prism,prism_pdist,...
    'ds2_id_field','cell_id','nperms',1000,'make_plot',true)



