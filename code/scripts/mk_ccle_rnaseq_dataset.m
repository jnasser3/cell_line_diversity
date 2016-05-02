%% Load big CCLE data set
ccle_all = parse_gctx('/cmap/projects/stacks/STK035_BASE_RNASEQ/ccle_may12_2015/ccle_rpkm_n1019x52576.gctx');

%% Merge with Cmap lines
cmap_lines = parse_gct('/cmap/projects/stacks/STK035_BASE_RNASEQ/cmap/SCset1.rpkm.gct');
ccle_cmap = merge_two(ccle_all,cmap_lines);

%% Subset to only protein coding genes
protein_coding_ids = parse_grp('/cmap/projects/cell_line_diversity/data/rnaseq/protein_coding_ensembl_id.grp');
ccle = ds_slice(ccle_cmap,'rid',protein_coding_ids);

%% Load the original ccle rna seq data set
% ds = parse_gctx('/cmap/projects/stacks/STK035_BASE_RNASEQ/STK035_BASE_RNASEQ_L1KFULL_LOG2_n1023x12450.gctx');

%% Remove the bad sample.
% C = fastcorr(ds.mat,'type','pearson');
% bad_idx = find((sum(C) == 0));
% ds = ds_slice(ds,'cid',ds.cid(bad_idx),'exclude_cid',true);
bad_sample = {'fh_CH157MN_CENTRAL_NERVOUS_SYSTEM-Tumor'};
ccle = ds_slice(ccle,'cid',bad_sample,'exclude_cid',true);

%% Transform expression to Log2(RPKM + 1)
ccle.mat = log2(ccle.mat + 1);

%% Annotate
ccle = annotate_ds(ccle,'/cmap/projects/stacks/STK035_BASE_RNASEQ/colmeta_n1023.txt',...
    'dim','column',...
    'skipmissing',true);

%% Get the cell_id right for cmap lines
keyset = {'CD34pos','HA1E_1','NEU_1','NPC_1'};
valset = {'CD34','HA1E','NEU','NPC'};
id_map = containers.Map(keyset,valset);
for ii = 1:4
    idx = strcmp(ccle.cid,keyset(ii));
    ccle.cdesc(idx,ccle.cdict('cell_id')) = {id_map(keyset{ii})};
end

%% Change cid to be the cell_id
original_id = ccle.cid;
ccle.cid = ccle.cdesc(:,ccle.cdict('cell_id'));
ccle = annotate_ds(ccle, struct('cid',ccle.cid,'original_id', original_id), 'dim', 'column', 'keyfield', 'cid');

%% Get and remove column data from ds
% we need to do this because we do not have prism annotations for all
% samples
% metastruct = meta2struct(ds, 'dim', 'column');
% ds = ds_delete_meta(ds, 'column', ds.chd);

%% Add PRISM annotations
%
%**ISSUE** Stripped cell line name is not unique (TT)
%
% prism_annot = parse_tbl('/cmap/projects/cell_line_diversity/data/PRISM_database.tsv','outfmt','record');
% 
% %subset to common ids
% common_ids = intersect({prism_annot.strippedname},good_ds.cid);
% prism_annot = prism_annot(ismember({prism_annot.strippedname},common_ids));
% 
% %Add prism annotations
% good_ds = annotate_ds(good_ds,prism_annot,...
%     'dim','column',...
%     'keyfield','strippedname',...
%     'skipmissing',true);

%The 'strippedname' field in prism_annot matches the 'cell_id' field in ds.cdesc
%So annotate the prism annot with the appropriate cids. Then use annotate_ds
% name2id_map = containers.Map({metastruct.cell_id},{metastruct.id});
% ids_for_prism_table = {};
% for ii = 1:numel({prism_annot.strippedname})
%     prism_annot(ii).temp_id = name2id_map(prism_annot(ii).strippedname);
% end

%% Get cell line annotations
% cell_annot = parse_tbl('/cmap/projects/cell_line_diversity/data/cell_line_metadata.tsv','outfmt','record');
% 
% %Extract the tissue type from primary_site and add to annotations.
% n = numel(cell_annot);
% for ii = 1:n
%     sites = strsplit(cell_annot(ii).primary_site);
%     tissue_type(ii) = sites(1); 
% end
% [cell_annot.tissue_type] = deal(tissue_type);
% good_ds = annotate_ds(good_ds,cell_annot,'keyfield','stripped_cell_line_name','skipmissing',true);

%% Add back in column meta
% ds = annotate_ds(ds,metastruct);

%% annotate tissue type
tissue_type = {};
backup_annot = parse_tbl('/cmap/projects/cell_line_diversity/data/ccle_lineage_backup_annot.txt','outfmt','record');
for ii = 1:numel(ccle.cid)
    try
        str = ccle.cdesc(ii,ccle.cdict('original_id'));
        underbar_loc = cell2mat(strfind(str,'_'));
        t1 = underbar_loc(2);
        dash_loc = cell2mat(strfind(str,'-'));
        t2 = dash_loc(end);
        str = cell2str(str);
        tissue_type(ii) = {str((t1+1):(t2-1))};
    catch
        tissue_type(ii) = {'-666'};
    end
    
end
tissue_type_struct = struct('cid',ccle.cid,'tissue_type',tissue_type');
ccle = annotate_ds(ccle,...
    tissue_type_struct,...
    'dim','column',...
    'keyfield','cid');

%If we can't extract tissue from the id, look in the backup file
for jj = 1:numel(ccle.cid)
    if ismember(ccle.cdesc(jj,ccle.cdict('tissue_type')),{'2','hT_DD','-666'})
        idx = find(strcmp({backup_annot.cell_id},ccle.cid(jj)));
        ccle.cdesc(jj,ccle.cdict('tissue_type')) = {backup_annot(idx).tissue_type};
        
%         if ~strcmp({backup_annot(idx).tissue_type},'')
%             good_ds.cdesc(jj,good_ds.cdict('tissue_type')) = {backup_annot(idx).tissue_type};
%         end
    end
end

%% Make the gctx
mkgctx('/cmap/projects/cell_line_diversity/data/rnaseq/CCLE_BASELINE_RNASEQ_PROTEIN_CODING_RPKM_LOG2_LABELED.gctx',ccle)