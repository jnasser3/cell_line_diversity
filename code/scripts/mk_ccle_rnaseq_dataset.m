%% Load the original ccle rna seq data set
ds = parse_gctx('/cmap/projects/stacks/STK035_BASE_RNASEQ/STK035_BASE_RNASEQ_L1KFULL_LOG2_n1023x12450.gctx');

%% Remove the bad sample.
C = fastcorr(ds.mat,'type','pearson');
bad_idx = find((sum(C) == 0));
ds = ds_slice(ds,'cid',ds.cid(bad_idx),'exclude_cid',true);

%% Change cid to be the cell_id
original_id = ds.cid;
ds.cid = ds.cdesc(:,ds.cdict('cell_id'));
ds = annotate_ds(ds, struct('cid',ds.cid,'original_id', original_id), 'dim', 'column', 'keyfield', 'cid');

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
for ii = 1:numel(ds.cid)
    try
        str = ds.cdesc(ii,ds.cdict('original_id'));
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
tissue_type_struct = struct('cid',ds.cid,'tissue_type',tissue_type');
ds = annotate_ds(ds,...
    tissue_type_struct,...
    'dim','column',...
    'keyfield','cid');

%If we can't extract tissue from the id, look in the backup file
for jj = 1:numel(ds.cid)
    if ismember(ds.cdesc(jj,ds.cdict('tissue_type')),{'2','hT_DD','-666'})
        idx = find(strcmp({backup_annot.cell_id},ds.cid(jj)));
        ds.cdesc(jj,ds.cdict('tissue_type')) = {backup_annot(idx).tissue_type};
        
%         if ~strcmp({backup_annot(idx).tissue_type},'')
%             good_ds.cdesc(jj,good_ds.cdict('tissue_type')) = {backup_annot(idx).tissue_type};
%         end
    end
end

%% Make the gctx
mkgctx('/cmap/projects/cell_line_diversity/data/CCLE_BASELINE_RNASEQ_L1KFULL_RPKM_LOG2_LABELED.gctx',ds)