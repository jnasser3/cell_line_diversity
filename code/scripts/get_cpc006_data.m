%Get all signatures from the CPC006 plates.
cpc_info = sig_info('{"sig_id":{"$regex":"^CPC006"}}');

sigs = {cpc_info.sig_id};
lm_path = mortar.common.Spaces.probe_space('lm');
lm_genes = parse_grp(lm_path{1});
cpc_ds = parse_gctx('/cmap/data/build/a2y13q1/modzs.gctx',...
    'rid',lm_genes,...
    'cid',sigs);

%Annotate them
cpc_ds = annotate_ds(cpc_ds,cpc_info,...
    'dim','column',...
    'keyfield','sig_id');

%Add the pert_id_dose_time combo field.
combo = strcat(cpc_ds.cdesc(:,cpc_ds.cdict('pert_id')),'_',...
    cpc_ds.cdesc(:,cpc_ds.cdict('pert_idose')),'_',...
    cpc_ds.cdesc(:,cpc_ds.cdict('pert_itime')));

combo_annot = struct('sig_id',cpc_ds.cid,'pert_id_dose_time',combo);

cpc_ds = annotate_ds(cpc_ds,combo_annot,...
    'dim','column',...
    'keyfield','sig_id');

%save the gctx
mkgctx('/cmap/projects/cell_line_diversity/data/cmap_training_data/CPC006.gctx',cpc_ds)