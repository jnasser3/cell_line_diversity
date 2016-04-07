%Get all A2 cp signatures
sigs = mongo_info('sig_info','{pert_type:"trt_cp"}');

%Filter for core lines


%get A2 signatures in lm space
lm = mortar.common.Spaces.probe('lm').asCell;

ds = parse_gctx('/cmap/data/build/a2y13q1/modzs.gctx',...
    'cid',{sigs.sig_id},...
    'rid',lm);

%annotate ds
ds = annotate_ds(ds,sigs,'keyfield','sig_id');

%annotate pert_id_dose_time
ds = annotate_ds_concatenated_field(ds,...
    'column',...
    {'pert_id','pert_idose','pert_itime'},...
    'pert_id_dose_time');

%Get core lines only
lines = parse_grp('/cmap/projects/cell_line_diversity/data/spaces/phase1_core_lines.grp');
ds_core = select_training_data_function(ds,...
    'cell_lines',lines,...
    'match_field','pert_id_dose_time');