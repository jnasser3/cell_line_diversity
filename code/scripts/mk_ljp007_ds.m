%Annotate LJP79 dataset

ds = parse_gctx('/cmap/obelix/pod/custom/LJP/builds/cohort_name_LJP007-9_new_pc/modzs_n20122x22268.gctx');

%Subset to lm genes
lm_genes = mortar.common.Spaces.probe('lm').asCell;
ds = ds_slice(ds, 'rid', lm_genes);

%Annotate ds
annot = parse_tbl('/cmap/projects/M1/LJP/annot/LJP79_siginfo.txt');
ds = annotate_ds(ds, annot, 'dim', 'column', 'keyfield', 'sig_id');

%Add pert id, dose, time combo
combo = strcat(ds.cdesc(:,ds.cdict('pert_id')),'_',...
    ds.cdesc(:,ds.cdict('pert_idose')),'_',...
    ds.cdesc(:,ds.cdict('pert_itime')));

combo_annot = struct('sig_id',ds.cid,'pert_id_dose_time',combo);

ds = annotate_ds(ds,combo_annot,...
    'dim','column',...
    'keyfield','sig_id');