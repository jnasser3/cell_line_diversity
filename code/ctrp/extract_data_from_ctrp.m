% Extract data from the ctrp into a gct struct. The output should be a 
% gct of perts x cells with the AUC measure

% Get meta info
cell_info = parse_tbl('/cmap/projects/cell_line_diversity/data/ctrp/v20.meta.per_cell_line.txt','outfmt','record','detect_numeric',false);
pert_info = parse_tbl('/cmap/projects/cell_line_diversity/data/ctrp/v20.meta.per_compound.txt','outfmt','record','detect_numeric',false);
experiment_info = parse_tbl('/cmap/projects/cell_line_diversity/data/ctrp/v20.meta.per_experiment.txt','outfmt','record','detect_numeric',false);

% Get curve info
curve_info = parse_tbl('/cmap/projects/cell_line_diversity/data/ctrp/v20.data.curves_post_qc_CLEAN_small_test.txt','outfmt','record');

% Get experiment ids (correspond to cell ids) and master_cpd_id. Form dimensions
unique_exp = unique([curve_info.experiment_id]);
unique_cpd = unique([curve_info.master_cpd_id]);

%Load data from flat file into matrix format
exp_map = containers.Map(unique_exp,1:numel(unique_exp));
cpd_map = containers.Map(unique_cpd,1:numel(unique_cpd));

data = NaN(numel(unique_cpd),numel(unique_exp));

column_name = 'apparent_ec50_umol';
file_suffix = 'ec50';

for ii = 1:numel(curve_info)
	   data(cpd_map(curve_info(ii).master_cpd_id),exp_map(curve_info(ii).experiment_id)) = curve_info(ii).(column_name);
end

%these are only temporary
temp_rid = unique_cpd';
temp_cid = unique_exp';

%ds = mkgctstruct(data,'rid',cellstr(int2str(temp_rid)),'cid',cellstr(int2str(temp_cid)));

%Make data set
ds = mkgctstruct(data,'rid',strtrim(cellstr(int2str(temp_rid))),'cid',strtrim(cellstr(int2str(temp_cid))));

%Annotate rows (compunds)
  ds = annotate_ds(ds,pert_info,'dim','row','keyfield','master_cpd_id')

%%Annotate columns (cell lines)
%Map experiment ids to cell ids
  exp_id_2_cell_id_map = containers.Map({experiment_info.experiment_id},{experiment_info.master_ccl_id});

%ds = annotate_ds(ds,experiment_info,'dim','column',keyfield,'experiment_id');

%change the cid in the ds to reflect cell ids
ds.cid = cellfun(@(s) exp_id_2_cell_id_map(s),ds.cid,'UniformOutput',false);

%now annotate based on cell id
ds = annotate_ds(ds,cell_info,'dim','column','keyfield','master_ccl_id','skipmissing',true);

%save the gctx
filename = strcat('/cmap/projects/cell_line_diversity/data/ctrp/ctrp_',file_suffix,'.gctx')
%mkgctx(filename,ds)
