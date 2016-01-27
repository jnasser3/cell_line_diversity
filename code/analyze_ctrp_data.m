ds = parse_gctx('/cmap/projects/cell_line_diversity/data/ctrp/ctrp_auc_n907x545.gctx');

%% How many NaN are there?
num_missing = nnz(isnan(ds.mat))
spy(isnan(ds.mat))
  xlabel('Cell Lines')
  ylabel('Compounds')
  title_str = sprintf('Missing AUC in CTRP data - %.2f percent missing', num_missing/numel(ds.mat));
  title(title_str)
 
%% Remove perts and cells with lots of missing data
  cell_cutoff_tol = .20; %cells with less than cutoff are included
  cp_cutoff_tol = .20; %cp with less than cutoff are included
  
  num_cp = size(ds.mat,1);
  num_cell = size(ds.mat,2);  
  cp_nan_count = sum(isnan(ds.mat),2);
  cp_nan_good_idx = find(cp_nan_count/num_cp < cp_cutoff_tol);
  cell_nan_count = sum(isnan(ds.mat),1);
  cell_nan_good_idx = find(cell_nan_count/num_cell < cell_cutoff_tol);

good_ds = ds_slice(ds,'ridx',cp_nan_good_idx,'cidx',cell_nan_good_idx);

%% Look at data for phase 1 core lines
phase1_core_lines = parse_grp('/cmap/projects/cell_line_diversity/data/phase1_core_lines.grp');
idx = find(ismember(ds.cdesc(:,ds.cdict('ccl_name')),phase1_core_lines))
ds_phase1_core = ds_slice(ds,'cidx',idx);

spy(ds_phase1_core.mat)


