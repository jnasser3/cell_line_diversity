ds = parse_gctx('/cmap/projects/cell_line_diversity/data/ctrp_paper/ctrp_auc_paper_n664x481.gctx');

%% How many NaN are there?
num_missing = nnz(isnan(ds.mat));
spy(isnan(ds.mat))
  xlabel('Cell Lines')
  ylabel('Compounds')
  title_str = sprintf('Missing AUC in CTRP data - %.2f percent missing', num_missing/numel(ds.mat));
  title(title_str)
  
%% Look at data for phase 1 core lines
phase1_core_lines = parse_grp('/cmap/projects/cell_line_diversity/data/phase1_core_lines.grp');
idx = find(ismember(ds.cdesc(:,ds.cdict('cell_line_name')),phase1_core_lines));
ds_phase1_core = ds_slice(ds,'cidx',idx);

spy(ds_phase1_core.mat)

%% 


