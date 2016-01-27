ds = parse_gctx('/cmap/projects/cell_line_diversity/data/ctrp/ctrp_auc_n907x545.gctx');

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

%% Use matrix completion to fill in the NaNs
%small_test = good_ds.mat(1:16,1:16);
missing_mat = good_ds.mat;

[num_cp,num_cell] = size(missing_mat);
mu = .001;
omega = sort(find(~isnan(missing_mat))');
obs = double(missing_mat(omega));

completed_mat = solver_sNuclearBP( {num_cp,num_cell,omega}, obs, mu);


