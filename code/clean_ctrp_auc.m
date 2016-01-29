function ds_out = clean_ctrp_auc(ds,varargin)
%Cleans up the AUC data
%
%Input: ds - the data set
%       cp_cutoff - removes compounds profiled in less than this fraction of
%                   cells. Default: .2
%       cell_cutoff - removes cells profiled in less than this fraction of
%                   compounds. Default: .2
%       max_auc - set all auc's greater than max_auc to max_auc
%       min_auc - set all auc's less than min_auc to min_auc
%
%Output: ds_out - the cleaned ds

params = {'cp_cutoff',...
          'cell_cutoff',...
          'max_auc',...
          'min_auc'};
dflts = {.2,...
    .2,...
    20,...
    0};
args = parse_args(params,dflts,varargin{:});

%% Remove cells/cpds with low coverage
num_cp = size(ds.mat,1);
num_cell = size(ds.mat,2);  
cp_nan_count = sum(isnan(ds.mat),2);
cp_nan_good_idx = find(cp_nan_count/num_cp < args.cp_cutoff);
cell_nan_count = sum(isnan(ds.mat),1);
cell_nan_good_idx = find(cell_nan_count/num_cell < args.cell_cutoff);

ds_out = ds_slice(ds,'ridx',cp_nan_good_idx,'cidx',cell_nan_good_idx);

%% Set extreme auc values to floor/ceiling
low_idx = (ds_out.mat < args.min_auc);
ds_out.mat(low_idx) = args.min_auc;
high_idx = (ds_out.mat > args.max_auc);
ds_out.mat(high_idx) = args.max_auc;
