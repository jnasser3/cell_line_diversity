function mdl = simple_distance_matrix_regression(ds_ind,ds_dep,varargin)
% simple_distance_matrix_regression(ds_dep,ds_ind,varargin)
%
% Runs a simple linear regression of one distance matrix against another.
%
% Input: ds_dep - distance we are trying to predict
%        ds_ind - Independant distance (Feature)

params = {'holdout_pct',...
          'ind_id_field',...
          'dep_id_field',...
          'ids'};
dflts = {0,...
        'cid',...
        'cid',...
        ''};
args = parse_args(params,dflts,varargin{:});

%% Subset each ds to common samples
[ds_ind,ds_dep,common_objects] = subset_ds_common_ids(ds_ind,ds_dep,...
                                                      'ids',args.ids,...
                                                      'ds1_id_field','cell_id');
n = numel(common_objects);
assert(isequal(ds_ind.cid,ds_dep.cid),'ds cids not same order!')
assert(isequal(ds_ind.rid,ds_dep.rid),'ds rids not same order!')

%% Remove a subset of the samples as a holdout set
holdout_num = floor(args.holdout_pct * n);
holdout_idx = randperm(n,holdout_num);
holdout_ids = ds_ind.cid(holdout_idx);

ds_ind_train = ds_slice(ds_ind, 'cid', holdout_ids, 'exclude_cid', true, ...
                          'rid', holdout_ids, 'exculde_rid', true);
ds_dep_train = ds_slice(ds_dep, 'cid', holdout_ids, 'exclude_cid', true, ...
                          'rid', holdout_ids, 'exculde_rid', true);
ds_ind_test = ds_slice(ds_ind, 'cid', holdout_ids, 'rid', holdout_ids); 
ds_dep_test = ds_slice(ds_dep, 'cid', holdout_ids, 'rid', holdout_ids);
                      
assert(isequal(ds_ind_train.cid,ds_dep_train.cid),'training cids not same order!')
assert(isequal(ds_ind_train.rid,ds_dep_train.rid),'training rids not same order!')
assert(isequal(ds_ind_test.cid,ds_dep_test.cid),'test cids not same order!')
assert(isequal(ds_ind_test.rid,ds_dep_test.rid),'test rids not same order!')

% Run the regression
ind_train = tri2vec(ds_ind_train.mat);
dep_train = tri2vec(ds_dep_train.mat);      
mdl = fitlm(ind_train,dep_train);

% Check on hold out
ind_test = tri2vec(ds_ind_test.mat);
dep_pred = feval(mdl,ind_test);
dep_test = tri2vec(ds_dep_test.mat);

end