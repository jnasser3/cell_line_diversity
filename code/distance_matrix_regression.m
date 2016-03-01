function mdl = distance_matrix_regression(feat,resp,varargin)
% simple_distance_matrix_regression(ds_dep,ds_ind,varargin)
%
% Runs a simple linear regression of one distance matrix against another.
%
% Input: feat - feature matrix
%        resp - response vector

params = {};
dflts = {};
args = parse_args(params,dflts,varargin{:});

%% Make features and responses into tables
feat_table = array2table(feat.mat,'RowNames',feat.rid,'VariableNames',feat.cid);
resp_table = array2table(resp.mat,'RowNames',resp.rid,'VariableNames',resp.cid);

lm_table = [feat_table resp_table];

mdl = fitlm(lm_table);

end