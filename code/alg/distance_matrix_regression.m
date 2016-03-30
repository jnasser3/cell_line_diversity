function [mdl,lm_table] = distance_matrix_regression(feat,resp,varargin)
% simple_distance_matrix_regression(ds_dep,ds_ind,varargin)
%
% Runs a linear regression of resp against feat. 
%
% Input: feat - feature matrix
%        resp - response vector

params = {};
dflts = {};
args = parse_args(params,dflts,varargin{:});

%% Make features and responses into tables
%[feat,resp] = order_alphabetically(feat,resp);

feat_table = array2table(feat.mat,'RowNames',feat.rid,'VariableNames',feat.cid);
% resp_table = array2table(resp.mat,'RowNames',resp.rid,'VariableNames',resp.cid);

lm_table = horzcat(feat_table,resp);

%assert(isequal(feat.rid,resp.rid),'Feature and response ds not ordered identically!');
mdl = fitlm(lm_table);

end

function [ds1,ds2] = order_alphabetically(ds1,ds2)

    [rid1,idx1] = sort(ds1.rid);
    [rid2,idx2] = sort(ds2.rid);
    
    ds1.rid = rid1;
    ds2.rid = rid2;
    ds1.mat = ds1.mat(idx1,:);
    ds2.mat = ds2.mat(idx2,:);


end