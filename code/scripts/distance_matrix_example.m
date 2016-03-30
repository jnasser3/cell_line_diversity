
% The goal is to answer the question:
% can we predict pairwise diversity using gene expression. Do this by
% regressing pairwise diversity index on PCs of baseline expression


%% Get, flatten pairwise distances and make into a table for fitlm
pdi = parse_gctx('/cmap/projects/cell_line_diversity/analysis/pairwise_diversity_index/cpc006/rel_bioa_corr/pdi_n51x51.gctx');

%% Get baseline expression
ccle_gex = parse_gctx('/cmap/projects/cell_line_diversity/data/CCLE_BASELINE_RNASEQ_L1KFULL_RPKM_LOG2_LABELED_n1022x12450.gctx');

%% Make repsone table
ccle_lines = ccle_gex.cdesc(:,ccle_gex.cdict('cell_id'));
common_lines = intersect(pdi.cid,ccle_lines);
resp = flatten_pw_distance_matrix(pdi,'ids',common_lines);
resp_table = array2table(resp.mat,'RowNames',resp.rid,'VariableNames',resp.cid);


%% Mk baseline gex into features
% blacklist = parse_grp('/cmap/projects/cell_line_diversity/data/cell_lines_blacklist.grp');
% blacklist_cell_ids = ccle_lines(ismember(ccle_gex.cid,blacklist));
% ccle_no_suspension = setdiff(ccle_lines, blacklist_cell_ids);
% common_lines = setdiff(common_lines,blacklist_cell_ids);
gex_feat = mk_gex_feature(ccle_gex,...
    'feature_name','GEX',...
    'ids_for_pca',ccle_lines,...
    'ids_for_features',common_lines,...
    'num_features',10);

% %% Get Copy number data
% ccle_cn = parse_gctx('/cmap/data/vdb/cline/ccle_copynumber_bygene_20120929_n899x11774.gctx');
% cn_feat = mk_gex_feature(ccle_cn,...
%     'feature_name','COPYNUM',...
%     'ids_for_pca',ccle_no_suspension,...
%     'ids_for_features',common_lines,...
%     'num_features',5);
% 
% %% Combine features
% all_feat = combine_features({gex_feat,cn_feat});

%% Run the regression
all_feat = gex_feat;
[mdl,tbl] = distance_matrix_regression(all_feat,resp_table);