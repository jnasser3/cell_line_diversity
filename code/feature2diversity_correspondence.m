function [rho, pval] = feature2diversity_correspondence(feature_ds,diversity_ds,varargin )
% Computes the correlation between the pairwise distances of cell lines in
% feature space and the pairwise ditances of cell lines in diversity space.
%
% Input: feature_ds: a cell_lines x cell_lines struct where each entry
%        represents the distance between the lines in feature space.
%        diversity_ds: a cell_lines x cell_lines struct where each entry
%        represents the distance between the lines in diversity space.
%        feature_cell_name_identifier: The feature_ds column identifier
%        corresponding the cell name. Must be either cid or
%        a member of feature_ds.chd. Default is cid.
%        diversity_cell_name_identifier: The diversity_ds column identifier
%        corresponding the cell name. Must be either cid or
%        a member of diversity_ds.chd. Default is cid.
%        corr_type: Measure of the correlation between feature and
%        diversity distances. Default is spearman.
%
% Output: rho: the correlation between the distances.
%         pval = pvalue of this correlation

params = {'feature_cell_name_identifier',...
          'diversity_cell_name_identifier',...
          'corr_type'};
dflts = {'cid',...
         'cid',...
         'spearman'};
args = parse_args(params,dflts,varargin{:});

% Get the ds from gctx if necessary
if exist(feature_ds,'file')
    feature_ds = parse_gctx(feature_ds);
end
if exist(diversity_ds,'file')
    diversity_ds = parse_gctx(diversity_ds);
end

%make the cid in each ds equal to the common identifier
if ~strcmp(args.feature_cell_name_identifier,'cid')
    feature_ds.cid = feature_ds.cdesc(:,feature_ds.cdict(args.feature_cell_name_identifier));
    feature_ds.rid = feature_ds.cid;
end
if ~strcmp(args.diversity_cell_name_identifier,'cid')
    diversity_ds.cid = diversity_ds.cdesc(:,diversity.cdict(args.diversity_cell_name_col));
    diversity_ds.rid = diversity_ds.cid;
end

%Subset each ds to common lines
%Note that ds_slice will order the structs alphabetically.
common_lines = intersect(feature_ds.cid,diversity_ds.cid);
feat_ds = ds_slice(feature_ds,'rid',common_lines,'cid',common_lines);
div_ds = ds_slice(diversity_ds,'rid',common_lines,'cid',common_lines);
fprintf('Feature ds contains %d cell lines \n',numel(feature_ds.cid));
fprintf('Diversity ds contains %d cell lines \n',numel(diversity_ds.cid));
fprintf('There are %d cell lines in common \n',numel(common_lines));

assert(isequal(feat_ds.cid,div_ds.cid),'ds cids not same order!')
assert(isequal(feat_ds.rid,div_ds.rid),'ds rids not same order!')
[rho, pval] = corr(vec_upper_triangle(feat_ds.mat),vec_upper_triangle(div_ds.mat),...
    'type','spearman',...
    'tail','right');
end

function v = vec_upper_triangle(A)
%returns the strictly upper triangle of A in vectorized form. A should be square

[n,m] = size(A);
assert(n==m, 'A is not square!')

%make the index vector
T = triu(ones(n),1);
T = T(:);

%index the matrix
temp = A(:);
v = temp(T == 1);

end
