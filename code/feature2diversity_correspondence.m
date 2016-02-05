function [rho, pval, feat_vec, div_vec] = feature2diversity_correspondence(feature_ds,diversity_ds,varargin )
% Computes the correlation between the pairwise distances of cell lines in
% feature space and the pairwise ditances of cell lines in diversity space.
%
% Input: feature_ds: a cell_lines x cell_lines struct where each entry
%                   represents the distance between the lines in feature space.
%        diversity_ds: a cell_lines x cell_lines struct where each entry
%                   represents the distance between the lines in diversity space.
%        feature_cell_name_identifier: The feature_ds column identifier
%                   corresponding the cell name. Must be either cid or
%                   a member of feature_ds.chd. Default is cid.
%        diversity_cell_name_identifier: The diversity_ds column identifier
%                   corresponding the cell name. Must be either cid or
%                   a member of diversity_ds.chd. Default is cid.
%        corr_type: Measure of the correlation between feature and
%                   diversity distances. Default is spearman.
%
% Output: rho: the correlation between the distances.
%         pval: pvalue of this correlation
%         feat_vec: vector of distances in feature space
%         div_vec: vector of distances in diversity space

params = {'feature_cell_name_identifier',...
          'diversity_cell_name_identifier',...
          'corr_type'};
dflts = {'cid',...
         'cid',...
         'Spearman'};
args = parse_args(params,dflts,varargin{:});

% Get the ds from gctx if necessary
if isstr(feature_ds)     %exist(feature_ds,'file')
    feature_ds = parse_gctx(feature_ds);
end
if isstr(diversity_ds)   %exist(diversity_ds,'file')
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

feat_vec = tri2vec(feat_ds.mat);
div_vec = tri2vec(div_ds.mat);
[rho, pval] = corr(feat_vec,div_vec,...
    'type',args.corr_type,...
    'tail','right');

%Plotting
scatter(feat_vec,div_vec, 250, '.')
xlabel('Distance in feature space')
ylabel('Distance in diversity space')
grid on;
title_str = sprintf('N = %d cell lines, %d comparisons \n%s correlation is %.2f, pval is %0.4f', ...
    numel(common_lines),nchoosek(numel(common_lines),2),args.corr_type,rho,pval);
title(title_str)

end

