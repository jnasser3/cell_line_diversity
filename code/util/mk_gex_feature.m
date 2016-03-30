function [feat,gex_pca,gex_pca_feat_subset] = mk_gex_feature(gex, varargin)
% Given a genes x cells struct of baseline gex data of p genes and N cells,
% returns a NC2 x m struct where the rows represent pairs of cells and the
% columns are a set of (transformed/engineered) baseline gex features.
%
% Inputs:
%        gex: genes x samples expression dataset
%        feature_name: Name prefix for the feature. Default is GEX
%        ids_for_pca: cids in the ds to be used in PCA calculation
%        ids_for_features: cids in the ds to return as pairwise features
%        num_features: number of PCs to compute
%        method: method for combining PC coefficients
%
% Outputs: feat

params = {'feature_name',...
          'ids_for_pca',...
          'ids_for_features',...
          'num_features',...
          'method'};
dflts = {'GEX',...
         '',...
         '',...
         10,...
         'abs_subtract'};
args = parse_args(params,dflts,varargin{:});

%% Make the cid to the cell id
%gex.cid = gex.cdesc(:,gex.cdict('cell_id'));
gex = ds_slice(gex, 'cid', args.ids_for_pca, 'ignore_missing',true);

%% Run pca and annotate result
gex_tr = transpose_gct(gex);
[~, pc_score] = pca(gex_tr.mat);
PCA_labels = gen_labels(size(pc_score, 2), 'prefix', strcat(args.feature_name,'_PCA_'));
row_meta = gctmeta(gex_tr, 'row');

gex_pca = mkgctstruct(pc_score, 'rid', gex_tr.rid, 'cid', PCA_labels);
gex_pca = annotate_ds(gex_pca, row_meta, 'dim', 'row');

gex_pca_feat_subset = ds_slice(gex_pca, 'rid', args.ids_for_features,...
    'ignore_missing',true);

%% Mk feature struct
p = args.num_features;
n = numel(gex_pca_feat_subset.rid);
mat = zeros(nchoosek(n,2),p);
rid = cell(nchoosek(n,2),1);

kk = 1;
for ii = 1:n
    for jj = (ii+1):n
        cell1 = gex_pca_feat_subset.rid(ii);
        cell2 = gex_pca_feat_subset.rid(jj);
        idx1 = strcmp(gex_pca_feat_subset.rid,cell1);
        idx2 = strcmp(gex_pca_feat_subset.rid,cell2);
        
        switch args.method
            case 'abs_subtract'
                mat(kk,:) = abs(gex_pca_feat_subset.mat(idx1,1:p) - gex_pca_feat_subset.mat(idx2,1:p));
        end
        rid(kk) = strcat(cell1,'_',cell2);
        kk = kk + 1;
    end
end

feat = mkgctstruct(mat, 'rid', rid, 'cid', gex_pca.cid(1:p));
end

