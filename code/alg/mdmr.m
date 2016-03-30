function [Fstat, pval, null0] = mdmr(dist, feat, varargin)
%[Fstat, pval] = mdmr(dist, feat, varargin)
%
%Runs Multivariance Distance Matrix Regression. Given an n x n matrix of
%pairwise distances between n objects and an m x n matrix of features for
%each of these objects. Tests whether variance in the pairwise distance
%matrix can be explained by variance in the features.
%
%Zapala, Schork PNAS 2006
%
%Input: dist - an n x n struct of pairwise distances between objects.
%              dist.mat(i,j) is the distance between cell lines i and j
%       feat - an m x n struct of features of these objects. feat.mat(i,j) 
%              is the value for feature i in cell line j. m must satisfy 1 < m <= n
%       cell_lines - a cell array or .grp of cell lines to include. Default
%              is to include all common lines b/w dist and feat.
%       feature_ids - a cell array or .grp of features to include. Must
%              match feat.rid. Default is all.
%       nperms - number of permutations to calculate significance. Default
%              1000
%
%Output:
%       Fstat - Test statistic
%       Pval - pvalue under permutation null
%       null0 - null distribution

pnames = {'cell_lines',...
          'feature_ids',...
          'nperms',...
          'dist_cell_field',...
          'feat_cell_field'};
dflts = {'',...
         '',...
         1000,...
         'cid',...
         'cid'};
args = parse_args(pnames, dflts, varargin{:});

% Get the ds from gctx if necessary
% if ischar(dist)     
%     dist = parse_gctx(dist);
% end
% if ischar(feat)   
%     feat = parse_gctx(feat);
% end
% 
% %Change id's if necessary
% [dist, feat] = change_id(dist,...
%                          feat,...
%                          args.dist_cell_field,...
%                          args.feat_cell_field);
                     
%% Subset ds to common cell lines
common_lines = intersect(dist.cid, feat.cid);
cell_lines = get_array_input(args.cell_lines,common_lines);
common_lines = intersect(common_lines, cell_lines);

feat = ds_slice(feat, 'cid', common_lines);
dist = ds_slice(dist, 'cid', common_lines, 'rid', common_lines);
assert(isequal(feat.cid,dist.cid),'cell lines not ordered equally')

%% Subset feature ds to features.
feature_ids = get_array_input(args.feature_ids,feat.rid);
feat = ds_slice(feat, 'rid', feature_ids);

%% Run MDMR
X = feat.mat';
D = dist.mat;
[n,m] = size(X);
assert(1 < m & m <= n, 'm must satisfy 1 < m <= n')

H = X*inv(X'*X)*X';
A = -.5 * (D .^ 2);
G = (eye(n) - (1/n)*ones(n))*A*(eye(n) - (1/n)*ones(n));

figure;
Fstat = compute_Fstat(G,H,n,m);
null0 = permute_Fstat(G,H,n,m,args.nperms);
pval = nnz(null0 > Fstat)/numel(null0);
histogram(null0)
hold on
vline(Fstat)
end

function [ds1, ds2] = change_id(ds1,ds2,field1,field2)

%make the cid in each ds equal to the common identifier
if ~strcmp(field1,'cid')
    ds1.cid = ds1.cdesc(:,ds1.cdict(field1));
end
if ~strcmp(field2,'cid')
    ds2.cid = ds2.cdesc(:,ds2.cdict(field2));
end

end

function permuted_Fstat = permute_Fstat(G,H,n,m,nperms)

permuted_Fstat = zeros(1,nperms);

for ii = 1:nperms
    idx = randperm(n);
    thisG = G(idx,idx);
    
    permuted_Fstat(ii) = compute_Fstat(thisG,H,n,m);
end

end

function F = compute_Fstat(G,H,n,m)

F = trace(H*G*H)/(m - 1) / trace((eye(n) - H)*G*(eye(n) - H)) / (n - m);

end