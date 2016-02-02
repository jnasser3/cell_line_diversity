function out_ds = pairwise_auc_distance(ds,varargin)
% Computes the pairwise distance matrix from auc pharmacological profiling data
% 
% Inputs:
%        ds - a data structure where the rows of ds.mat are the perts and the columns are cell lines
%        method - method to use to compute distances. 
% Outputs:
%        out_ds - A struct. out_ds.mat will have the pairwise distances

params = {'method'};
dflts = {'euclidean'};
args = parse_args(params, dflts, varargin{:});



switch args.method
    case 'euclidean'
        D = apply_euclidean(ds.mat);
    case 'transformed_fisher'
        D = apply_transformed_fisher(ds.mat);
end

out_ds = mkgctstruct(D,'rid',ds.cid,...
        'cid',ds.cid,...
        'cdesc',ds.cdesc,...
        'chd',ds.chd,...
        'cdict',ds.cdict);

end

function D = apply_euclidean(A)

[~,num_cell] = size(A);
D = zeros(num_cell);

for ii = 1:num_cell
       for jj = ii:num_cell
              nan_idx = isnan(A(:,ii)) | isnan(A(:,jj));
              D(ii,jj) = norm(A(~nan_idx,ii) - A(~nan_idx,jj),2) / nnz(~nan_idx);
       end
end

%fill in lower triangle
D = D + D';
end

function D = apply_transformed_fisher(A)
%see Seashore-Ludlow et al 2015 and reference (16) therein.
%
%How it works:
%   1. Compute pearson on mutual support
%   2. Transform to Fisher score
%   3. Z-score the fisher score
%   4. Tune alpha
%   5. Apply double sigmoid transform.

[~,num_cell] = size(A);
rho = zeros(num_cell);
F = zeros(num_cell);
Z = zeros(num_cell);

for ii = 1:num_cell
       for jj = (ii+1):num_cell
              nan_idx = isnan(A(:,ii)) | isnan(A(:,jj));
              rho(ii,jj) = corr(A(~nan_idx,ii),A(~nan_idx,jj),'type','pearson');
              F(ii,jj) = log((1+rho(ii,jj))/(1-rho(ii,jj)))/2;
              Z(ii,jj) = F(ii,jj)*sqrt(numel(~nan_idx) - 3);
       end
end

%Compute alpha
alpha = tune_alpha(Z, 3, -.95);

%compute Distances using double sigmoid transform
D = arrayfun(@(z) 1 + ((z/alpha)^K) / sqrt(1 + (z/alpha)^(2*K)), Z);

%fill in lower triangle
D = D + D';
end

function alpha = tune_alpha(Z, K, min_lim)

pmin = normcdf(min(min(tri2vec(Z))));
plim = 10^(-1*(min_lim*(-1*log10(pmin) + log10(1+min_lim)) - log10(1+min_lim)));
zlim = norminv(plim);
alpha = abs(zlim / (((min_lim - 1)^2/(2 - min_lim))^(1/(2*K))));

end