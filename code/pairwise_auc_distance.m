function out_ds = pairwise_auc_distance(ds,varargin)
% Computes the pairwise distance matrix from auc pharmacological profiling data
% 
% Inputs:
%        ds - a data structure where the rows of ds.mat are the perts and the columns are cell lines
%        method - method to use to compute distances
% Outputs:
%        out_ds - A struct. out_ds.mat will have the pairwise distances

params = {'method'};
dflts = {'euclidean'};
args = parse_args(params, dflts, varargin{:});

[~,num_cell] = size(ds.mat);
D = zeros(num_cell);

switch args.method
    case 'euclidean'
        for ii = 1:num_cell
               for jj = ii:num_cell
                      nan_idx = isnan(ds.mat(:,ii)) | isnan(ds.mat(:,jj));
                      D(ii,jj) = norm(ds.mat(~nan_idx,ii) - ds.mat(~nan_idx,jj),2) / nnz(~nan_idx);
               end
        end

        %fill in lower triangle
        D = D + D';
end

out_ds = mkgctstruct(D,'rid',ds.cid,...
        'cid',ds.cid,...
        'cdesc',ds.cdesc,...
        'chd',ds.chd,...
        'cdict',ds.cdict);

end