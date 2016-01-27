function out_ds = pairwise_auc_distance(ds,varargin);
% Computes the pairwise distance matrix from auc pharmacological profiling data
% 
% Inputs:
%        ds - a data structure where the rows of ds.mat are the perts and the columns are cell lines
%        
% Outputs:
%        out_ds - A struct. out_ds.mat will have the pairwise distances (note only the upper triangle i%                           s computed)

[num_cp,num_cell] = size(ds.mat);
D = zeros(num_cell);

for ii = 1:num_cell
	   for jj = ii:num_cell
	  	      nan_idx = isnan(ds.mat(:,ii)) | isnan(ds.mat(:,jj));
                      D(ii,jj) = norm(ds.mat(~nan_idx,ii) - ds.mat(~nan_idx,jj),2) / nnz(~nan_idx); 
           end
end

	   out_ds = mkgctstruct(D,'rid',ds.cid,...
				'cid',ds.cid,...
				'cdesc',ds.cdesc,...
				'chd',ds.chd,...
				'cdict',ds.cdict);
