function ds_out = mk_conditional_pdist(ds_in,known_exemplar)
%ds_out = mk_conditional_pdist(ds_in,known_exemplar)
%
%Transforms a pairwise distance matrix to account for known examplars.
%Solving the conditional p-medians problem on the original matrix is
%equivalent to solving the unconditional problem on the transfored matrix.
%
%Based on the reference
%O. Berman, Z. Drezner, A new formulation for the conditional p-median and pcenter
%problems, Operations Research Letters 36 (2008) 481?483.
%
%Input: 
%       ds_in: a square struct of pairwise distances
%       known_ex: A list of known exemplars
%
%Output: 
%       ds_out: A square struct with a modified distance matrix. 
%               The known exemplars are not included. Note that
%               that ds_out.mat need not be symmetric

known_idx = ismember(ds_in.rid,known_exemplar);
n = numel(ds_in.rid);

ds_out = ds_in;
D = ds_in.mat;

%transform the distance matrix
for ii = 1:n
    for jj = 1:n
        D(ii,jj) = min(D(ii,jj),min(D(ii,known_idx)));
    end
end

%Return a struct excluding known lines
non_exemplar = setdiff(ds_in.rid,known_exemplar);
ds_out.mat = D;
ds_out = ds_slice(ds_out,'rid',non_exemplar,'cid',non_exemplar);

end

