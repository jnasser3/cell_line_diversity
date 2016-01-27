function DI = compute_diversity_index(BIOA,CORR,varargin)
%Given an expression matrix, computes the diversity index. 
%
%Input: BIOA - a struct containing a matrix measuring bioactivity.
% BIOA.mat is a matrix of num_perts x num_lines with BIOA_ij being a
% measure of the bioactivity of perturbagen i in cell line j
%
% CORR - a correlation matrix. 
%
% Output: The diversity index.

pnames = {'bioa_quantile'};
dflts = {1};
args = parse_args(pnames, dflts, varargin{:});

perts = BIOA.cid;

DI = 0;
for i = 1:numel(perts)
    this_bioa = quantile(BIOA.mat(i,:),args.bioa_quantile);
    
    corr_idx = ismember(BIOA.cid(i),CORR.cid);
    C = CORR(corr_idx,corr_idx)
    
    DI = DI + this_bioa + this_corr;
end

end

