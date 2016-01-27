function DI = compute_diversity_index_two(ds1, ds2)
%Given expression matrices for two cell lines, computes their diversity index.
%
%

pnames = {'bioa_quantile','corr_quantile'};
dflts = {1,.5};
args = parse_args(pnames, dflts, varargin{:});

% Data validation
assert(isequal(ds1.rid,ds2.rid), 'Gene space not the same')

perts = ds1.cid;

DI = 0;
for i = 1:numel(perts)
    
    this_bioa = 
    

    
    
    DI = DI + this_bioa + this_corr;
end

end

