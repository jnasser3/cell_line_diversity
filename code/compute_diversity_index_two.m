function DI = compute_diversity_index_two(ds1, ds2)
%Given expression matrices for two cell lines, computes their diversity index.

pnames = {'bioa_quantile','corr_quantile'};
dflts = {1,.5};
args = parse_args(pnames, dflts, varargin{:});

% Data validation
assert(isequal(ds1.rid,ds2.rid), 'Gene space not the same')

% Order the perts identically b/w ds1 and ds2

perts = ds1.cid;

DI = 0;
for ii = 1:numel(perts)
    this_bioa = max(ds1.cdesc(ii,ds1.cdict('distil_cc_q75')),ds2.cdesc(ii,ds1.cdict('distil_cc_q75')));
    this_corr = corr(ds1.mat(:,ii),ds2.mat(:,ii);
    DI = DI + this_bioa*(1-this_corr);
end

end

