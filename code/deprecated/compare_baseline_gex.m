function compare_baseline_gex(ds1,ds2,varargin)
%Compare similarities between cell lines in different baseline GEX datasets

params = {'cells',...
    'ds1_cell_id',...
    'ds2_cell_id',...
    'ds1_gene_id',...
    'ds2_gene_id'};
dflts = {'',...
    'cid',...
    'cid',...
    'rid',...
    'rid'};
args = parse_args(params, dflts, varargin{:});

%make the cid in each ds equal to the common identifier
if ~strcmp(args.ds1_cell_id,'cid')
    ds1.cid = ds1.cdesc(:,ds1.cdict(args.ds1_cell_id));
end
if ~strcmp(args.ds2_cell_id,'cid')
    ds2.cid = ds2.cdesc(:,ds2.cdict(args.ds2_cell_id));
end

%make the rid in each ds equal to the common identifier
if ~strcmp(args.ds1_gene_id,'rid')
    ds1.rid = ds1.rdesc(:,ds1.rdict(args.ds1_gene_id));
end
if ~strcmp(args.ds2_gene_id,'rid')
    ds2.rid = ds2.rdesc(:,ds2.rdict(args.ds2_gene_id));
end

common_lines = intersect(ds1.cid, ds2.cid);

if isempty(args.cells)
    args.cells = common_lines;
elseif ~iscell(args.cells)
    args.cells = parse_grp(args.cells);
end

%Subset data sets to common cell lines
common_lines = intersect(common_lines, args.cells);
ds1 = ds_slice(ds1, 'cid', common_lines);
ds2 = ds_slice(ds2, 'cid', common_lines);
assert(isequal(ds1.cid,ds2.cid),'cids not in same order!')

%Subset data sets to common genes
common_genes = intersect(ds1.rid, ds2.rid);
ds1 = ds_slice(ds1, 'rid', common_genes);
ds2 = ds_slice(ds2, 'rid', common_genes);
assert(isequal(ds1.rid,ds2.rid),'rids not in same order!')

%compute pairwise correlations and null correlations
corrs = fastcorr(ds1.mat, ds2.mat);
pw_corrs = diag(corrs);
null_corrs = tri2vec(corrs);

histogram(pw_corrs, 'DisplayName', 'Matched Lines')
hold on
histogram(null_corrs, 'DisplayName', 'Unmatched Lines')
legend show

end

