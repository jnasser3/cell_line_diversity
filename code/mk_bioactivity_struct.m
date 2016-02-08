function bioa = mk_bioactivity_struct(ds, varargin)
%Given an annotated data set, create a perts x cell lines matrix whose
%(i,j) entry is the cc_q75 of pert i in cell line j.
%
%Input: cell_lines: Which cell lines to include. Default all in ds.
%       make_plot: Boolean. 

%Parse input dataset
if ~isstruct(ds) 
    ds = parse_gctx(ds);
end

cell_lines = unique(ds.cdesc(:,ds.cdict('cell_id')));

%Parse inputs
params = {'cell_lines',...
          'make_plot'};
dflts = {cell_lines,...
        false};
args = parse_args(params, dflts, varargin{:});

if ~iscell(args.cell_lines)
    args.cell_lines = parse_grp(args.cell_lines);
end

%Remove signatures that don't have a cc_q75
bad_idx = cell2mat(ds.cdesc(:,ds.cdict('distil_cc_q75'))) == -666;
bad_cids = ds.cid(bad_idx);
ds = ds_slice(ds, 'cid', bad_cids, 'exclude_cid', true);

%Subset ds to the cell lines of interest
cell_idx = find(ismember(ds.cdesc(:,ds.cdict('cell_id')),args.cell_lines));
ds = ds_slice(ds, 'cidx', cell_idx);

%Get perts that appear in each cell line
perts = get_unique_perts(ds, args.cell_lines);
perts_idx = find(ismember(ds.cdesc(:,ds.cdict('pert_id_dose_time')),perts));
ds = ds_slice(ds, 'cidx', perts_idx);

%Get the cc_q75 for each pert x cell lines
mat = zeros(numel(perts),numel(args.cell_lines));
for ii = 1:numel(perts)
    pert_idx = strcmp(ds.cdesc(:,ds.cdict('pert_id_dose_time')),perts(ii));
    for jj = 1:numel(args.cell_lines)
        cell_idx = strcmp(ds.cdesc(:,ds.cdict('cell_id')),args.cell_lines(jj));
        temp = ds.cdesc(pert_idx & cell_idx, ds.cdict('distil_cc_q75'));
        
        mat(ii,jj) = max(cell2mat(temp(1)),0);      
    end
end

bioa = mkgctstruct(mat,'rid',perts,'cid',args.cell_lines);
inames = get_inames_for_combos(ds,perts);
bioa = annotate_ds(bioa, struct('id',perts,'pert_iname',inames'), 'dim', 'row');

if args.make_plot
    plot_bioa(bioa)
end
end

function uperts = get_unique_perts(ds,lines)
%Finds all pert_id_dose_time combos appearing in all cell lines of ds.

uperts = unique(ds.cdesc(:,ds.cdict('pert_id_dose_time')));

for ii = 1:numel(lines)
    this_cell_idx = strcmp(ds.cdesc(:,ds.cdict('cell_id')),lines(ii));
    this_cell_perts = unique(ds.cdesc(this_cell_idx,ds.cdict('pert_id_dose_time')));
    uperts = intersect(uperts,this_cell_perts);
end

end

function inames = get_inames_for_combos(ds,ids)
%For each pert_id_dose_time combo, get the pert_iname
for ii = 1:numel(ids)
    idx = strcmp(ds.cdesc(:, ds.cdict('pert_id_dose_time')),ids(ii));
    temp = unique(ds.cdesc(idx, ds.cdict('pert_iname')));
    if numel(temp) > 1
        error('ERROR: Multiple pert inames for a combo!')
    else
        inames(ii) = temp;
    end
end
end

function plot_bioa(bioa)
imagesc(bioa.mat);
%namefig(fig, 'Replicate_reproducibility')
set(gca, 'XTick', 1:numel(bioa.cid))
set(gca, 'XTickLabel', bioa.cid)

ax = gca;
ax.XTickLabelRotation=45;

set(gca, 'YTick', 1:numel(bioa.rid))
set(gca, 'YTickLabel', bioa.rdesc(:,bioa.rdict('pert_iname')))
set(gca, 'TickLength',[0 0])
colormap(cool)
colorbar 

title_str = sprintf('Replicate Reproducability for %d cell lines and %d perturbagens',...
    numel(bioa.cid),numel(bioa.rid));
title(title_str)

end