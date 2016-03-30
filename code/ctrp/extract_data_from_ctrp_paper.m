% Extract data from the CTRP paper into a ds struct

% Get meta info
cp_meta = readtable('/cmap/projects/cell_line_diversity/data/ctrp_paper/ctrp_v2_data_paper.xlsx',...
    'FileType','spreadsheet',...
    'sheet','S1');

cell_meta = readtable('/cmap/projects/cell_line_diversity/data/ctrp_paper/ctrp_v2_data_paper.xlsx',...
    'FileType','spreadsheet',...
    'sheet','S2');

% Get AUC data
auc_data = readtable('/cmap/projects/cell_line_diversity/data/ctrp_paper/ctrp_v2_data_paper.xlsx',...
    'FileType','spreadsheet',...
    'sheet','S3');

% Load in the data
data = NaN(size(cp_meta,1),size(cell_meta,1));
for ii = 1:size(auc_data,1)
    data(auc_data{ii,1},auc_data{ii,2}) = auc_data{ii,3};  
end

%make ds
ds = mkgctstruct(data,'rid',strtrim(cellstr(int2str(cp_meta{:,1}))),'cid',strtrim(cellstr(int2str(cell_meta{:,1}))));

%annotate ds - warning: hacky
cp_meta_struct = table2struct(cp_meta);
cp_meta_annot = parse_tbl(cp_meta_struct,'outfmt','record');
temp = cellfun(@(s) num2str(s), {cp_meta_annot.('index_cpd')},'UniformOutput',false);
[cp_meta_annot.('index_cpd')] = temp{:};
ds = annotate_ds(ds, cp_meta_annot, 'dim', 'row', 'keyfield', 'index_cpd');

cell_meta_struct = table2struct(cell_meta);
cell_meta_annot = parse_tbl(cell_meta_struct,'outfmt','record');
temp = cellfun(@(s) num2str(s), {cell_meta_annot.('index_ccl')},'UniformOutput',false);
[cell_meta_annot.('index_ccl')] = temp{:};
ds = annotate_ds(ds, cell_meta_annot, 'dim', 'column', 'keyfield', 'index_ccl');

mkgctx('/cmap/projects/cell_line_diversity/data/ctrp_paper/ctrp_auc_paper.gctx',ds)