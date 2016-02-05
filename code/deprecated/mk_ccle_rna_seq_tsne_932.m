%Make ccle RNA-Seq tsne. Annotate. Tableauify

% data_file = '/cmap/data/vdb/cline/ccle_rnaseq_rpkm_n932x23686.gctx';
% 
% sig_tsne_tool('--ds',data_file,...
%     '--create_subdir',0,...
%     '--out','/cmap/projects/cell_line_diversity/analysis/ccle_rna_seq_viz/tsne',...
%     '--perplexity',10)

tsne_ds = parse_gctx('/cmap/projects/cell_line_diversity/analysis/ccle_rna_seq_viz/tsne/tsne_n2x932.gctx');

%annotate w/ cline db info
annot_table = parse_tbl('/cmap/data/vdb/cline/clinedb.txt','outfmt','record');
tsne_ds = annotate_ds(tsne_ds,annot_table,...
    'keyfield','cell_id',...
    'dim','row',...
    'skipmissing',true);

%annotate w/ prism info
prism_table = readtable('/cmap/projects/cell_line_diversity/data/CalicoTranche1PrimaryMetaData_01212016.xlsx',...
    'FileType','spreadsheet',...
    'Range','A1:N579');
prism_struct = table2struct(prism_table);
[prism_struct.name] = prism_struct.StrippedName;
prism_struct = rmfield(prism_struct,'StrippedName');
tsne_ds = annotate_ds(tsne_ds,prism_struct,...
    'keyfield','name',...
    'dim','row',...
    'skipmissing',true);

%Get tissue type
tissue_type = {};
for ii = 1:numel(tsne_ds.rid)
    try
        str = tsne_ds.rdesc{ii,tsne_ds.rdict('cell_lineage')};
        comma_loc = strfind(str,',');
        t1 = comma_loc(1);
        str(1:(t1-1));
        tissue_type(ii) = {str(1:(t1-1))};
    catch
        tissue_type(ii) = tsne_ds.rdesc(ii,tsne_ds.rdict('cell_lineage'));
    end
end

%Tableauify
temp1 = table(tsne_ds.mat(:,1),tsne_ds.mat(:,2),tissue_type',...
    'VariableNames',{'TSNE_1','TSNE_2','tissue_type'});
temp2 = cell2table(tsne_ds.rdesc,'VariableNames',tsne_ds.rhd);
tsne_table = [temp1 temp2];

%Write table
writetable(tsne_table,'/cmap/projects/cell_line_diversity/analysis/ccle_rna_seq_viz/tsne/tsne_tableau_932.txt',...
    'FileType','text',...
    'Delimiter','\t');