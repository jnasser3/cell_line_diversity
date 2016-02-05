%Make ccle RNA-Seq tsne. Annotate. Tableauify

data_file = '/cmap/projects/cell_line_diversity/data/CCLE_BASELINE_RNASEQ_L1KFULL_RPKM_LOG2_n1022x12450.gctx';

sig_tsne_tool('--ds',good_ds,...
    '--create_subdir',0,...
    '--out','/cmap/projects/cell_line_diversity/analysis/ccle_rna_seq_viz/tsne',...
    '--initial_dim',200,...
    '--perplexity',15)

tsne_ds = parse_gctx('/cmap/projects/cell_line_diversity/analysis/ccle_rna_seq_viz/tsne/tsne_n2x1022.gctx');

%annotate
%extract the tissue type from the rid
tissue_type = {};
for ii = 1:numel(tsne_ds.rid)
    try
        str = tsne_ds.rid{ii};
        underbar_loc = strfind(str,'_');
        t1 = underbar_loc(2);
        dash_loc = strfind(str,'-');
        t2 = dash_loc(end);
        tissue_type(ii) = {str((t1+1):(t2-1))};
    catch
        tissue_type(ii) = {'-666'};
    end
end

%Tableauify
temp1 = table(tsne_ds.mat(:,1),tsne_ds.mat(:,2),tsne_ds.rid,tsne_ds.rdesc(:,tsne_ds.rdict('cell_id')),tissue_type',...
    'VariableNames',{'TSNE_1','TSNE_2','cell_name','cell_id','tissue_type'});
temp2 = cell2table(tsne_ds.rdesc,'VariableNames',tsne_ds.rhd);
tsne_table = join(temp1,temp2,'Keys','cell_id');

%Write table
writetable(tsne_table,'/cmap/projects/cell_line_diversity/analysis/ccle_rna_seq_viz/tsne/tsne_tableau_1022.txt',...
    'FileType','text',...
    'Delimiter','\t');