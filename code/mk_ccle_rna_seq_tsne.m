%Make ccle RNA-Seq tsne. Annotate. Tableauify
PERPLEXITY = 15;
INITIAL_DIM = 200;
data_file = '/cmap/projects/cell_line_diversity/data/CCLE_BASELINE_RNASEQ_L1KFULL_RPKM_LOG2_n1022x12450.gctx';

outdir = '/cmap/projects/cell_line_diversity/analysis/ccle_rna_seq_viz/tsne/';
outfile = strcat('ccle_rnaseq_tsne_annotated_perp_',num2str(PERPLEXITY),'_','initdim_',num2str(INITIAL_DIM));

%%Run tSNE
sig_tsne_tool('--ds',data_file,...
    '--create_subdir',0,...
    '--out','/cmap/projects/cell_line_diversity/analysis/ccle_rna_seq_viz/tsne',...
    '--initial_dim',INITIAL_DIM,...
    '--perplexity',PERPLEXITY)

tsne_ds = parse_gctx('/cmap/projects/cell_line_diversity/analysis/ccle_rna_seq_viz/tsne/tsne_n2x1022.gctx');

%Tableauify
temp1 = table(tsne_ds.mat(:,1),tsne_ds.mat(:,2),tsne_ds.rid,...
    'VariableNames',{'TSNE_1','TSNE_2','cell_name'});
temp2 = cell2table(tsne_ds.rdesc,'VariableNames',tsne_ds.rhd);
tsne_table = [temp1 temp2];

%Write table
writetable(tsne_table,fullfile(outdir,outfile),...
    'FileType','text',...
    'Delimiter','\t');