%Make ccle RNA-Seq tsne. Annotate. Tableauify
PERPLEXITY = 30;
INITIAL_DIM = 50;
THETA = .5;
data_file = '/cmap/projects/cell_line_diversity/data/cmap_training_data/CPC006_n22098x978.gctx';

outdir = '/cmap/projects/cell_line_diversity/analysis/cpc006/tsne/';
outfile = strcat('cpc006_tsne_annotated_perp_',num2str(PERPLEXITY),'_','initdim_',num2str(INITIAL_DIM),'_theta_',num2str(THETA),'.txt');

%%Run tSNE
sig_tsne_tool('--ds',data_file,...
    '--create_subdir',0,...
    '--out','/cmap/projects/cell_line_diversity/analysis/cpc006/tsne',...
    '--initial_dim',INITIAL_DIM,...
    '--perplexity',PERPLEXITY,...
    '--algorithm','barnes-hut',...
    '--theta',THETA)

tsne_ds = parse_gctx('/cmap/projects/cell_line_diversity/analysis/cpc006/tsne/tsne_n2x22098.gctx');

%Tableauify
temp1 = table(tsne_ds.mat(:,1),tsne_ds.mat(:,2),tsne_ds.rid,...
    'VariableNames',{'TSNE_1','TSNE_2','sig_id'});
temp2 = cell2table(tsne_ds.rdesc,'VariableNames',tsne_ds.rhd);
tsne_table = [temp1 temp2];

%Write table
writetable(tsne_table,fullfile(outdir,outfile),...
    'FileType','text',...
    'Delimiter','\t');