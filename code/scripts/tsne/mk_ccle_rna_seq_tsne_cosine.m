%Make ccle RNA-Seq tsne. Annotate. Tableauify
PERPLEXITY = 15;
ccle_cosine = parse_gctx('/cmap/projects/cell_line_diversity/analysis/baseline_gex_distance/CCLE_RNA_SEQ_DISTANCE_COSINE_n1022x1022.gctx');

% blacklist = parse_grp('/cmap/projects/cell_line_diversity/data/cell_lines_blacklist.grp');
% ccle_cosine = ds_slice(ccle_cosine,'cid',setdiff(ccle_cosine.cid,blacklist),'rid',setdiff(ccle_cosine.rid,blacklist));

outdir = '/cmap/projects/cell_line_diversity/analysis/ccle_rna_seq_viz/tsne/';
outfile = strcat('ccle_rnaseq_tsne_cosine_annotated_perp_',num2str(PERPLEXITY));

%%Run tSNE
tsne = tsne_d(ccle_cosine.mat,[],2,PERPLEXITY);
tsne_ds = mkgctstruct(tsne,'rid',ccle_cosine.rid,'cid',{'TSNE1','TSNE2'});
annot = meta2struct(ccle_cosine,'dim','row');
tsne_ds = annotate_ds(tsne_ds,annot,'dim','row');

%Tableauify
temp1 = table(tsne_ds.mat(:,1),tsne_ds.mat(:,2),tsne_ds.rid,...
    'VariableNames',{'TSNE_1','TSNE_2','cell_name'});
temp2 = cell2table(tsne_ds.rdesc,'VariableNames',tsne_ds.rhd);
tsne_table = [temp1 temp2];

%Write table
writetable(tsne_table,fullfile(outdir,outfile),...
    'FileType','text',...
    'Delimiter','\t');