function tsne_ds = mk_ccle_rna_seq_tsne_function(varargin)
%Runs tsne on the ccle rnaseq data set and writes an annotated text file
%which can be passed to tableau.

params = {'metric',...
    'normalize',...
    'perplexity',...
    'num_pcs',...
    'write_text_file',...
    'outdir'};
dflts = {'euclidean',...
    1,...
    10,...
    30,...
    false,...
    '/cmap/projects/cell_line_diversity/analysis/ccle_rna_seq_viz/tsne/'};
args = parse_args(params,dflts,varargin{:});

ccle_pdist = compute_genex_sim('metric',args.metric,...
    'norm',args.normalize,...
    'num_pcs',args.num_pcs);

%%Run tSNE
tsne = tsne_d(ccle_pdist.mat,[],2,args.perplexity);
tsne_ds = mkgctstruct(tsne,'rid',ccle_pdist.rid,'cid',{'TSNE1','TSNE2'});
annot = meta2struct(ccle_pdist,'dim','row');
tsne_ds = annotate_ds(tsne_ds,annot,'dim','row');

scatter(tsne(:,1),tsne(:,2));

%Tableauify
temp1 = table(tsne_ds.mat(:,1),tsne_ds.mat(:,2),tsne_ds.rid,...
    'VariableNames',{'TSNE_1','TSNE_2','cell_name'});
temp2 = cell2table(tsne_ds.rdesc,'VariableNames',tsne_ds.rhd);
tsne_table = [temp1 temp2];

%Write table
if args.write_text_file
    outfile = sprintf('ccle_rnaseq_tsne_%s_norm_%d_perp_%d_num_pcs_%d',...
        args.metric,...
        args.normalize,...
        args.perplexity,...
        args.num_pcs);

    writetable(tsne_table,fullfile(args.outdir,outfile),...
        'FileType','text',...
        'Delimiter','\t');
end
end