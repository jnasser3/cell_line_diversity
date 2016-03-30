%Make ccle RNA-Seq tsne. Annotate. Tableauify
PERPLEXITY = 15;
ccle_cosine = parse_gctx('/cmap/projects/cell_line_diversity/analysis/baseline_gex_distance/CCLE_RNA_SEQ_DISTANCE_COSINE_n1022x1022.gctx');

blacklist = parse_grp('/cmap/projects/cell_line_diversity/data/cell_lines_blacklist.grp');
good_lines = setdiff(ccle_cosine.rid,blacklist);
ds = ds_slice(ccle_cosine,'rid',good_lines,'cid',good_lines);
tsne = tsne_d(ds.mat,[],3,PERPLEXITY);

%stupid hack to make plot since gscatter sucks
tissue = ds.cdesc(:,ds.cdict('tissue_type'));
utissue = unique(tissue);
for ii = 1:numel(utissue)
    idx = strcmp(tissue,utissue(ii));
    plot3(tsne(idx,1),tsne(idx,2),tsne(idx,3),'.','DisplayName',utissue{ii})
    hold on
end

%stupid hack to make colors since gscatter sucks
% tissue = ds.cdesc(:,ds.cdict('tissue_type'));
% utissue = unique(tissue);
% C = zeros(1,numel(utissue));
% for ii = 1:numel(utissue);
%     idx = strcmp(tissue,utissue(ii));
%     C(idx) = ii;
% end
% figure;
% scatter3(tsne(:,1),tsne(:,2),tsne(:,3),[],C)

