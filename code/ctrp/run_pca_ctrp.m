ds = parse_gctx('/cmap/projects/cell_line_diversity/data/ctrp_paper/ctrp_auc_paper_clean_n619x349.gctx');

[COEFF, SCORE, LATENT, ~, EXPLAINED, MU] = pca(ds.mat','algorithm','als','Options',statset('Display','iter'));
