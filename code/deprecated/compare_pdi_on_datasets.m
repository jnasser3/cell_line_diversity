
core_lines = parse_grp('/cmap/projects/cell_line_diversity/data/phase1_core_lines.grp');
ds_ts = parse_gctx('/cmap/projects/cell_line_diversity/data/cmap_training_data/cp_pruned_n1890x978.gctx');
ds_cpc = parse_gctx('/cmap/projects/cell_line_diversity/data/cmap_training_data/CPC006_n22098x978.gctx');

pdi_ts = compute_pairwise_diversity_index(ds_ts, 'cell_lines', core_lines);
pdi_cpc = compute_pairwise_diversity_index(ds_cpc, 'cell_lines', core_lines);

corr(tri2vec(pdi_ts.mat),tri2vec(pdi_cpc.mat))