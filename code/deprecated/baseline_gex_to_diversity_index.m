
core_lines = parse_grp('/cmap/projects/cell_line_diversity/data/phase1_core_lines.grp');
pdi_ts = parse_gctx('/cmap/projects/cell_line_diversity/analysis/pairwise_diversity_index/core_ts/pdi_n9x9.gctx');
pdi_cpc = parse_gctx('/cmap/projects/cell_line_diversity/analysis/pairwise_diversity_index/cpc006/pdi_n51x51.gctx');
sim_old = compute_genex_sim('ccleds','rnaseq','norm',1);
sim_2015 = compute_genex_sim('ccleds','rnaseq2015','norm',1);

%% Baseline to old rnaseq - core
% feature2diversity_correspondence(sim_old,pdi_ts,'cells',core_lines)
% feature2diversity_correspondence(sim_old,pdi_cpc,'cells',core_lines)

%% Baseline to new rnaseq - core
feature2diversity_correspondence(sim_2015,pdi_ts,'cells',core_lines,'feature_cell_name_identifier','cell_id')
feature2diversity_correspondence(sim_2015,pdi_cpc,'cells',core_lines,'feature_cell_name_identifier','cell_id')

%%
feature2diversity_correspondence(sim_2015,pdi_ts,'feature_cell_name_identifier','cell_id')
feature2diversity_correspondence(sim_2015,pdi_cpc,'feature_cell_name_identifier','cell_id')