%Create a .grp of blacklisted cell lines

ds = parse_gctx('/cmap/projects/cell_line_diversity/data/CCLE_BASELINE_RNASEQ_L1KFULL_RPKM_LOG2_n1022x12450.gctx');


%Add HAEMATOPOETIC lines to blacklist
idx = strcmp(ds.cdesc(:,ds.cdict('tissue_type')),'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE');
blacklist = ds.cid(idx);

mkgrp('/cmap/projects/cell_line_diversity/data/cell_lines_blacklist_pca.grp',blacklist)