lines = parse_grp('/cmap/projects/cell_line_diversity/data/phase2_core_lines.grp');

ds_affy = parse_gctx('/cmap/data/vdb/cline/cline_gene_n1515x12716.gctx');
ds_rnaseq = parse_gctx('/cmap/projects/stacks/STK035_BASE_RNASEQ/STK035_BASE_RNASEQ_L1KFULL_n1023x12450.gctx');

%Lines missing Affy data
lines(~ismember(lines,ds_affy.cid))

%Lines missing RNASeq data
lines(~ismember(lines,ds_rnaseq.cdesc(:,ds_rnaseq.cdict('cell_id'))))
