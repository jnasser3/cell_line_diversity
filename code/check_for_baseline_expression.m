lines1 = parse_grp('/cmap/projects/cell_line_diversity/data/phase1_core_lines.grp');
lines2 = parse_grp('/cmap/projects/cell_line_diversity/data/phase2_core_lines.grp');

ds_affy = parse_gctx('/cmap/data/vdb/cline/cline_gene_n1515x12716.gctx');
ds_rnaseq = parse_gctx('/cmap/projects/stacks/STK035_BASE_RNASEQ/STK035_BASE_RNASEQ_L1KFULL_n1023x12450.gctx');
ds_rnaseq2 = parse_gctx('/cmap/data/vdb/cline/ccle_rnaseq_rpkm_n932x23686.gctx');

%Phase 2 Lines missing Affy data
fprintf('Phase 2 lines missing Affy: \n')
lines2(~ismember(lines2,ds_affy.cid))

%Phase 2 Lines missing RNASeq data
fprintf('Phase 2 lines missing RNA-Seq: \n')
lines2(~ismember(lines2,ds_rnaseq.cdesc(:,ds_rnaseq.cdict('cell_id'))))

%Phase 2 Lines missing RNASeq data
fprintf('Phase 2 lines missing RNA-Seq 2: \n')
lines2(~ismember(lines2,ds_rnaseq2.cid))

%Phase 1 Lines missing Affy data
fprintf('Phase 1 lines missing Affy: \n')
lines1(~ismember(lines1,ds_affy.cid))

%Phase 1 Lines missing RNASeq data
fprintf('Phase 1 lines missing RNA-Seq: \n')
lines1(~ismember(lines1,ds_rnaseq.cdesc(:,ds_rnaseq.cdict('cell_id'))))

%Phase 1 Lines missing RNASeq data
fprintf('Phase 1 lines missing RNA-Seq 2: \n')
lines1(~ismember(lines1,ds_rnaseq2.cid))

