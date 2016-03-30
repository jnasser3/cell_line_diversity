%Just get some basic stats from the CPC006 data

cpc006 = parse_gctx('/cmap/projects/cell_line_diversity/data/cmap_training_data/CPC006_n22098x978.gctx');

%% Number of replicates
cp_idx = strcmp(cpc006.cdesc(:,cpc006.cdict('pert_type')),'trt_cp');
nrep = cell2mat(cpc006.cdesc(cp_idx,cpc006.cdict('distil_nsample')));
histogram(nrep)
xlabel('Number of replicates')
ylabel('Count')
title_str = sprintf('CPC006 Compounds - Distribution of Number of Replicates\n %d Compounds', sum(cp_idx));
title(title_str)
grid on
