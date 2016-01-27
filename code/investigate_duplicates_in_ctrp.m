ds = parse_gctx('/cmap/projects/cell_line_diversity/data/ctrp/ctrp_auc_n907x545.gctx');

dup_cell_id = duplicates(ds.cid);

for ii = 1:numel(dup_cell_id)
	   idx = find(strcmp(ds.cid,dup_cell_id(ii)))
	   ds.cdesc(idx,:)
	   ds.cid(idx)
	   ds.mat(:,idx)

end
