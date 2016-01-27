% In the cell diversity section of M1 we describe that certain compounds are 
% bioactive in all cell lines whereas others are context specfic. The examples
% given are HDAC inhibitors and BRAF KD. Lets take a look at them

sig = sig_info('{"pert_iname":"BRAF","pert_type":"trt_sh.cgs","cell_id":{"$in":["MCF7","PC3","HT29","A375","A549","HEPG2","HA1E","VCAP","HCC515"]}}');

lm = mortar.common.Spaces.probe_space('lm');
lm_genes = parse_grp(lm{1});
ds = parse_gctx('/cmap/data/build/a2y13q1/modzs.gctx',...
		'rid',lm_genes,...
		'cid',{sig.sig_id})
		
ds = annotate_ds(ds,sig,'dim','column','keyfield','sig_id')
