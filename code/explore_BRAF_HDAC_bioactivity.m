% In the cell diversity section of M1 we describe that certain compounds are 
% bioactive in all cell lines whereas others are context specfic. The examples
% given are HDAC inhibitors and BRAF KD. Lets take a look at them

vem_sig = sig_info('{"pert_iname":"vemurafenib"}');

lm = mortar.common.Spaces.probe_space('lm');
lm_genes = parse_grp(lm{1});
vem_ds = parse_gctx('/cmap/data/build/a2y13q1/modzs.gctx',...
		'rid',lm_genes,...
		'cid',{vem_sig.sig_id})
		

