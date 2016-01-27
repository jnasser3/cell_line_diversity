
lm_path = mortar.common.Spaces.probe_space('lm');
lm_genes = parse_grp(lm_path{1});

%Get some replicates
meta = parse_tbl('/cmap/projects/M1/TA.OE/data/build_LUAD_A549/instinfo.txt',...
    'outfmt','record');

