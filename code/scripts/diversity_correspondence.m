%% See whether lines that are diverse in ccle are diverse in cpc006
cpc006_lines = parse_grp('/cmap/projects/cell_line_diversity/data/spaces/cpc006_lines.grp');

ccle_pdist = compute_genex_sim('metric','mahal',...
    'norm',0,...
    'num_pcs',30);



%% Subset to cpc006 lines
ccle_pdist_temp = ds_subset(ccle_pdist,'column','cell_id',cpc006_lines);
ccle_pdist_temp = ds_subset(ccle_pdist_temp,'row','cell_id',cpc006_lines);

%% Run affinity propogation on ccle data
K = 5;
ccle_sim_temp = ccle_pdist_temp;
ccle_sim_temp.mat = -1*ccle_pdist_temp.mat;
ccle_exemplar = run_apcluster(ccle_sim_temp,...
    'num_exemplar',K);

%%Plot
run_mds(ccle_pdist_temp,'highlight',ccle_exemplar)
run_mds(ccle_pdist,'highlight',cpc006_lines)

evaluate_total_diversity(ccle_pdist_temp,...
    {ccle_exemplar},...
    'lines_labels',{'ccle_exemplar'},...
    'objective','repr',...
    'ds_type','distance');

% Test whether the CCLE exemplars are good for cpc006
cpc006_pdi = parse_gctx('/cmap/projects/cell_line_diversity/analysis/pairwise_diversity_index/cpc006/pdi_n51x51.gctx');
cpc006_sim = cpc006_pdi;
cpc006_sim.mat = -1*cpc006_pdi.mat;

run_mds(cpc006_pdi,'highlight',ccle_exemplar)

evaluate_total_diversity(cpc006_pdi,...
    {ccle_exemplar},...
    'lines_labels',{'ccle_exemplar'},...
    'objective','repr',...
    'ds_type','distance');