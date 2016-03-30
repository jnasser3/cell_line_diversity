%%Get the data
ccle_cosine = parse_gctx('/cmap/projects/cell_line_diversity/analysis/baseline_gex_distance/CCLE_RNA_SEQ_DISTANCE_COSINE_n1022x1022.gctx');
blacklist = parse_grp('/cmap/projects/cell_line_diversity/data/cell_lines_blacklist.grp');
phase1_core_lines = parse_grp('/cmap/projects/cell_line_diversity/data/phase1_core_lines.grp');

%% Parameters
p = 20;

%% Run each of the algos
%Exemplars
ex_unconditional = compute_pmedian(ccle_cosine,p,...
    'exclude',blacklist);
ex_conditional = compute_pmedian(ccle_cosine,p,...
    'exclude',blacklist,'known_ex',phase1_core_lines);

% Maximally diverse
div_unconditional = compute_maximally_diverse(ccle_cosine,p,...
    'exclude',blacklist);
div_conditional = compute_maximally_diverse(ccle_cosine,p,...
    'exclude',blacklist,'known',phase1_core_lines);

%% Make some awesome plots
evaluate_total_diversity(ccle_cosine,ex_unconditional,'exclude',blacklist);
namefig('diversity_objective_ex_unconditional');

evaluate_total_diversity(ccle_cosine,div_unconditional,'exclude',blacklist);
namefig('diversity_objective_div');

run_mds(ccle_cosine,'exclude',blacklist,'highlight',ex_unconditional);
namefig('mds_scatter_ex_unconditional');

run_mds(ccle_cosine,'exclude',blacklist,'highlight',div_unconditional);
namefig('mds_scatter_div_unconditional');

run_mds(ccle_cosine,'exclude',blacklist);
namefig('mds_scatter_no_highlight');