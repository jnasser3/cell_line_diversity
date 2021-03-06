%% Get the data
%ccle = parse_gctx('/cmap/projects/cell_line_diversity/analysis/baseline_gex_distance/CCLE_RNA_SEQ_DISTANCE_COSINE_n1022x1022.gctx');
ccle_pdist = compute_genex_sim('metric','cosine','norm',1);

blacklist = parse_grp('/cmap/projects/cell_line_diversity/data/spaces/cell_lines_blacklist.grp');
phase1_core_lines = parse_grp('/cmap/projects/cell_line_diversity/data/spaces/phase1_core_lines.grp');
prism_lines = parse_grp('/cmap/projects/cell_line_diversity/data/spaces/prism_lines.grp');
atcc_lines = parse_grp('/cmap/projects/cell_line_diversity/data/spaces/atcc_lines.grp');

%% Make subsetted datasets
% ccle_pdist_prism_only = ds_slice(ccle_pdist,...
%     'cid',prism_lines,...
%     'rid',prism_lines,...
%     'ignore_missing',true);

%% Run each of the algos
%Exemplars - only choose prism lines
ex_unconditional_prism_candidate = run_apcluster(ccle_pdist,...
    'include',prism_lines,...
    'exclude',blacklist,...
    'pref_quantile',.00435);

ex_unconditional_prism_atcc_candidate = run_apcluster(ccle_pdist,...
    'include',intersect(atcc_lines, prism_lines),...
    'exclude',blacklist,...
    'pref_quantile',.0073);

ex_unconditional_no_candidate = run_apcluster(ccle_pdist,...
    'exclude',blacklist,...
    'candidate','',...
    'pref_quantile',.0090);

%% Make some plots
evaluate_total_diversity(ccle_pdist,...
    ex_unconditional_prism_candidate,...
    'objective','repr',...
    'exclude',blacklist);
namefig('exemplar_objective_unconditional_prism_candidate');

evaluate_total_diversity(ccle_pdist,...
    ex_unconditional_prism_atcc_candidate,...
    'objective','repr',...
    'exclude',blacklist);
namefig('exemplar_objective_unconditional_prism_atcc_candidate');

evaluate_total_diversity(ccle_pdist,...
    ex_unconditional_no_candidate,...
    'objective','repr',...
    'exclude',blacklist);
namefig('exemplar_objective_unconditional_no_candidate');

mk_pdist_choice_diagnostic_plot(ccle_pdist,...
    ex_unconditional_prism_candidate,...
    'exclude',blacklist)
namefig('exemplar_unconditional_prism_candidate_diagnostic');

mk_pdist_choice_diagnostic_plot(ccle_pdist,...
    ex_unconditional_prism_atcc_candidate,...
    'exclude',blacklist)
namefig('exemplar_unconditional_prism_atcc_candidate_diagnostic');

mk_pdist_choice_diagnostic_plot(ccle_pdist,...
    ex_unconditional_no_candidate,...
    'exclude',blacklist)
namefig('exemplar_unconditional_no_candidate_diagnostic');

%% mk tables