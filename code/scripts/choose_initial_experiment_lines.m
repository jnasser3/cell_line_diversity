%% Get the data
%ccle = parse_gctx('/cmap/projects/cell_line_diversity/analysis/baseline_gex_distance/CCLE_RNA_SEQ_DISTANCE_COSINE_n1022x1022.gctx');
ccle_pdist = compute_genex_sim('metric','euclidean',...
    'norm',1,...
    'num_pcs',30);

blacklist = parse_grp('/cmap/projects/cell_line_diversity/data/spaces/cell_lines_blacklist.grp');
phase1_core_lines = parse_grp('/cmap/projects/cell_line_diversity/data/spaces/phase1_core_lines.grp');
prism_lines = parse_grp('/cmap/projects/cell_line_diversity/data/spaces/prism_lines.grp');
atcc_lines = parse_grp('/cmap/projects/cell_line_diversity/data/spaces/atcc_lines.grp');
achilles_lines = parse_grp('/cmap/projects/cell_line_diversity/data/spaces/achilles_lines.grp');
cpc006_lines = parse_grp('/cmap/projects/cell_line_diversity/data/spaces/cpc006_lines.grp');

%% Run affinity propogation
ex_unconditional_prism_only_atcc_preference = run_apcluster(ccle_pdist,...
    'include',union(prism_lines,phase1_core_lines),...
    'candidate',intersect(atcc_lines,achilles_lines),...
    'known',phase1_core_lines,...
    'exclude',blacklist,...
    'pref_quantile',.0125);%cosine for 20: .00435

%% Find some close ones
cell = 'TOV112D';
good_lines = prism_lines;
good_idx = ismember(ccle_pdist.cid,good_lines);
cell_idx = find(strcmp(ccle_pdist.cid,cell));
T = array2table(ccle_pdist.mat(good_idx,cell_idx),'RowNames',ccle_pdist.rid(good_idx),'VariableNames',{'Distance_from_TOV112D'});
T_sort = sortrows(T);
close_ones = T_sort.Properties.RowNames(2:3);
close1_idx = strcmp(ccle_pdist.rid,close_ones(1));
close2_idx = strcmp(ccle_pdist.rid,close_ones(2));
ccle_pdist.mat(close1_idx,close2_idx)

%% Find some far ones
good_lines = intersect(atcc_lines,prism_lines);
good_idx = ismember(ccle_pdist.cid,good_lines);
good_ds = ds_slice(ccle_pdist,'rid',good_lines,'cid',good_lines,'ignore_missing',true);
chosen_prism = intersect(ex_unconditional_prism_only_atcc_preference,prism_lines);
chosen_prism_idx = ismember(good_ds.cid,chosen_prism);
T = array2table([sum(ccle_pdist.mat(good_idx,chosen_prism_idx),2), ccle_pdist.mat(good_idx,chosen_prism_idx)],...
    'RowNames',ccle_pdist.rid(good_idx),...
    'VariableNames',['Sum'; chosen_prism]);
T_sort = sortrows(T,-1);
far_ones = T_sort.Properties.RowNames(1);
    
%% Make some plots
evaluate_total_diversity(ccle_pdist,...
    ex_unconditional_prism_only_atcc_preference,...
    'include',union(prism_lines,phase1_core_lines),...
    'objective','repr',...
    'exclude',blacklist);
namefig('exemplar_objective_unconditional_prism_candidate');

evaluate_total_diversity(ccle_pdist,...
    union(ex_unconditional_prism_only_atcc_preference,[close_ones; far_ones]),...
    'include',union(prism_lines,phase1_core_lines),...
    'objective','repr',...
    'exclude',blacklist);
namefig('exemplar_objective_unconditional_prism_candidate');

mk_pdist_choice_diagnostic_plot(ccle_pdist,...
    union(ex_unconditional_prism_only_atcc_preference,[close_ones; far_ones]),...
    'include',union(prism_lines,phase1_core_lines),...
    'exclude',blacklist)
namefig('exemplar_unconditional_prism_atcc_candidate_diagnostic');

%% Mk tables
chosen = setdiff(ex_unconditional_prism_only_atcc_preference,phase1_core_lines);
chosen_plus = union(chosen,[close_ones; far_ones]);
idx = ismember(ccle_pdist.cid,chosen_plus);

lineage = ccle_pdist.cdesc(idx,ccle_pdist.cdict('tissue_type'));
source = ccle_pdist.cdesc(idx,ccle_pdist.cdict('cell_source'));
histology = ccle_pdist.cdesc(idx,ccle_pdist.cdict('cell_histology'));
is_prism = ismember(chosen_plus,prism_lines);
is_achilles = ismember(chosen_plus,achilles_lines);
is_cpc006 = ismember(chosen_plus,cpc006_lines);

% Mutations
%mut = mongo_info('cell_info', chosen, 'fields',{'cell_id','mutations'});

% cell_table = cell2table([chosen_plus';lineage';histology';source';num2cell(is_prism');num2cell(is_achilles');num2cell(is_cpc006')]',...
%     'VariableNames',{'cell_id','lineage','histology','source','is_prism','is_achilles','is_cpc006'});
% writetable(cell_table,'/cmap/projects/cell_line_diversity/analysis/initial_experiment_candidates.txt',...
%     'Delimiter','|')
