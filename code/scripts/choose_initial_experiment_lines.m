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

KNN = 3;

%% Find some close ones
% cell = 'TOV112D';
% good_lines = prism_lines;
% good_idx = ismember(ccle_pdist.cid,good_lines);
% cell_idx = find(strcmp(ccle_pdist.cid,cell));
% T = array2table(ccle_pdist.mat(good_idx,cell_idx),'RowNames',ccle_pdist.rid(good_idx),'VariableNames',{'Distance_from_TOV112D'});
% T_sort = sortrows(T);
% close_ones = T_sort.Properties.RowNames(2:3);
% close1_idx = strcmp(ccle_pdist.rid,close_ones(1));
% close2_idx = strcmp(ccle_pdist.rid,close_ones(2));
% ccle_pdist.mat(cell_idx,close1_idx)
% ccle_pdist.mat(cell_idx,close2_idx)
% ccle_pdist.mat(close1_idx,close2_idx)

%% Clique finding
% good_lines = setdiff(prism_lines,blacklist);
% good_ds = ds_slice(ccle_pdist,'cid',good_lines,'rid',good_lines,'ignore_missing',true);
% ccle_sim = good_ds;
% ccle_sim.mat = -1*ccle_sim.mat;
% cutoff = prctile(tri2vec(ccle_sim.mat),99);
% max_clique = clique_heuristic(ccle_sim,'cutoff',cutoff,'num_cliques',20);

%% Run affinity propogation
preference_lines = intersect(prism_lines,atcc_lines);
preference_lines = intersect(preference_lines,achilles_lines);
ex_unconditional_prism_achilles_atcc_preference = run_apcluster(ccle_pdist,...
    'preference',preference_lines,...
    'exclude',blacklist,...
    'num_exemplar',9); %'pref_quantile',.03

ex_unconditional_no_preference = run_apcluster(ccle_pdist,...
    'exclude',blacklist,...
    'num_exemplar',9);

% ex_unconditional_prism_achilles_atcc_preference = ...
% {        'HS940T'
%     'HT29'
%     'MCF7'
%     'NCIH1648'
%     'NCIH2009'
%     'NCIH2081'
%     'SKMEL5'
%     'T98G'
%     'TOV112D'};


%% Find close lines
close_map = containers.Map;
cands = ex_unconditional_prism_achilles_atcc_preference;
for jj = 1:numel(cands)
    nn = find_nearest_neighbors(ccle_pdist,cands(jj),KNN,...
        'exclude',blacklist);
    close_map(cands{jj}) = nn;
end

cl1 = close_map('SKMEL5');
cl2 = close_map('TOV112D');
close_ones = vertcat(cl1,cl2);

%% Find some far ones
%good_lines = intersect(atcc_lines,prism_lines);
%good_idx = ismember(ccle_pdist.cid,good_lines);
%good_ds = ds_slice(ccle_pdist,'rid',good_lines,'cid',good_lines,'ignore_missing',true);
chosen_prism = intersect(ex_unconditional_prism_achilles_atcc_preference,prism_lines);
chosen_prism_idx = ismember(ccle_pdist.cid,chosen_prism);
T = array2table([sum(ccle_pdist.mat(good_idx,chosen_prism_idx),2), ccle_pdist.mat(good_idx,chosen_prism_idx)],...
    'RowNames',ccle_pdist.rid(good_idx),...
    'VariableNames',['Sum'; chosen_prism]);
T_sort = sortrows(T,-1);
far_ones = T_sort.Properties.RowNames(1:5);

far_ones = [];
    
%% Make some plots

evaluate_total_diversity(ccle_pdist,...
    {ex_unconditional_prism_achilles_atcc_preference,ex_unconditional_no_preference,phase1_core_lines},...
    'lines_labels',{'exemplar_preference','exemplar_no_preference','phase1_core'},...
    'objective','repr',...
    'exclude',blacklist);

% evaluate_total_diversity(ccle_pdist,...
%     ex_unconditional_prism_achilles_atcc_preference,...
%     'objective','repr',...
%     'exclude',blacklist);
% namefig('exemplar_objective_unconditional_prism_achilles_atcc_candidate');
% 
% evaluate_total_diversity(ccle_pdist,...
%     phase1_core_lines,...
%     'objective','repr',...
%     'exclude',blacklist);
% namefig('exemplar_objective_phase1_core');
% 
% evaluate_total_diversity(ccle_pdist,...
%     ex_unconditional_no_preference,...
%     'objective','repr',...
%     'exclude',blacklist);
% namefig('exemplar_objective_unconditional_no_preference');

% evaluate_total_diversity(ccle_pdist,...
%     union(ex_unconditional_prism_achilles_atcc_preference,[close_ones; far_ones]),...
%     'objective','repr',...
%     'exclude',blacklist);
% namefig('exemplar_objective_unconditional_prism_candidate');

mk_pdist_choice_diagnostic_plot(ccle_pdist,...
    union(ex_unconditional_prism_achilles_atcc_preference,[close_ones; far_ones]),...
    'exclude',blacklist)
namefig('exemplar_unconditional_prism_atcc_candidate_diagnostic');

%% Mk tables
chosen = setdiff(ex_unconditional_prism_achilles_atcc_preference,[]);
chosen_plus = union(chosen,[close_ones; far_ones]);
idx = ismember(ccle_pdist.cid,chosen_plus);

lineage = ccle_pdist.cdesc(idx,ccle_pdist.cdict('tissue_type'));
source = ccle_pdist.cdesc(idx,ccle_pdist.cdict('cell_source'));
histology = ccle_pdist.cdesc(idx,ccle_pdist.cdict('cell_histology'));
is_prism = ismember(chosen_plus,prism_lines);
is_achilles = ismember(chosen_plus,achilles_lines);
is_cpc006 = ismember(chosen_plus,cpc006_lines);

% Mutations
% mut = mongo_info('cell_info', chosen, 'fields',{'cell_id','mutations'});

cell_table = cell2table([chosen_plus';lineage';histology';source';num2cell(is_prism');num2cell(is_achilles');num2cell(is_cpc006')]',...
    'VariableNames',{'cell_id','lineage','histology','source','is_prism','is_achilles','is_cpc006'});
writetable(cell_table,'/cmap/projects/cell_line_diversity/analysis/initial_experiment_candidates.txt',...
    'Delimiter','|')

%% Mk Plots
%chosen_plus_phase1 = union(chosen_plus,phase1_core_lines);
run_mds(ccle_pdist,'exclude',blacklist,'highlight',chosen,'highlight_label','Exemplar_preference')
run_mds(ccle_pdist,'exclude',blacklist,'highlight',ex_unconditional_no_preference,'highlight_label','Exemplar_no_preference')
run_mds(ccle_pdist,'exclude',blacklist,'highlight',phase1_core_lines,'highlight_label','Phase1_core')
%run_mds(ccle_pdist,'exclude',blacklist,'highlight',chosen_plus)

%% Get cluster ID's
cl_id = assign_cluster_ids(ccle_pdist,union(chosen,phase1_core_lines),'exclude',blacklist,'write_table',true)

