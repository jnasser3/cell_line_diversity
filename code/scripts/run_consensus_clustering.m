%% Load ccle data
ccle = parse_gctx('/cmap/projects/cell_line_diversity/data/rnaseq/CCLE_BASELINE_RNASEQ_PROTEIN_CODING_RPKM_LOG2_LABELED_n1022x18772.gctx');

%Remove hematopoietic lines
% blacklist = parse_grp('/cmap/projects/cell_line_diversity/data/spaces/cell_lines_blacklist.grp');
% ccle = ds_slice(ccle,'cid',blacklist,'exclude_cid',true);

% cpc006_lines = parse_grp('/cmap/projects/cell_line_diversity/data/spaces/cpc006_lines.grp')

%ccle = ds_slice(ccle,'cidx',1:10)
% ccle = ds_subset(ccle,'column','cell_id',cpc006_lines)

ccle.cid

%Sort ccle based on cids
[~,sort_idx] = sort(ccle.cid);
ccle = ds_slice(ccle,'cidx',sort_idx);

ccle.cid

ccle_all_genes = ccle.rid;

%% Setup parameters
num_ex = 50;

%% Run consensus clustering by varying genes
%Make some gene sets
lm = parse_grp('/cmap/data/vdb/spaces/lm_epsilon_symbols_n978.grp');
bing = parse_grp('/cmap/data/vdb/spaces/bing_symbols_n10638.grp');
temp1 = randperm(length(ccle_all_genes),5000); rand1 = ccle_all_genes(temp1);
temp2 = randperm(length(ccle_all_genes),5000); rand2 = ccle_all_genes(temp2);
temp3 = randperm(length(ccle_all_genes),5000); rand3 = ccle_all_genes(temp3);
temp4 = randperm(length(ccle_all_genes),5000); rand4 = ccle_all_genes(temp4);
temp5 = randperm(length(ccle_all_genes),5000); rand5 = ccle_all_genes(temp5);
temp6 = randperm(length(ccle_all_genes),5000); rand6 = ccle_all_genes(temp6);

gene_sets = {lm, bing, ccle_all_genes, rand1, rand2, rand3};

cluster_mat = zeros(numel(ccle.cid),numel(gene_sets));
for ii = 1:size(cluster_mat,2);
    
    %Compute pdist
    this_ccle = ds_subset(ccle,'row','desc',gene_sets{ii});
    this_pdist = compute_genex_sim2(this_ccle,...
        'metric','euclidean',...
        'num_pcs',30);
    this_psim = this_pdist;
    this_psim.mat = -1*this_pdist.mat;
    
    %Run clustering
    [this_exemplar,this_membership] = run_apcluster(this_psim,...   
        'num_exemplar',num_ex,...
        'num_exemplar_tol',0,...
        'high_quantile_start',.01);
    
    %convert membership to integer based clustering
    [~,this_clustering] = ismember(this_membership,this_exemplar);
    
    cluster_mat(:,ii) = this_clustering;
end

% Compute similarity matrix
sim_mat = evaluate_consensus_clustering(cluster_mat,'Gene Space');
namefig('Gene Space');

%% Vary number of pcs
num_pc_sets = [20,30,50];
cluster_mat = zeros(numel(ccle.cid),numel(num_pc_sets));
for ii = 1:size(cluster_mat,2);
    
    %Compute pdist
    this_pdist = compute_genex_sim2(ccle,...
        'metric','mahal',...
        'norm',0,...
        'num_pcs',num_pc_sets(ii));
    this_psim = this_pdist;
    this_psim.mat = -1*this_pdist.mat;
    
    %Run clustering
    [this_exemplar,this_membership] = run_apcluster(this_psim,...   
        'num_exemplar',num_ex,...
        'num_exemplar_tol',0,...
        'high_quantile_start',.01);
    
    %convert membership to integer based clustering
    [~,this_clustering] = ismember(this_membership,this_exemplar);
    
    cluster_mat(:,ii) = this_clustering;
end

% Compute similarity matrix
sim_mat = evaluate_consensus_clustering(cluster_mat,'Number of PCs');
namefig('Num PCs')

%% Vary distance metric
metric_choices = {'mahal','euclidean','cosine','spearman'};
num_pcs = [30,30,0,0];
norms = [0,1,1,0];
cluster_mat = zeros(numel(ccle.cid),numel(num_pcs));
for ii = 1:size(cluster_mat,2);
    
    %Compute pdist
    this_pdist = compute_genex_sim2(ccle,...
        'metric',metric_choices{ii},...
        'num_pcs',num_pcs(ii),...
        'norm',norms(ii));
    this_psim = this_pdist;
    this_psim.mat = -1*this_pdist.mat;
    
    %Run clustering
    [this_exemplar,this_membership] = run_apcluster(this_psim,...   
        'num_exemplar',num_ex,...
        'num_exemplar_tol',0,...
        'high_quantile_start',.01);
    
    %convert membership to integer based clustering
    [~,this_clustering] = ismember(this_membership,this_exemplar);
    
    cluster_mat(:,ii) = this_clustering;
end

% Compute similarity matrix
sim_mat = evaluate_consensus_clustering(cluster_mat,'Distance metrics');
namefig('Distance Metrics')

%% Subsample
subsample_rate = .8;
num_subsample = 5;
cluster_mat = zeros(numel(ccle.cid),num_subsample);

for ii = 1:size(cluster_mat,2);
    
    %Compute pdist
    this_cidx = randperm(numel(ccle.cid),ceil(numel(ccle.cid)*subsample_rate));
    this_cidx = sort(this_cidx);
    this_excluded_cidx = setdiff(1:numel(ccle.cid),this_cidx);
    this_ccle = ds_slice(ccle,'cidx',this_cidx);
    
    this_pdist = compute_genex_sim2(this_ccle,...
        'metric','mahal',...
        'num_pcs',30,...
        'norm',0);
    this_psim = this_pdist;
    this_psim.mat = -1*this_pdist.mat;
    
    %Run clustering
    [this_exemplar,this_membership,ids] = run_apcluster(this_psim,...   
        'num_exemplar',num_ex,...
        'num_exemplar_tol',0,...
        'high_quantile_start',.05);
    
%     this_membership_padded = union(this_membership,ccle_temp(this_excluded_cidx));
%     this_membership_padded = sort(this_membership_padded)
    
    %convert membership to integer based clustering
    [~,this_clustering_temp] = ismember(this_membership,this_exemplar);
    
    this_clustering(this_cidx) = this_clustering_temp;
    this_clustering(this_excluded_cidx) = 0;
    
    cluster_mat(:,ii) = this_clustering;
end

% Compute similarity matrix
sim_mat = evaluate_consensus_clustering(cluster_mat,'Sample Holdout');
namefig('Sample holdout')