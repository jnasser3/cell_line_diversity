%Make a plot of ap objective loss as a function of number of exemplars

%% Compute pairwise distances
ccle_pdist = compute_genex_sim('metric','euclidean',...
    'norm',1,...
    'num_pcs',30);

ccle_psim = ccle_pdist;
ccle_psim.mat = -1*ccle_pdist.mat;
%%
nperms = 1000;
% MIN = 10;
% MAX = 160;
% GRID = 50;
% num_ex_target = MIN:GRID:MAX;
num_ex_target = [10,25,50,75,100,150,200,300];
loss = zeros(1,length(num_ex_target));
null_loss = zeros(nperms,length(num_ex_target));
num_ex_actual = zeros(1,length(num_ex_target));

for ii = 1:length(num_ex_target)
    ex = run_apcluster(ccle_psim,...
        'exclude',blacklist,...
        'num_exemplar',num_ex_target(ii),...
        'num_exemplar_tol',.1,...
        'high_quantile_start',.02);
    
    [this_loss, this_bkg] = evaluate_total_diversity(ccle_pdist,...
        {ex},...
        'exclude',blacklist,...
        'ds_type','distance',...
        'lines_labels',{'exemplar_test'},...
        'nperms',nperms);
    
    loss(ii) = this_loss;
    null_loss(:,ii) = this_bkg;
    
    num_ex_actual(ii) = numel(ex);
end

%% Plot
figure;
scatter(num_ex_actual,loss,100,'k','*',...
    'DisplayName','Exemplar Loss')
grid on
title('Number of exemplars vs AP loss')
xlabel('Number of exemplars')
ylabel('Loss')
hold on
boxplot(null_loss,...
    'positions',num_ex_actual,...
    'labels',num_ex_actual,...
    'whisker',10)
ylim auto
legend show