function evaluate_total_diversity(ds, lines, varargin )
% evaluate_total_diversity(ds, lines, varargin )
%
% Start with a universe of objects and pairwise distances between them.
% Given a subset of these objects, computes the sum of their  mutual 
% pairwise distances and compares this to random choices.
%
% Input
%       ds: A square struct of pairwise distances
%       lines: The set of objects to evaluate
%       objective: The objective function to evaluate. Options are divr
%           (sum of distances among chosen lines) or repr (sum of minimum
%           distances from all lines to nearest chosen line. repr is optimized
%           by affinity propogation. Default is 'repr'
%       exclude: A set of id's to remove from the analysis
%       fixed_lines: A subset of objects to include in the calculation but
%                   not to evaluate
%       nperms: Number of permutations to compare against. Default 1000
%

params = {'fixed_lines',...
    'objective',...
    'exclude',...
    'nperms'};
dflts = {'',...
    'repr',...
    '',...
    1000};
args = parse_args(params,dflts,varargin{:});

%Subset ds if needed
exclude = get_array_input(args.exclude,'');
if ~isempty(exclude)
    objects = setdiff(ds.rid,exclude);
    ds = ds_slice(ds,'cid',objects,'rid',objects);
end

%get cell line indices
fixed_lines = get_array_input(args.fixed_lines,'');
fixed_lines_idx = find(ismember(ds.rid,fixed_lines));
non_fixed_idx = find(~ismember(ds.rid,fixed_lines));
lines = get_array_input(lines,'');
lines_idx = find(ismember(ds.rid,lines));

%compute the objective
switch lower(args.objective)
    case 'divr'
        test_statistic = compute_total_diversity(ds,[lines_idx; fixed_lines_idx]);
    case 'repr'
        test_statistic = evaluate_pmedian_loss(ds, [lines; fixed_lines]);
end

%compute a background
n = numel(setdiff(ds.rid,fixed_lines));
k = numel(lines);
bkg = zeros(1,args.nperms);
for ii = 1:args.nperms
    temp = randperm(n,k);
    this_lines_idx = non_fixed_idx(temp);
    this_lines = ds.rid(this_lines_idx);
    
    switch lower(args.objective)
        case 'divr'
            bkg(ii) = compute_total_diversity(ds,[this_lines_idx; fixed_lines_idx]);  
        case 'repr'
            bkg(ii) = evaluate_pmedian_loss(ds, [this_lines; fixed_lines]);
    end
end

%plot
figure;
histogram(bkg,'DisplayName','Background');
hold on
g = gca;
plot([test_statistic test_statistic], g.YLim, 'DisplayName', 'Test Statistic')
xlabel('Sum of pairwise distances')
pct = sum(test_statistic > bkg)/numel(bkg);
title_str = sprintf(['Objective with respect to random background\n',...
    'Evaluting lines: %s\n',...
    'Objective type is: %s\n',...
    'Objective = %.2f\n',...
    'nperms = %d\n',...
    'Percentile rank = %.3f'],...
    inputname(2),args.objective,test_statistic,args.nperms,pct);
title(title_str,'Interpreter','none')
grid on
legend show

end

function d = compute_total_diversity(ds, idx)
%Computes the sum of pairwise distances between objects

d = sum(tri2vec(ds.mat(idx,idx)));

end