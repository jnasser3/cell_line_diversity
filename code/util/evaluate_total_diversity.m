function evaluate_total_diversity(ds, lines, varargin )
%
% Start with a universe of objects and pairwise distances between them.
% Given a subset of these objects, computes the sum of their  mutual 
% pairwise distances and compares this to random choices.
%
% Input
%       ds: A square struct of pairwise distances
%       exclude: A set of id's to remove from the analysis
%       lines: The set of objects to evaluate
%       fixed_lines: A subset of objects to include in the calculation but
%                   not to evaluate
%       nperms: Number of permutations to compare against. Default 1000
%

params = {'fixed_lines',...
    'exclude',...
    'nperms'};
dflts = {'',...
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
total_diversity = compute_total_diversity(ds,[lines_idx; fixed_lines_idx]);

%compute a background
n = numel(setdiff(ds.rid,fixed_lines));
k = numel(lines);
bkg = zeros(1,args.nperms);
for ii = 1:args.nperms
    temp = randperm(n,k);
    this_lines_idx = non_fixed_idx(temp);
    bkg(ii) = compute_total_diversity(ds,[this_lines_idx; fixed_lines_idx]);  
end

%plot
figure;
histogram(bkg,'DisplayName','Background')
hold on
vline(total_diversity,'r-');
xlabel('Sum of pairwise distances')
ylabel('Frequency')
pct = sum(total_diversity > bkg)/numel(bkg);
title_str = sprintf('Diversity objective with respect to background. Percentile rank is: %.2f', pct);
title(title_str)
grid on

end

function d = compute_total_diversity(ds, idx)
%Computes the sum of pairwise distances between objects

d = sum(tri2vec(ds.mat(idx,idx)));

end