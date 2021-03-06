function [test_statistic, bkg] = evaluate_total_diversity(ds, lines, varargin)
% test_statistic = evaluate_total_diversity(ds, lines, varargin)
%
% Start with a universe of objects and pairwise distances between them.
% Given a subset of these objects, computes the sum of their mutual 
% pairwise distances and compares this to random choices.
%
% Input
%       ds: A square struct of pairwise distances
%       lines: An array of cell arrays. Each array contains a set of lines to evaluate.
%           Each array must have the same number of elements.
%           Note: If only a single array is present it should be passed in
%           as an array containing a single array. Ie
%           evaluate_total_diversity(ds, {lines_array})
%       ds_type: String. Either 'similarity' or 'distance', denoting
%           ds contains either similarities or distances. Default
%           'similarity'
%       include: The set of lines to include in the analysis
%       lines_labels: An array of labels for lines for plotting purposes.
%       objective: The objective function to evaluate. Options are divr
%           (sum of distances among chosen lines) or repr (sum of minimum
%           distances from all lines to nearest chosen line. repr is optimized
%           by affinity propogation. Default is 'repr'
%       exclude: A set of id's to remove from the analysis
%       fixed_lines: A subset of objects to include in the calculation but
%                   not to evaluate
%       nperms: Number of permutations to compare against. Default 1000
%       make_plot: Boolean. Default true
%
% Output
%       test_statistic: array of objective loss for each array in lines
%       bkg: Objective loss under random selection of exemplars.

params = {'include',...
    'fixed_lines',...
    'objective',...
    'exclude',...
    'lines_labels',...
    'ds_type',...
    'nperms',...
    'make_plot'};
dflts = {'',...
    '',...
    'repr',...
    '',...
    '',...
    'similarity',...
    1000,...
    true};
args = parse_args(params,dflts,varargin{:});

%Check ds_type
assert(ismember(args.ds_type,{'distance','similarity'}),'Invalid ds_type entry')

%Make sure all arrays in lines have the same number of elements.
sizes = cellfun('length',lines);
assert(range(sizes) == 0, 'All arrays must have the same number of elements');

%only need to check this if we are making plots.
if args.make_plot
    assert(length(lines) == length(args.lines_labels), 'Must have same number of lines and labels')
end

%Subset ds if needed
exclude = get_array_input(args.exclude,'');
if ~isempty(exclude)
    objects = setdiff(ds.rid,exclude);
    ds = ds_slice(ds,'cid',objects,'rid',objects);
end

%Subset ds if needed
include = get_array_input(args.include,ds.cid);
include = intersect(include,ds.cid);
ds = ds_slice(ds,'cid',include,'rid',include,'ingore_missing',true);

%get fixed cell lines
fixed_lines = get_array_input(args.fixed_lines,'');
non_fixed_lines = setdiff(ds.rid,fixed_lines);

%compute the objectives
test_statistic = zeros(1,length(lines));
for ii = 1:length(lines)
    this_lines = get_array_input(lines{ii},'');
    switch lower(args.objective)
        case 'divr'
                test_statistic(ii) = compute_total_diversity(ds,[this_lines; fixed_lines]);
        case 'repr'
                test_statistic(ii) = evaluate_pmedian_loss(ds, [this_lines; fixed_lines], args.ds_type);
    end
end

%compute a background
n = numel(setdiff(ds.rid,fixed_lines));
k = numel(lines{1});
bkg = zeros(1,args.nperms);
for ii = 1:args.nperms
    temp = randperm(n,k);
    this_lines = non_fixed_lines(temp);

    switch lower(args.objective)
        case 'divr'
            bkg(ii) = compute_total_diversity(ds,[this_lines; fixed_lines]);  
        case 'repr'
            bkg(ii) = evaluate_pmedian_loss(ds, [this_lines; fixed_lines],  args.ds_type);
    end
end

%plot
if args.make_plot
    figure;
    histogram(bkg,'DisplayName','Background');
    hold on
    g = gca;
    for ii = 1:length(lines)
        plot([test_statistic(ii) test_statistic(ii)], g.YLim,...
            'DisplayName', args.lines_labels{ii},...
            'LineWidth',2)
        hold on
    end
    xlabel('Objective')
    ylabel('Frequency')
    title_str = sprintf(['Objective with respect to random background\n',...
        'Pairwise Distance Dataset: %s\n',...
        'Num Exemplar: %d\n',...
        'Objective type is: %s\n',...
        'nperms = %d'],...
        inputname(1),numel([lines{1}; fixed_lines]),args.objective,args.nperms);
    title(title_str,'Interpreter','none')
    grid on
    legend('show')
    leg = legend;
    set(leg, 'Interpreter', 'None')
end

end

function d = compute_total_diversity(ds, lines)
%Computes the sum of pairwise distances between objects

idx = ismember(ds.rid,lines);
d = sum(tri2vec(ds.mat(idx,idx)));

%Average for easier interpretation
d = d / nchoosek(nnz(idx),2);

end

