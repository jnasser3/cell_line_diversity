function  mk_pdist_choice_diagnostic_plot(ds,chosen,varargin)
%Given a matrix of pairwise distances and a set of objects, makes some
%diagnostics to describe the relative diversity of the set.
%
%Input
%       ds: square struct of pairwise distances
%       chosen: A set of id's. The chosen elements
%       exclude: A set of id's to ignore from ds
%       known: A set of fixed id's

params = {'exclude',...
    'known'};
dflts = {'',...
    ''};
args = parse_args(params,dflts,varargin{:});

%% Exclude ids if necessary
exclude = get_array_input(args.exclude,'');
if ~isempty(exclude)
    objects = setdiff(ds.rid,exclude);
    ds = ds_slice(ds,'cid',objects,'rid',objects);
end

%% Get distances
chosen = get_array_input(chosen,'');
all_dists = tri2vec(ds.mat);
choice_idx = ismember(ds.rid,chosen);
chosen_dists = tri2vec(ds.mat(choice_idx,choice_idx));

%% Plot
mesh = 2^8;
MIN = min(min(ds.mat));
MAX = max(max(ds.mat));
[~,f_all,x_all] = kde(all_dists,mesh,MIN,MAX);
[~,f_chosen,x_chosen] = kde(chosen_dists,mesh,MIN,MAX);

figure;
plot(x_all, f_all, 'DisplayName', 'All Lines')
hold on
plot(x_chosen, f_chosen, 'DisplayName', 'Chosen Lines')
hold on
scatter(chosen_dists,zeros(1,numel(chosen_dists)),'DisplayName','Chosen Lines')
xlabel('Distances')
ylabel('Density')
grid on
legend show
title_str = sprintf(['Pairwise distances chosen vs all elements\n',...
    'Chosen Lines: %s'],...
    inputname(2));
title(title_str,'interpreter','none')


end

