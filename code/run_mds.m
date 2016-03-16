function run_mds(ds,varargin)
%
%Run multidimensional scaling and make plots
%
%Input
%       ds: a square struct of pairwise distances
%       exclude: A set of id's to remove from the analysis
%       highlight: A set of id's to highlight in the plots
%
%Output
%       

params = {'exclude',...
    'highlight'};
dflts = {'',...
    ''};
args = parse_args(params,dflts,varargin{:});

highlight = get_array_input(args.highlight,'');
    
%Subset ds if needed
exclude = get_array_input(args.exclude,'');
if ~isempty(exclude)
    objects = setdiff(ds.rid,exclude);
    ds = ds_slice(ds,'cid',objects,'rid',objects);
end

%Run MDS
Y = mdscale(ds.mat,2);

%plotting
if ~isempty(highlight)
    figure;
    idx = ismember(ds.rid,highlight);
    scatter(Y(~idx,1),Y(~idx,2),'k','filled');
    hold on
    scatter(Y(idx,1),Y(idx,2),'r','filled');
    text(Y(idx,1),Y(idx,2),...
        ds.rid(idx),...
        'horizontal','left',...
        'vertical','bottom')
    grid on
    title('MDS Plot of pairwise distances')
else
    figure;
    scatter(Y(:,1),Y(:,2),'k','filled');
    grid on
    title('MDS Plot of pairwise distances')
end

