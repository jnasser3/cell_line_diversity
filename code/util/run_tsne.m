function Y = run_tsne(ds,varargin)
%Run tsne and make highlight plots
%
%Input
%       ds: a square struct of pairwise distances
%       perplexity: Perplexity for tsne. Default 10
%       exclude: A set of id's to remove from the analysis
%       highlight: A set of id's to highlight in the plots
%       highlight_label: Highlight label for plotting.
%
%Output

params = {'perplexity',...
    'exclude',...
    'include',...
    'highlight',...
    'highlight_label'};
dflts = {10,...
    '',...
    '',...
    '',...
    ''};
args = parse_args(params,dflts,varargin{:});

highlight = get_array_input(args.highlight,'');
    
%Subset ds if needed
exclude = get_array_input(args.exclude,'');
include = get_array_input(args.include,ds.cid);
objects = setdiff(include,exclude);
ds = ds_slice(ds,'cid',objects,'rid',objects);

%Run tsne
Y = tsne_d(ds.mat,[],2,args.perplexity);

%plotting
if ~isempty(highlight)
    figure;
    idx = ismember(ds.rid,highlight);
    scatter(Y(~idx,1),Y(~idx,2),36,[.7,.7,.7],'filled');
    hold on
    scatter(Y(idx,1),Y(idx,2),100,'r','filled');
    text(Y(idx,1),Y(idx,2),...
        ds.rid(idx),...
        'FontSize',10,...
        'Color','r',...
        'horizontal','left',...
        'vertical','bottom');
    grid on
    title_str = sprintf(['t-SNE Plot of pairwise distances\n'...
        'Dataset is: %s\n'...
        'Highlight is: %s'],...
        inputname(1),...
        args.highlight_label);
    title(title_str,'Interpreter','none')
else
    figure;
    scatter(Y(:,1),Y(:,2),'k','filled');
    grid on
    title('t-SNE Plot of pairwise distances')
end


