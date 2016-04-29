function [Y,STRESS] = run_mds(ds,varargin)
%Run multidimensional scaling and make plots
%
%Input
%       ds: a square struct of pairwise distances
%       exclude: A set of id's to remove from the analysis
%       highlight: A set of id's to highlight in the plots
%       highlight_label: Highlight label for plotting.
%       mds_type: classical or non-classical mds. runs cmdscale or mdscale.
%               Default classical.
%
%Output
%       [Y,STRESS], see help mdscale

params = {'exclude',...
    'include',...
    'highlight',...
    'highlight_label',...
    'mds_type'};
dflts = {'',...
    '',...
    '',...
    '',...
    'classical'};
args = parse_args(params,dflts,varargin{:});

highlight = get_array_input(args.highlight,'');
    
%Subset ds if needed
exclude = get_array_input(args.exclude,'');
include = get_array_input(args.include,ds.cid);
objects = setdiff(include,exclude);
ds = ds_slice(ds,'cid',objects,'rid',objects);

%Run MDS
switch args.mds_type
    case 'non-classical'
        [Y,STRESS] = mdscale(double(ds.mat),2);
    case 'classical'
        [Y,STRESS] = cmdscale(double(ds.mat),2);
end

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
    title_str = sprintf(['MDS Plot of pairwise distances\n'...
        'Dataset is: %s\n'...
        'Highlight is: %s'],...
        inputname(1),...
        args.highlight_label);
    title(title_str,'Interpreter','none')
else
    figure;
    scatter(Y(:,1),Y(:,2),'k','filled');
    grid on
    title('MDS Plot of pairwise distances')
end


