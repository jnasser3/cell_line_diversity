function [DI,all_corrs] = compute_diversity_index_two(ds1, ds2, args)
%Given expression matrices for two cell lines, computes their diversity index.

% pnames = {'metric',...
%     'kde_method',...
%     'make_plot',...
%     'show_plot',...
%     'save_plot',...
%     'plot_dir'};
% dflts = {'rel_bioa_corr',...
%     'matlab',...
%     false,...
%     false,...
%     false,...
%     '/cmap/projects/cell_line_diversity/analysis/pairwise_diversity_index/core_ts/rel_bioa_corr/fig'};
% args = parse_args(pnames, dflts, varargin{:});

% Data validation
assert(isequal(ds1.rid,ds2.rid), 'Gene space not the same')

% Order the perts identically b/w ds1 and ds2
[ds1,ds2] = order_data_sets(ds1, ds2);

assert(isequal(ds1.cdesc(:,ds1.cdict('pert_id_dose_time')),...
    ds2.cdesc(:,ds2.cdict('pert_id_dose_time'))), 'Sample space not ordered the same')

% Load bioactivity background
load /cmap/projects/cell_line_diversity/data/ts_cp_ccq75.mat

switch lower(args.metric)
    case 'basic_bioa_corr'
        DI = basic_bioa_corr(ds1,ds2,args);
    case 'rel_bioa_corr'
        [DI, all_corrs] = rel_wtd_corr(ds1,ds2,ts_cp_ccq75,args);
    case 'cov_fro'
        DI = covariance_frobenius(ds1,ds2);  
    case 'cov_eig'
end

end

function [ds1, ds2] = order_data_sets(ds1, ds2)
    %order data sets by pert_id_dose_time combo
    ds1 = remove_duplicates(ds1,'pert_id_dose_time');
    ds2 = remove_duplicates(ds2,'pert_id_dose_time');
    [ds1, ds2] = subset_data_sets(ds1,ds2);

    %Order both datasets alphabetically
    %[~,idx] = orderas(ds1.cdesc(:,ds1.cdict('pert_id_dose_time')),ds2.cdesc(:,ds2.cdict('pert_id_dose_time')));
    [~,idx1] = sort(ds1.cdesc(:,ds1.cdict('pert_id_dose_time')));
    [~,idx2] = sort(ds2.cdesc(:,ds2.cdict('pert_id_dose_time')));
    
    ds1.cdesc = ds1.cdesc(idx1,:);
    ds1.mat = ds1.mat(:,idx1);
    
    ds2.cdesc = ds2.cdesc(idx2,:);
    ds2.mat = ds2.mat(:,idx2);
end

function [ds1, ds2] = subset_data_sets(ds1, ds2)
    %subset data sets to common combos
    common_combos = intersect(ds1.cdesc(:,ds1.cdict('pert_id_dose_time')), ds2.cdesc(:,ds2.cdict('pert_id_dose_time')));
    idx1 = ismember(ds1.cdesc(:,ds1.cdict('pert_id_dose_time')),common_combos);
    idx2 = ismember(ds2.cdesc(:,ds2.cdict('pert_id_dose_time')),common_combos);
    ds1 = ds_slice(ds1, 'cidx', find(idx1));
    ds2 = ds_slice(ds2, 'cidx', find(idx2));
end

function out_ds = remove_duplicates(ds,field)
%Given a ds, returns ds_out such that ds_out.cdesc(,ds_out.cdict('field'))
%contains only unique entries

dups = duplicates(ds.cdesc(:,ds.cdict(field)));
idx_to_remove = [];
for ii = 1:numel(dups)
    temp = find(strcmp(ds.cdesc(:,ds.cdict(field)),dups(ii)));
    bad_idx = temp(2:end);
    idx_to_remove = [idx_to_remove; bad_idx];
end

out_ds = ds_slice(ds,'cid',ds.cid(idx_to_remove),'exclude_cid',true);

end

function DI = basic_bioa_corr(ds1,ds2,args)
%compute diversity index using a basic sum over max(bioa)*(1 - corr)
%metric.

%compute bioactivity weights
bioactivity_mat = [...
    cell2mat(ds1.cdesc(:,ds1.cdict('distil_cc_q75'))),...
    cell2mat(ds2.cdesc(:,ds2.cdict('distil_cc_q75'))),...
    zeros(size(ds1.mat,2),1)
    ];
bioa = max(bioactivity_mat,[],2);

%compute pairwise correlations
pwcorrs = pwcorr(ds1.mat,ds2.mat);

%compute diversity index
num_perts = size(ds1.mat,2);
DI = sum(bioa .* (1 - pwcorrs)) / num_perts;

if args.make_plot
    if args.show_plot
        h = figure;
    else
        h = figure('Visible','Off');
    end
    name = strcat(ds1.cdesc{1,ds1.cdict('cell_id')},...
        '_',ds2.cdesc{1,ds2.cdict('cell_id')});
    namefig(h,name);
    scatter(bioa,pwcorrs,[],bioa .* (1 - pwcorrs))
    axis([0 1 -1 1])
    colormap(cool)
    colorbar
    caxis([0, 1])
    grid on
    xlabel('Max Bioactivity')
    ylabel('Correlation')
    title_str = sprintf('Bioactivity vs correlation for perts in %s and %s,\n Num perts is: %d \n Diversity index is: %.3f',...
        ds1.cdesc{1,ds1.cdict('cell_id')},...
        ds2.cdesc{1,ds2.cdict('cell_id')},...
        num_perts,...
        DI);
    title(title_str)
    
    if args.save_plot
        saveas(h,strcat(fullfile(args.plot_dir,name),'.png'),'png')
        close;
    end
end

end

function DI = covariance_frobenius(ds1,ds2)
% Compute the frobenius norm of the difference in the covariance matrices.

cov1 = cov(ds1.mat');
cov2 = cov(ds2.mat');
DI = norm(cov1 - cov2, 'fro');
end

function [DI, all_corrs] = rel_wtd_corr(ds1,ds2,bioa_bkg,args)

%get corr contributions
all_corrs = fastcorr(ds1.mat,ds2.mat,...
    'type',args.corr_type);
[matched_corrs, unmatched_corrs] = mat2vec_diag(corrs);
corr_contribution = xform_corr_diversity(matched_corrs,unmatched_corrs,...
    'cdf',args.kde_method);

%get bioactiviy cdfs and convert to cdf contributions
bioa1 = cell2mat(ds1.cdesc(:,ds1.cdict('distil_cc_q75')));
bioa2 = cell2mat(ds2.cdesc(:,ds2.cdict('distil_cc_q75')));
bioa1_rel = xform_bioa_diversity(bioa1, bioa_bkg);
bioa2_rel = xform_bioa_diversity(bioa2, bioa_bkg);

bioa_contribution = max([bioa1_rel, bioa2_rel], [], 2);
DI = sum(corr_contribution .* bioa_contribution) / sum(bioa_contribution);
sum(bioa_contribution);

if args.make_plot
    
    cell1 = ds1.cdesc{1,ds1.cdict('cell_id')};
    cell2 = ds2.cdesc{1,ds2.cdict('cell_id')};

    mk_di_scatter2(max([bioa1, bioa2], [], 2),...
        diag(all_corrs),...
        corr_contribution .* bioa_contribution,...
        DI,...
        cell1,...
        cell2,...
        args);
    
%     mk_di_scatter3(bioa1,...
%         bioa2,...
%         diag(all_corrs),...
%         corr_contribution .* bioa_contribution,...
%         DI,...
%         cell1,...
%         cell2,...
%         args);
    
    mk_corr_plot(all_corrs,...
        cell1,...
        cell2,...
        args);

end
end

function h = mk_di_scatter2(bioa,corrs,DI_contribution,DI,cell1,cell2,args)
    if args.show_plot
        h = figure;
    else
        h = figure('Visible','Off');
    end
    name = strcat(cell1,...
        '_',cell2,...
        '_scatter2');
    namefig(h,name);
    scatter(bioa,corrs,[],DI_contribution)
    axis([0 1 -1 1])
    colormap(cool)
    c = colorbar;
    c.Label.String = 'Contribution to Diversity';
    caxis([0, 1])
    grid on
    xlabel('Maximum Replicate Reproducibility')
    ylabel('Correlation')
    title_str = sprintf(['Reproducibility vs correlation for signatures in %s and %s \n',...
        'Number of signatures = %d \n',...
        'Correlation type: %s \n',...
        'Diversity index = %.3f'],...
        cell1,...
        cell2,...
        numel(bioa),...
        args.corr_type,...
        DI);
    title(title_str,'interpreter','none')
    
    if args.save_plot
        saveas(h,strcat(fullfile(args.plot_dir,name),'.png'),'png')
        close;
    end
end

function h = mk_di_scatter3(bioa1,bioa2,corrs,DI_contribution,DI,cell1,cell2,args)
    if args.show_plot
        h = figure;
    else
        h = figure('Visible','Off');
    end
    name = strcat(cell1,...
        '_',cell2,...
        '_scatter3');
    namefig(h,name);
    scatter3(bioa1,bioa2,corrs,[],DI_contribution)
    axis([0 1 0 1 -1 1])
    colormap(cool)
    colorbar
    caxis([0, 1])
    grid on
    xlabel([cell1 ' ccq75'])
    ylabel([cell2 ' ccq75'])
    zlabel('Correlation')
    title_str = sprintf(['Reproducibility vs correlation for signatures in %s and %s\n',...
        'Number of signatures is: %d \n'...
        'metric: %s \n'...
        'Diversity index is: %.3f'],...
        cell1,...
        cell2,...
        numel(bioa1),...
        args.metric,...
        DI);
    title(title_str,'Interpreter','none')
    
    if args.save_plot
        saveas(h,strcat(fullfile(args.plot_dir,name),'.png'),'png')
        close;
    end
end

function h = mk_corr_plot(all_corrs,cell1,cell2,args)

    if args.show_plot
        h = figure;
    else
        h = figure('Visible','Off');
    end
    
    name = strcat(cell1,...
        '_',cell2,...
        '_corr');
    namefig(h,name);
    
    [fnull, xnull] = ksdensity(tri2vec(all_corrs));
    [fmatch, xmatch] = ksdensity(diag(all_corrs));
    plot(xnull,fnull,'DisplayName','Null - unmatched pairs');
    hold on
    plot(xmatch,fmatch,'DisplayName','Matched pairs');
    title_str = sprintf(['Correlation of signatures in %s and %s\n',...
        'Correlation type is: %s\n',...
        'Number of signatures = %d'],...
        cell1,cell2,...
        args.corr_type,...
        size(all_corrs,1));
    title(title_str)
    grid on
    axis([-1 1 0 6])
    xlabel('Correlation')
    ylabel('pdf')
    legend('show')
    
    if args.save_plot
        saveas(h,strcat(fullfile(args.plot_dir,name),'.png'),'png')
        close;
    end
end