function DI = compute_diversity_index_two(ds1, ds2, varargin)
%Given expression matrices for two cell lines, computes their diversity index.

pnames = {'metric',...
    'make_plot',...
    'show_plot',...
    'save_plot',...
    'plot_dir'};
dflts = {'basic_bioa_corr',...
    false,...
    false,...
    false,...
    ''};
args = parse_args(pnames, dflts, varargin{:});

% Data validation
assert(isequal(ds1.rid,ds2.rid), 'Gene space not the same')

% Order the perts identically b/w ds1 and ds2
[ds1,ds2] = order_data_sets(ds1, ds2);

switch lower(args.metric)
    case 'basic_bioa_corr'
        DI = basic_bioa_corr(ds1,ds2,args);
    case 'rel_wtd_corr'
        DI = rel_wtd_corr(ds1,ds2,args);
    case 'cov_fro'
        DI = covariance_frobenius(ds1,ds2);  
    case 'cov_eig'
end

end

function [ds1, ds2] = order_data_sets(ds1, ds2)
    %order data sets by pert_id_dose_time combo
    [ds1, ds2] = subset_data_sets(ds1,ds2);
    ds1 = remove_duplicates(ds1,'pert_id_dose_time');
    ds2 = remove_duplicates(ds2,'pert_id_dose_time');

    %keep ds1 fixed. Order ds2 in the same order as ds1
    [~,idx] = orderas(ds1.cdesc(:,ds1.cdict('pert_id_dose_time')),ds2.cdesc(:,ds2.cdict('pert_id_dose_time')));
    ds2 = ds_order_meta(ds2, 'column', ds1.cdesc(:,ds1.cdict('pert_id_dose_time')));
    ds2.mat = ds2.mat(:,idx);
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

function DI = rel_wtd_corr(ds1,ds2,args)

%get corr contributions
all_corrs = fastcorr(ds1.mat,ds2.mat);
corr_contributions = xform_corr_diversity(diag(all_corrs),tri2vec(all_corrs));

%get bioactiviy cdfs and convert to cdf contributions
bioa1 = cell2mat(ds1.cdesc(:,ds1.cdict('distil_cc_q75')));
bioa1 = rankorder(bioa1)/numel(bioa1);
bioa2 = cell2mat(ds1.cdesc(:,ds1.cdict('distil_cc_q75')));
bioa2 = rankorder(bioa2)/numel(bioa2);

max_bioa = max([bioa1, bioa2], [], 2);
DI = sum(corr_contributions .* max_bioa) / sum(max_bioa);

end