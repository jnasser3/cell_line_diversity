function DI = compute_diversity_index_two(ds1, ds2, varargin)
%Given expression matrices for two cell lines, computes their diversity index.

pnames = {'corr_type',...
    'make_plot',...
    'show_plot',...
    'save_plot',...
    'plot_dir'};
dflts = {'pearson',...
    false,...
    false,...
    false,...
    '/cmap/projects/cell_line_diversity/analysis/pairwise_diversity_index/fig'};
args = parse_args(pnames, dflts, varargin{:});

% Data validation
assert(isequal(ds1.rid,ds2.rid), 'Gene space not the same')

% Order the perts identically b/w ds1 and ds2
[ds1,ds2] = order_data_sets(ds1, ds2);

%compute bioactivity weights
max_ccq_75 = max(...
        [cell2mat(ds1.cdesc(:,ds1.cdict('distil_cc_q75')))
        cell2mat(ds2.cdesc(:,ds1.cdict('distil_cc_q75')))]...
        );
bioactivity_mat = [...
    cell2mat(ds1.cdesc(:,ds1.cdict('distil_cc_q75'))),...
    cell2mat(ds2.cdesc(:,ds2.cdict('distil_cc_q75'))),...
    zeros(size(ds1.mat,2),1)
    ];
%bioa = max(bioactivity_mat,[],2)/max_ccq_75;
bioa = max(bioactivity_mat,[],2);

%compute pairwise correlations
corrs = pwcorr(ds1.mat,ds2.mat);

%compute diversity index
num_perts = size(ds1.mat,2);
DI = sum(bioa .* (1 - corrs)) / num_perts;

%old calculation
% for ii = 1:num_perts
%     bioa(ii) = max([cell2mat(ds1.cdesc(ii,ds1.cdict('distil_cc_q75')))
%         cell2mat(ds2.cdesc(ii,ds1.cdict('distil_cc_q75')))
%         0])/max_ccq_75;
%     corrs(ii) = corr(ds1.mat(:,ii),ds2.mat(:,ii),...
%                     'type',args.corr_type);
%     DI = DI + bioa(ii)*(1 - corrs(ii));
% end

if args.make_plot
    if args.show_plot
        figure;
    else
        figure('Visible','Off');
    end
    namefig(strcat(ds1.cdesc{1,ds1.cdict('cell_id')},'_',ds2.cdesc{1,ds2.cdict('cell_id')}));
    scatter(bioa,corrs,[],bioa .* (1 - corrs))
    axis([0 1 0 1])
    colormap(cool)
    colorbar
    caxis([0, 1])
    xlabel('Max Bioactivity')
    ylabel('Correlation')
    title_str = sprintf('Bioactivity vs correlation for perts in %s and %s,\n Num perts is: %d \n Diversity index is: %.3f',...
        ds1.cdesc{1,ds1.cdict('cell_id')},...
        ds2.cdesc{1,ds2.cdict('cell_id')},...
        num_perts,...
        DI);
    title(title_str)
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