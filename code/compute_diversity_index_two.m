function DI = compute_diversity_index_two(ds1, ds2, varargin)
%Given expression matrices for two cell lines, computes their diversity index.

pnames = {'corr_type',...
    'make_plot',...
    'save_plot',...
    'plot_dir'};
dflts = {'pearson',...
    false,...
    false,...
    '.'};
args = parse_args(pnames, dflts, varargin{:});

% Data validation
assert(isequal(ds1.rid,ds2.rid), 'Gene space not the same')

% Order the perts identically b/w ds1 and ds2
[ds1,ds2] = order_data_sets(ds1, ds2);

max_ccq_75 = max(...
        [cell2mat(ds1.cdesc(:,ds1.cdict('distil_cc_q75')))
        cell2mat(ds2.cdesc(:,ds1.cdict('distil_cc_q75')))]...
        );

num_perts = numel(ds1.cid);
bioa = zeros(1,num_perts);
corrs = zeros(1,num_perts);

DI = 0;
for ii = 1:num_perts
    bioa(ii) = max([cell2mat(ds1.cdesc(ii,ds1.cdict('distil_cc_q75')))
        cell2mat(ds2.cdesc(ii,ds1.cdict('distil_cc_q75')))
        0])/max_ccq_75;
    corrs(ii) = corr(ds1.mat(:,ii),ds2.mat(:,ii),...
                    'args.corr_type');
    DI = DI + bioa(ii)*(1 - corrs(ii));
end

if args.make_plot
    figure;
    namefig(strcat(ds1.cdesc{1,ds1.cdict('cell_id')},'_',ds2.cdesc{1,ds2.cdict('cell_id')}))
    scatter(bioa,corrs,[],bioa .* (1 - corrs))
    axis([0 1 0 1])
    colormap(cool)
    colorbar
    caxis([0, 1])
    xlabel('Max Bioactivity')
    ylabel('Correlation')
    title_str = sprintf('Bioactivity vs correlation for perts in %s and %s,\n Diversity index is: %.2f',...
        ds1.cdesc{1,ds1.cdict('cell_id')},...
        ds2.cdesc{1,ds2.cdict('cell_id')},...
        DI);
    title(title_str)
    
    if args.save_plot
        savefigures('out',args.plot_dir,...
                    'closefig',true)
    end
end

end

function [ds1, ds2] = order_data_sets(ds1, ds2)
    %order data sets by pert_id_dose_time combo
    [ds1, ds2] = subset_data_sets(ds1,ds2);
    
    %keep ds1 fixed. Order ds2 in the same order as ds1
    [~,idx] = orderas(ds1.cid,ds2.cid);
    ds2 = ds_order_meta(ds2, 'column', ds1.cid);
    ds2.mat = ds2.mat(:,idx);
end

function [ds1, ds2] = subset_data_sets(ds1, ds2)
    %subset data sets to common sig ids
    ids = intersect(ds1.cid, ds2.cid);
    ds1 = ds_slice(ds1, 'cid', ids);
    ds2 = ds_slice(ds2, 'cid', ids);
end