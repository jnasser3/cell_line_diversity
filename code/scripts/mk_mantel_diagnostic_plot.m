function [test_statistic, mantel_pval, naive_pval] = mk_mantel_diagnostic_plot(ds1,ds2,varargin)
% [test_statistic, mantel_pval, naive_pval] = mk_mantel_diagnostic_plot(ds1,ds2,varargin)
% Compute permutation pvalues for correlations of upper triangles of two
% distance matrices using naive and mantel tests. Make a nice plot to show
% the null distributions.
%
% Input:
%       ds1: square pairwise distance gct struct
%       ds2: square pairwise distance gct struct
%       nperms: number of permutations. Default 10000
%       corr_type: type of correlations. Default 'spearman'
%
% Output:
%       test_statistic: The correlation between the pairwise distance
%           matrices
%       mantel_pval: Pvalue from mantel permutation test
%       naive_pval: Pvalue from naive permutation test

params = {'nperms',...
    'corr_type'};
dflts = {10000,...
    'spearman'};
args = parse_args(params,dflts,varargin{:});

%% Subset data sets to common object id's
[ds1,ds2,~] = subset_ds_common_ids(ds1,ds2);

assert(isequal(ds1.cid,ds2.cid),'ds cids not same order!')
assert(isequal(ds1.rid,ds2.rid),'ds rids not same order!')
num_objects = numel(ds1.rid);

%% Compute statistic
ds1_vec = tri2vec(ds1.mat);
ds2_vec = tri2vec(ds2.mat);
test_statistic = corr(ds1_vec,ds2_vec,'type',args.corr_type);

%% Mantel test
mantel_null = zeros(1,args.nperms);
for ii = 1:args.nperms
    idx = randperm(num_objects);
    this_div_vec = tri2vec(ds2.mat(idx,idx));
    mantel_null(ii) = corr(ds1_vec,this_div_vec,...
        'type',args.corr_type);
end
mantel_pval = nnz(mantel_null > test_statistic)/(1 + args.nperms);

%% Naive test
naive_null = zeros(1,args.nperms);
for ii = 1:args.nperms
    idx = randperm(numel(ds1_vec));
    this_div_vec = ds2_vec(idx);
    naive_null(ii) = corr(ds1_vec,this_div_vec,...
        'type',args.corr_type);
end
naive_pval = nnz(naive_null > test_statistic)/(1 + args.nperms);


%% make pdfs
numpoints = 2^14;
min_lim = min([mantel_null naive_null]);
max_lim = max([mantel_null naive_null]);
[~,fmantel,xmantel] = kde(mantel_null,numpoints,...
    min_lim,max_lim);
[~,fnaive,xnaive] = kde(naive_null,numpoints,...
    min_lim,max_lim);

%% Plotting
figure;
plot(xmantel,fmantel,'-r','DisplayName','Mantel Null')
hold on
plot(xnaive,fnaive,'-k','DisplayName','Naive Null')
grid on
hold on
y=get(gca,'ylim');
hold on
plot([test_statistic test_statistic],y,'b:','DisplayName','Test Statistic')
xlabel('Correlation')
ylabel('pdf')
title_str = sprintf(['Null correlations of upper triangles of distance matrices\n',...
    'Distance matrices are: %s and %s\n',...
    'Number of elements = %d, Number of comparisons = %d\n',...
    'statistic = %.3f, nperms = %d, naive pval = %.3f, mantel pval = %.3f'],...
    inputname(1),inputname(2),...
    num_objects,nchoosek(num_objects,2),...
    test_statistic,args.nperms,naive_pval,mantel_pval);
title(title_str,'Interpreter','none')
legend('show')
shg

end