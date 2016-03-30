
ds1 = cpc006_pdi;
ds2 = ccle_cosine;

nperms = 10000;
corr_type = 'spearman';

%% Subset data sets to common object id's
[ds1,ds2,common_objects] = subset_ds_common_ids(ds1,ds2);

assert(isequal(ds1.cid,ds2.cid),'ds cids not same order!')
assert(isequal(ds1.rid,ds2.rid),'ds rids not same order!')
num_objects = numel(ds1.rid);

%% Mantel test
ds1_vec = tri2vec(ds1.mat);
mantel_null = zeros(1,nperms);
for ii = 1:nperms
    idx = randperm(num_objects);
    this_div_vec = tri2vec(ds2.mat(idx,idx));
    mantel_null(ii) = corr(ds1_vec,this_div_vec,...
        'type',corr_type);
end

%% Naive test
ds1_vec = tri2vec(ds1.mat);
ds2_vec = tri2vec(ds2.mat);
naive_null = zeros(1,nperms);
for ii = 1:nperms
    idx = randperm(numel(ds1_vec));
    this_div_vec = ds2_vec(idx);
    naive_null(ii) = corr(ds1_vec,this_div_vec,...
        'type',corr_type);
end

%% make pdfs
numpoints = 2^14;
[~,fmantel,xmantel] = kde(mantel_null,numpoints,...
    min(mantel_null),max(mantel_null));
[~,fnaive,xnaive] = kde(naive_null,numpoints,...
    min(naive_null),max(naive_null));

%% Plotting
figure;
plot(xmantel,fmantel,'-r','DisplayName','Mantel')
hold on
plot(xnaive,fnaive,'-b','DisplayName','Naive')
grid on
xlabel('Correlation')
ylabel('pdf')
title_str = sprintf(['Null correlations of upper triangles of distance matrices\n',...
    'Distance matrices are: %s and %s\n nperms = %d'])
title(title_str,'Interpreter','none')
legend('show')
shg