function [rho, pval, ds1_vec, ds2_vec] = distance_matrix_correspondence(ds1,ds2,varargin)
%[rho, pval, feat_vec, div_vec] = distance_matrix_correspondence(ds1,ds2,varargin)
%
% Computes the correlation between the pairwise distances of a set of
% objects measured using two different distance metrics
%
% Input: ds1: an square struct where each entry. ds1(i,j) is the distance
%                   between objects i and j under distance 1
%        ds2: an square struct where each entry. ds2(i,j) is the distance
%                   between objects i and j under distance 2
%        ds1_id_field: The common identifier field. Default cid
%        ds2_id_field: The common identifier field. Default is cid
%        exclude: cell array or grp of objects to exclude.
%        corr_type: Measure of the correlation between ds1 and
%                   ds2 distances. Default is spearman.
%        method: How to compute significance. Default mantel.
%        nperms: Number of permutations for significance test. Default
%               1000
%        make_plot: Boolean. Make/show plots. Default false.
%        show_label: Boolean. Show object label on plots. Default false.
%        make_dist_graph: Boolean. Make/Show mds/tsne plots. Default false.
%        perplexity: Perplexity for tSNE. Default is 5
%
% Output: rho: the correlation between the distances.
%         pval: pvalue of this correlation
%         feat_vec: vector of distances in d1
%         div_vec: vector of distances in d2

params = {'ds1_id_field',...
          'ds2_id_field',...
          'exclude',...
          'corr_type',...
          'method',...
          'nperms',...
          'make_plot',...
          'show_label',...
          'perplexity',...
          'make_dist_graph'};
dflts = {'cid',...
         'cid',...
         '',...
         'Spearman',...
         'mantel',...
         1000,...
         false,...
         false,...
         5,...
         false};
args = parse_args(params,dflts,varargin{:});

%% Exclude objects if needed
exclude = get_array_input(args.exclude,'');
if ~isempty(exclude)
    objects1 = setdiff(ds1.rid,exclude);
    ds1 = ds_slice(ds1,'cid',objects1,'rid',objects1);
    
    objects2 = setdiff(ds2.rid,exclude);
    ds2 = ds_slice(ds2,'cid',objects2,'rid',objects2);
end

%% Subset data sets to common object id's
[ds1,ds2,common_objects] = subset_ds_common_ids(ds1,ds2,...
    'ds1_id_field', args.ds1_id_field,...
    'ds2_id_field', args.ds2_id_field);

assert(isequal(ds1.cid,ds2.cid),'ds cids not same order!')
assert(isequal(ds1.rid,ds2.rid),'ds rids not same order!')

num_objects = numel(common_objects);

%% Compute correspondence
switch args.method
    case 'naive'
        ds1_vec = tri2vec(ds1.mat);
        ds2_vec = tri2vec(ds2.mat);
        [rho, pval] = corr(ds1_vec,ds2_vec,...
            'type',args.corr_type,...
            'tail','right');
        
    case 'mantel'
        ds1_vec = tri2vec(ds1.mat);
        ds2_vec = tri2vec(ds2.mat);
        rho = corr(ds1_vec,ds2_vec,...
            'type',args.corr_type);
        
        %permutation null
        null0 = zeros(1,args.nperms);
        for ii = 1:args.nperms
            idx = randperm(num_objects);
            this_div_vec = tri2vec(ds2.mat(idx,idx));
            null0(ii) = corr(ds1_vec,this_div_vec,...
                'type',args.corr_type);
        end
        
        pval = nnz(null0 > rho)/args.nperms;
        
%         if args.make_plot
%             figure;
%             histogram(null0,'DisplayName','Permuted Null')
%             hold on
%             vline(rho,'r-')
%             xlabel('Correlation')
%             ylabel('Frequency')
%             title(sprintf('Test statistic vs null for correspondence of distance matrices.\nnperms = %d, pval = %.4f',...
%                 args.nperms,...
%                 pval));
%         end   
end

%% Plotting
if args.make_plot
    figure;
    scatter(ds1_vec,ds2_vec, 250, '.')
    xlabel(inputname(1),'Interpreter','none')
    ylabel(inputname(2),'Interpreter','none')
    grid on;
    title_str = sprintf('N = %d cell lines, %d comparisons \n%s correlation is %.2f, nperms = %d, pval = %0.3f \n pval method = %s', ...
        numel(common_objects),...
        nchoosek(numel(common_objects),2),...
        args.corr_type,...
        rho,...
        args.nperms,...
        pval,...
        args.method);
    title(title_str)
    
    if args.show_label
        labels = mk_pwlabels(ds1.rid);
        text(double(ds1_vec),double(ds2_vec),labels,'Interpreter','none');
    end
end

if args.make_dist_graph
    mk_tsne_d(ds1,inputname(1),args);
    mk_tsne_d(ds2,inputname(2),args);
    mk_mds(ds1,inputname(1),args);
    mk_mds(ds2,inputname(2),args);
end

end

function mk_tsne_d(ds,ds_name,args)

mapped = tsne_d(ds.mat,[],2,args.perplexity);

figure;
scatter(mapped(:,1),mapped(:,2))
hold on
grid on
if args.show_label
    text(mapped(:,1),mapped(:,2),ds.rid)
end
title(sprintf('tSNE - %s',ds_name),...
    'interpreter','none')
end

function mk_mds(ds,ds_name,args)

mapped = cmdscale(double(ds.mat),2);

figure;
scatter(mapped(:,1),mapped(:,2))
hold on
grid on
if args.show_label
    text(mapped(:,1),mapped(:,2),ds.rid)
end
title(sprintf('MDS - %s',ds_name),...
    'interpreter','none')
end


