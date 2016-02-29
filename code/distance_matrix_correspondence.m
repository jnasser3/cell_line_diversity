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
%        objects: cell array or grp of objects to include.
%        corr_type: Measure of the correlation between ds1 and
%                   ds2 distances. Default is spearman.
%        method: How to compute significance. Default mantel.
%        nperms: Number of permutations for significance test.
%        make_plot: Boolean. Make/show plots. Default false.
%        show_label: Boolean. Show object label on plots. Default false.
%
% Output: rho: the correlation between the distances.
%         pval: pvalue of this correlation
%         feat_vec: vector of distances in d1
%         div_vec: vector of distances in d2

params = {'ds1_id_field',...
          'ds2_id_field',...
          'objects',...
          'corr_type',...
          'method',...
          'nperms',...
          'make_plot',...
          'show_label'};
dflts = {'cid',...
         'cid',...
         '',...
         'Spearman',...
         'mantel',...
         5000,...
         false,...
         false};
args = parse_args(params,dflts,varargin{:});

% Get the ds from gctx if necessary
if ischar(ds1)    
    ds1 = parse_gctx(ds1);
end
if ischar(ds2)   
    ds2 = parse_gctx(ds2);
end

%make the cid in each ds equal to the common identifier
if ~strcmp(args.ds1_id_field,'cid')
    ds1.cid = ds1.cdesc(:,ds1.cdict(args.ds1_id_field));
    ds1.rid = ds1.cid;
end
if ~strcmp(args.ds2_id_field,'cid')
    ds2.cid = ds2.cdesc(:,ds2.cdict(args.ds2_id_field));
    ds2.rid = ds2.cid;
end

common_objects = intersect(ds1.cid, ds2.cid);
%get the cell lines to include
if isempty(args.objects)
    args.objects = common_objects;
else
    try
        args.objects = parse_grp(args.objects);
    end
end

%Subset each ds to common lines
%Note that ds_slice will order the structs alphabetically.
%common_lines = intersect(feature_ds.cid, diversity_ds.cid);
common_objects = intersect(common_objects,args.objects);
ds1_slice = ds_slice(ds1,'rid',common_objects,'cid',common_objects);
ds2_slice = ds_slice(ds2,'rid',common_objects,'cid',common_objects);
fprintf('ds1 contains %d objects \n',numel(ds1.cid));
fprintf('ds2 contains %d objects \n',numel(ds2.cid));
fprintf('There are %d objects in common \n',numel(common_objects));

assert(isequal(ds1_slice.cid,ds2_slice.cid),'ds cids not same order!')
assert(isequal(ds1_slice.rid,ds2_slice.rid),'ds rids not same order!')

num_objects = numel(common_objects);

%compute correspondence
switch args.method
    case 'naive'
        ds1_vec = tri2vec(ds1_slice.mat);
        ds2_vec = tri2vec(ds2_slice.mat);
        [rho, pval] = corr(ds1_vec,ds2_vec,...
            'type',args.corr_type,...
            'tail','right');
        
    case 'mantel'
        ds1_vec = tri2vec(ds1_slice.mat);
        ds2_vec = tri2vec(ds2_slice.mat);
        rho = corr(ds1_vec,ds2_vec,...
            'type',args.corr_type);
        
        %permutation null
        null0 = zeros(1,args.nperms);
        for ii = 1:args.nperms
            idx = randperm(num_objects);
            this_div_vec = tri2vec(ds2_slice.mat(idx,idx));
            null0(ii) = corr(ds1_vec,this_div_vec,...
                'type',args.corr_type);
        end
        
        pval = nnz(null0 > rho)/args.nperms;
        
        if args.make_plot
            figure;
            histogram(null0,'DisplayName','Permuted Null')
            hold on
            vline(rho,'r-')
            xlabel('Correlation')
            ylabel('Frequency')
            title(sprintf('Test statistic vs null for correspondence of distance matrices.\nnperms = %d, pval = %.4f',...
                args.nperms,...
                pval));
        end   
end

%Plotting
if args.make_plot
    figure;
    scatter(ds1_vec,ds2_vec, 250, '.')
    xlabel(inputname(1),'Interpreter','none')
    ylabel(inputname(2),'Interpreter','none')
    grid on;
    title_str = sprintf('N = %d cell lines, %d comparisons \n%s correlation is %.2f, nperms = %d, pval = %0.4f \n pval method = %s', ...
        numel(common_objects),...
        nchoosek(numel(common_objects),2),...
        args.corr_type,...
        rho,...
        args.nperms,...
        pval,...
        args.method);
    title(title_str)
    
    if args.show_label
        labels = mk_labels(ds1_slice.rid);
        text(double(ds1_vec),double(ds2_vec),labels,'Interpreter','none');
    end
end
end

function labels = mk_labels(ids)
%Given a set of n labels, return an choosek(n,2) set of labels by
%concatenating pairs of labels in the same ordering as produced by tri2vec

n = numel(ids);
labels = cell(nchoosek(n,2),1);

kk=1;
for ii = 2:n
    for jj = 1:(ii-1)
        labels(kk) = strcat(ids(jj),'_',ids(ii));
        kk = kk+1;
    end
end

end