function sim = compute_genex_pdist(ds,varargin)
% compute_genex_sim - given a genes x samples struct, compute the pairwise
% distance between pairs of samples.
% 
% [sim, excl] = compute_genex_sim(varargin)
% Returns:
%     sim       - a gctx with cid = rid = set of cell lines with element ij in mat given by the 
%                 distance between cell line i and j using the given metric in the given CCLE ds
%
% Parameters:
%     ds        - genes by cells struct
%     cells     - A cell array of cell lines.  
%     metric    - string; the distance function.  
%     num_pcs - Dimension of PCA to perform before computing distances.
%               Only supports Euclidean and mahalanbois distance. Default 30          
%     gene_space - .grp file or cell array. Space of genes in which to
%                  compute distances. Default is all genes in the ds.
%     norm      - normalize the stdev of all genes to 1. 

pnames = {'cells', ...
    'metric', ...
    'num_pcs',...
    'gene_space',...
    'norm'};
dflts = {{}, ...
    'euclidean', ...
    30,...
    '',...
    1};
args = parse_args(pnames, dflts, varargin{:});

%Subset to relevant cell lines and genes
cells = get_array_input(args.cells,ds.cid);
genes = get_array_input(args.gene_space,ds.rid);
ds = ds_slice(ds,'cid',cells,'rid',genes);

if args.norm
  ds = norm_data(ds,args);
end

sim = calc_gex_sim(ds, args);

end

function ds = norm_data(ds,args)
    switch args.norm_type
        case 'std'
            s = std(ds.mat, 0, 2);
        case 'mad'
            s = mad(ds.mat, 0, 2);
    end
    s(s == 0) = 1;
    ds.mat = ds.mat./repmat(s, 1, numel(ds.cid));
end

function ret = calc_gex_sim(ds, args)
  switch lower(args.metric)
    case 'euclidean'
        data = ds.mat';
        [~, pca_data] = pca(data);
        pca_data = pca_data(:,1:args.num_pcs);
        sim = squareform(pdist(pca_data, 'euclidean')); 
    
    case 'mahal'
        data = ds.mat';
        [~, pca_data] = pca(data);
        pca_data = pca_data(:,1:args.num_pcs);
        sim = squareform(pdist(pca_data, 'mahalanobis'));
        
    case 'cosine'
        sim = squareform(pdist(ds.mat', 'cosine'));
    
    case 'spearman'
        %Note: dont transpose matrix
        sim = 1 - fastcorr(ds.mat, 'type', 'spearman');
        
    case 'pearson'
        %Note: dont transpose matrix
        sim = 1 - fastcorr(ds.mat, 'type', 'pearson'); 
  end

  ret = mkgctstruct(sim, ...
      'rid', ds.cid, ...
      'cid', ds.cid, ...
      'rdesc', ds.cdesc, ...
      'cdesc', ds.cdesc, ...
      'chd', ds.chd, ...
      'rhd', ds.chd);
end