function flat_ds = flatten_pw_distance_matrix(ds,varargin)
%flat_ds = flatten_pw_distance_matrix(ds,varagin)
%
%Input: ds - an N x N dataset of pairwise distances
%       ids - an array or .grp of ids to include
%
%Output: flat_ds - an NC2 x 1 data set of pairwise distances where each row
%                  in flat_ds represents a pair of ids in ds.

params = {'ids'};
dflts = {''};
args = parse_args(params,dflts,varargin{:});


%% Subset ds to relevant ids
ids = ds.rid;

if iscell(args.ids)
    ids = args.ids;
elseif exist(args.ids,2)
    ids = parse_grp(args.ids);
end

ds = ds_slice(ds, 'rid', ids, 'cid', ids);

%% mk flattened ds
%n = numel(ds.rid);
%mat = zeros(nchoosek(n,2),1);
%rid = cell(nchoosek(n,2),1);

% kk = 1;
% for ii = 1:n
%     for jj = (ii+1):n
%         cell1 = ds.rid(ii);
%         cell2 = ds.rid(jj);
%         idx1 = strcmp(ds.rid,cell1);
%         idx2 = strcmp(ds.rid,cell2);
%         
%         mat(kk,:) = ds.mat(idx1,idx2);
%         rid(kk) = strcat(cell1,'_',cell2);
%         kk = kk + 1;
%     end
% end

flat_mat = tri2vec(ds.mat);
flat_rid = mk_pwlabels(ds.rid);

flat_ds = mkgctstruct(flat_mat, 'rid', flat_rid, 'cid', {'pw_distance'});



end

