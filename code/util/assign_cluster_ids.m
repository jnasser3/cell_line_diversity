function cluster_id = assign_cluster_ids(ds,exemplar,varargin )
%Assign cluster ids to ds
%
%Input:
%       ds: square struct of pairwise distances
%       exemplar: List of exemplars to assign from which to assign cluster
%           memberships
%       exclude: array of ids to exclude from analysis
%       write_table: Boolean. Whether to write membership id's to
%           filesystem. Default false.
%
%Output
%       cluster_id: array of cluster memberships

params = {'exclude',...
    'write_table'};
dflts = {'',...
    false};
args = parse_args(params,dflts,varargin{:});

%Subset ds if needed
exclude = get_array_input(args.exclude,'');
if ~isempty(exclude)
    objects = setdiff(ds.rid,exclude);
    ds = ds_slice(ds,'cid',objects,'rid',objects);
end

exemplar_idx = ismember(ds.rid,exemplar);
ds_exemplar = ds.rid(exemplar_idx);

for ii = 1:numel(ds.rid) 
%     ds.rid(ii)
%     ds.mat(ii,exemplar_idx)
    [~,cluster_idx] = min(ds.mat(ii,exemplar_idx));
    cluster_id(ii) = ds_exemplar(cluster_idx);
    ds_exemplar(cluster_idx);
end

if args.write_table
    T = cell2table([ds.rid, cluster_id'], 'VariableNames',{'cell_id','cluster_id'});
    writetable(T,'/cmap/projects/cell_line_diversity/analysis/ccle_rna_seq_viz/tsne/cluster_ids.txt',...
        'Delimiter','\t')
end

end

