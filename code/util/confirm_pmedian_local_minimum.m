function confirm_pmedian_local_minimum(ds,exemplar,exemplar_to_check,varargin)
%Sanity check to make sure we are in a local minimum.
%
%Checks to see if we can get a better optimum by making a single exemplar
%exchange
%
%Input:
%       ds: Square struct of pairwise distances
%       exemplar: Known/fixed exemplars
%       exemplar_to_check: which exemplars should be checked for local
%       optima
%       exclude: List of cell lines to exclude from the analysis

params = {'exclude'};
dflts = {''};
args = parse_args(params,dflts,varargin{:});

%Subset ds if needed
exclude = get_array_input(args.exclude,'');
if ~isempty(exclude)
    objects = setdiff(ds.rid,exclude);
    ds = ds_slice(ds,'cid',objects,'rid',objects);
end

objective_loss = evaluate_pmedian_loss(ds,exemplar);

for ii = 1:numel(exemplar_to_check)
    ii
    vertex_candidate = setdiff(ds.rid,exemplar);
    
    for jj = 1:numel(vertex_candidate)
        this_exemplar_list = replace_element(exemplar,vertex_candidate(jj),exemplar_to_check(ii));
        candidate_loss = evaluate_pmedian_loss(ds,this_exemplar_list);
        
        if candidate_loss < objective_loss
            exemplar_to_check(ii)
            vertex_candidate(jj)
            disp('Not the best!')
        end
    
    end
    
end


end

function newlist = replace_element(list,to_add,to_remove)
    newlist = list;
    idx = strcmp(list,to_remove);
    newlist(idx) = to_add;
end