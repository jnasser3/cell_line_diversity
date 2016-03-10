function confirm_pmedian_local_minimum(ds,exemplar,exemplar_to_check)
%Sanity check to make sure we are in a local minimum.
%
%Checks to see if we can get a better optimum by making a single exemplar
%exchange
%
%

objective_loss = evaluate_pmedian_loss(ds,exemplar)

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