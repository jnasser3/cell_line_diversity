function objective = evaluate_pmedian_loss(ds,exemplar,type)
%objective = evaluate_pmedian_loss(ds,exemplars)
%
%Calculate the objective function for a pairwise distance or similarity matrix and a set of
%exemplars. If type = 'distance', the objective is the average 
%distance from each point to its nearest exemplar. If type = 'similarity',
%the objective is the average similarity from each point to its
%nearest exemplar.

exemplar_idx = ismember(ds.rid,exemplar);

switch type
    case 'distance'
        objective = sum(min(ds.mat(:,exemplar_idx),[],2));
    case 'similarity'
        objective = sum(max(ds.mat(:,exemplar_idx),[],2));
end

%Average objective for easier interpretation
objective = objective / numel(ds.rid);


end

