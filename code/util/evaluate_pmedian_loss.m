function loss = evaluate_pmedian_loss(ds,exemplar)
%loss = evaluate_pmedian_loss(ds,exemplars)
%
%Calculate the loss function for a pairwise distance matrix and a set of
%exemplars. The loss is equal to the sum of the distances from each
%non-exemplar to it's nearest exemplar.

loss = 0;
exemplar_idx = ismember(ds.rid,exemplar);

for ii = 1:numel(ds.rid);
    loss = loss + min(ds.mat(ii,exemplar_idx));
end


end

