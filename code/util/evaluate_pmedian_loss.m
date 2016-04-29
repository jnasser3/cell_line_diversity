function loss = evaluate_pmedian_loss(ds,exemplar)
%loss = evaluate_pmedian_loss(ds,exemplars)
%
%Calculate the loss function for a pairwise distance matrix and a set of
%exemplars. The loss is equal to the sum of the distances from each
%non-exemplar to it's nearest exemplar.

exemplar_idx = ismember(ds.rid,exemplar);
loss = sum(min(ds.mat(:,exemplar_idx),[],2));

%Average loss for easier interpretation
loss = loss / numel(ds.rid);

% %Slow version of the vectorized code
% for ii = 1:numel(ds.rid);
%     loss = loss + min(ds.mat(ii,exemplar_idx));
% end


end

