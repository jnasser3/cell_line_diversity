function sim_mat = evaluate_consensus_clustering(mat,consensus_label)
%
%Input
%       mat: Clustering matrix produced by run_consensus_clustering
%       consensus_label: A string denoting the feature that varies among the
%           clustering. Eg 'Sample Holdout'.
%
%Output:
%       sim_mat: Similarity matrix. A square matrix denoting the frequency
%           of co-clustering. sim_mat(i,j) = x means that elements i and j
%           clustered together x percent of the time. 

[n,~] = size(mat);

%%Create consensus similarity matrix
sim_mat = zeros(n,n);

for ii = 1:n
    for jj = ii:n
        
        %subset to relevant part of consensus matrix
        temp = mat([ii jj], :);
        
        %Remove zeros
        temp1 = prod(temp);
        temp(:,temp1 == 0) = [];
        
        %Compute similarity
        if length(temp(1,:)) > 0
            sim_mat(ii,jj) = nnz(~(temp(1,:) - temp(2,:))) / length(temp(1,:));
        else
            sim_mat(ii,jj) = 1;
        end
        
        
    end
end

%Fill in lower triangle
sim_mat = sim_mat + sim_mat';

%Set diagonal to 1. The previous symmetrization step set diagonal to 2. 
sim_mat(logical(eye(size(sim_mat)))) = 1;

%% Plots
Z = linkage(1 - sim_mat);
idx = optimalleaforder(Z,1-sim_mat);

figure
imagesc(sim_mat(idx,idx))
colormap(flipud(gray))
title_str = sprintf(['Consensus Clustering - %s\n',...
    'Num Exemplar = %d\n',...
    'Num clusterings = %d'],...
    consensus_label,nnz(unique(mat(:,1)) > 0),size(mat,2));
title(title_str)

% ax = gca;
% 
% ax.XTick =1:length(row_labels);
% ax.XTickLabels = row_labels(idx);
% ax.FontSize = 1;
% rotateXLabels(ax,90)

end

