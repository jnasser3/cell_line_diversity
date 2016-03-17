



%Generate the data
cl1 = mvnrnd([3,3],.1*eye(2),100);
cl2 = mvnrnd([5,5],.1*eye(2),100);
out1 = [1,1];
out2 = [7,7];

data = vertcat(cl1,cl2,out1,out2);

%Make into a struct
dist = squareform(pdist(data));
labels = gen_labels(size(data,1));
dist_struct = mkgctstruct(dist,'rid',labels,'cid',labels);

%Compute the diverse elements
diverse_elements = compute_maximally_diverse(dist_struct,2);
exemplars = compute_pmedian(dist_struct,2);

%% Plot
figure;
scatter(data(str2num(cell2mat(diverse_elements)),1),...
    data(str2num(cell2mat(diverse_elements)),2),...
    100,...
    'r',...
    'filled',...
    'DisplayName','Maximally Diverse')
hold on
scatter(data(:,1),data(:,2),[],'b','filled','DisplayName','data')
hold on
scatter(data(str2num(cell2mat(exemplars)),1),...
    data(str2num(cell2mat(exemplars)),2),...
    100,...
    'g',...
    'filled',...
    'DisplayName','Most Representative')
grid on
legend('show','location','northwest')
axis([0 8 0 8])
title('Choosing Diverse vs Representative lines in the presence of clusters and outliers')