%Visualize the relationship between choosing diverse vs representative
%lines in the presence of outliers and clusters

n = 100;

%Generate the data
cl1 = mvnrnd([2 6],.025*eye(2),n);
cl2 = mvnrnd([3 5],.025*eye(2),n);
cl3 = mvnrnd([1,1],.01*eye(2),10);
out1 = [1,1];
out2 = [3,8];


data = vertcat(cl1,cl2,cl3,out1,out2);

%Make into a struct
dist = squareform(pdist(data));
labels = gen_labels(size(data,1));
dist_struct = mkgctstruct(dist,'rid',labels,'cid',labels);

%Compute the diverse elements
diverse_elements = compute_maximally_diverse(dist_struct,5)
exemplars = compute_pmedian(dist_struct,5)

%% Plot
figure;
scatter(data(:,1),data(:,2),[],'b','filled','DisplayName','data')
hold on
scatter(data(str2num(cell2mat(diverse_elements)),1),...
    data(str2num(cell2mat(diverse_elements)),2),...
    [],...
    'r',...
    'filled',...
    'DisplayName','Maximally Diverse')
hold on
scatter(data(str2num(cell2mat(exemplars)),1),...
    data(str2num(cell2mat(exemplars)),2),...
    [],...
    'g',...
    'filled',...
    'DisplayName','Most Representative')
grid on
legend('show','location','northwest')
axis([0 5 0 10])
title('Choosing Diverse vs Representative lines in the presence of clusters and outliers')