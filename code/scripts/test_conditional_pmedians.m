num_big = 100;
num_small = 20;

%%Create the data set
mat = [[-3,3]
    mvnrnd([-3,3],.5*eye(2),num_big-1)
    mvnrnd([-3,-3],.5*eye(2),num_small) 
    mvnrnd([3,3],.5*eye(2),num_small) 
    mvnrnd([3,-3],.5*eye(2),num_small) ];

labels = gen_labels(num_big + 3*num_small);
pd = squareform(pdist(mat));
ds = mkgctstruct(pd,'rid',labels,'cid',labels);

%%Run unconditional pmedians
known1 = '';
ex1 = compute_pmedian(ds,4);
sites1 = cell2mat(ex1);

GROUP_VAR = zeros(1,160);
GROUP_VAR(ismember(ds.rid,ex1)) = 1;
figure;
gscatter(mat(:,1),mat(:,2),GROUP_VAR)
title('Unconditional pmedians')
grid on

%%Run conditional pmedians
known2 = {'001'};
ex2 = compute_pmedian(ds,4,'known_ex',known2);
sites2 = cell2mat(ex2);

GROUP_VAR = zeros(1,160);
GROUP_VAR(ismember(ds.rid,known2)) = 1;
GROUP_VAR(ismember(ds.rid,ex2)) = 2;
figure;
gscatter(mat(:,1),mat(:,2),GROUP_VAR)
title('Conditional pmedians')
grid on