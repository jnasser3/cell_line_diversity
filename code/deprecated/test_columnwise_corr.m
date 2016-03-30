%% Test script to make sure vectorized columnwise corr functions implemented correctly

p = 978;
n = 7000;

A = randn(p,n);
B = randn(p,n);

%% Pearson
%compute the slow way
corr_pearson = zeros(1,n);
tic;
for ii = 1:n
    corr_pearson(ii) = corr(A(:,ii),B(:,ii));
end
toc;

%vectorized
tic;
ZA = zscore(A);
ZB = zscore(B);
corr_pearson_vec = sum(ZA .* ZB) / (p-1);
toc;

%mortar
tic;
corr_mortar = pwcorr(A,B);
toc;

%norm(corr_pearson - corr_pearson_vec)
norm(corr_mortar' - corr_pearson)

%% Spearman
%compute the slow way
% corr_spearman = zeros(1,n);
% tic;
% for ii = 1:n
%     corr_spearman(ii) = corr(A(:,ii),B(:,ii),'type','spearman');
% end
% toc;
% 
% %vectorized
% tic;
% ZA = zscore(rankorder(A));
% ZB = zscore(rankorder(B));
% corr_spearman_vec = sum(ZA .* ZB) / (p-1);
% toc;
% 
% norm(corr_spearman - corr_spearman_vec)
