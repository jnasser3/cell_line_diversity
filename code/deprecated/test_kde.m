
num_sample = 1e6;
mesh_size = 2 ^ 14;

sample = [randn(100,1);randn(100,1)*2+35 ;randn(100,1)+55];
%sample(sample < 0) = 0;

tic;[bw_kde,f2,x2] = kde(sample, mesh_size);toc
tic;[f,x,bw] = ksdensity(sample,'npoints',mesh_size,'bandwidth',bw_kde);toc;


bw
bw2

plot(x,f,'DisplayName','ksdensity')
hold on
plot(x2,f2,'DisplayName','kde')
legend('show')
grid on