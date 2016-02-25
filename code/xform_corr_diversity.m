function [xformcorrs] = xform_corr_diversity(corrs,null0)
%[xformcorrs] = xform_corr_diversity(corrs,null0)
%
%Given a null distribution and a set of observations. Computes the rank of each
%observation with respect to the null in a transformed way. 

[xform,x2] = compute_xform(null0);
xformcorrs = interp1(x2,xform,corrs);

end

function [xform,x2] = compute_xform(null0)
%Create the function that maps correlation to contribution.
%Less than the median null map to 1
%Greater than the max null map to 0
%In between follow the pdf of the null

numpoints = 1000;
[f, x] = ksdensity(null0,'npoints',numpoints);
f = f / max(f);

med0 = median(null0);
max0 = max(null0);

x2 = linspace(-1,1,numpoints);
xform = zeros(1,numpoints);
xform(x2 < med0) = 1;
xform(x2 >= med0 & x2 < max0) = interp1(x,f,x2(x2 >= med0 & x2 < max0));
xform(x2 > max0) = 0;
end