function [xformcorrs,xform,x2] = xform_corr_diversity(corrs,null0,xform_method,pdf_estimate_method)
%[xformcorrs] = xform_corr_diversity(corrs,null0)
%
%Given a null distribution of correlations and a set of observations,
%computes the strength of the observations relative to the null.

[xform,x2] = compute_xform(null0,xform_method,pdf_estimate_method);
xformcorrs = interp1(x2,xform,corrs);

end

function [xform,x2] = compute_xform(null0,xform_method,pdf_estimate_method)
%Create the function that maps correlation to contribution.
%Less than the median null map to 1
%Greater than the max null map to 0
%In between follow the shape of the null

switch xform_method
    case 'pdf'
        numpoints = 1024;
        switch pdf_estimate_method
            case 'matlab'
                [f, x] = ksdensity(null0,'npoints',numpoints);
            case 'fast_kde'
                [~,f, x] = kde(null0,numpoints);
        end
        f = f / max(f);

        med0 = median(null0);
        max0 = max(null0);

        x2 = linspace(-1,1,numpoints);
        xform = zeros(1,numpoints);
        xform(x2 < med0) = 1;
        xform(x2 >= med0 & x2 < max0) = interp1(x,f,x2(x2 >= med0 & x2 < max0));
        xform(x2 > max0) = 0;
    case 'cdf'        
        numpoints = 1024;
        med0 = median(null0);
        max0 = max(null0);
        
        %compute the empirical cdf of the null. We have a hard threshold at
        %the median.
        null0 = null0(null0 >= med0);
        null0 = round(null0,3);
        [f,x] = ecdf(null0);
        
        %If the first two entries of x are equal interp1 will fail. So just hack it
        if x(1) == x(2);
            f(1) = [];
            x(1) = [];
        end
        
        %The contribution to diversity should be 1 minus the cdf
        f = 1 - f;
        
        %Compute the transform.
        x2 = linspace(-1,1,numpoints);
        xform = zeros(1,numpoints);
        xform(x2 < med0) = 1;
        xform(x2 >= med0 & x2 < max0) = interp1(x,f,x2(x2 >= med0 & x2 < max0));
        xform(x2 > max0) = 0;
end

end