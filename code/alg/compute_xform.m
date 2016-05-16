function [xform,x2] = compute_xform(null0,xform_method,pdf_estimate_method)
%Create the function that maps correlation to contribution.
%Less than the median null map to 1
%Greater than the max null map to 0
%In between follow the shape of the null

switch xform_method
    case 'pdf'
        numpoints = 2^10;
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
        numpoints = 2^10;
        med0 = median(null0);
        max0 = max(null0);
        
        %compute the empirical cdf of the null. We have a soft threshold at
        %the median.
        null0 = null0(null0 >= med0);
        [~,null_pdf,x,f] = kde(null0,numpoints);
        f = clip(f,0,1);
        
%       THIS IS ONLY NEEDED IF USING MATLAB's ecdf.
%         If the first two entries of x are equal interp1 will fail. So just hack it
%         if x(1) == x(2);
%             f(1) = [];
%             x(1) = [];
%         end
        
        %The contribution to diversity should be 1 minus the cdf
        f = 1 - f;
        
        %Compute the transform.
        x2 = linspace(-1,1,numpoints);
        xform = zeros(1,numpoints);
        xform(x2 < med0) = 1;
        
        xform(x2 >= med0 & x2 <= max0) = interp1(x,f,x2(x2 >= med0 & x2 <= max0),'linear');
        xform(x2 > max0) = 0;
        
        %Diagnostic plot
        figure;
        plot(x2,xform,'DisplayName','xform')
        hold on
        plot(x2,null_pdf,'DisplayName','pdf')
        hold on
        plot(x2,f,'DisplayName','cdf')
        xlabel('Correlation')
        ylabel('Transform - contribution to Diversity')
        namefig('Correlation_contribution_transform');
        grid on
        
        
        
end

end