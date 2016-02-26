function xformbioa = xform_bioa_diversity(bioa,bkg)
%[xformbioa] = xform_bioa_diversity(bioa,bkg)
%
%Given a background of cc_q75's, computes the strength of the observed
%values relative to the background

[xform,x2] = compute_bioa_xform(bkg);
xformbioa = interp1(x2,xform,bioa);

end

function [xform,x2] = compute_bioa_xform(bkg)
%Create the function that maps bioa to contribution.
%Less than the minimum maps to 0
%Greater than the max maps to 1
%In between follow the cdf of the bkg

numpoints = 201;

[f,x] = ecdf(bkg);

%If the first two entries of x are equal interp1 will fail. So just hack it
if x(1) == x(2);
    f(1) = [];
    x(1) = [];
end

minbkg = min(bkg);
maxbkg = max(bkg);

x2 = linspace(-1,1,numpoints);
xform = zeros(1,numpoints);

xform(x2 < minbkg) = 0;
xform(x2 > minbkg & x2 < maxbkg) = interp1(x, f, x2(x2 > minbkg & x2 < maxbkg));
xform(x2 > maxbkg) = 1;
end