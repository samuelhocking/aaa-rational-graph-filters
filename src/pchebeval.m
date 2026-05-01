function yp = pchebeval(cc, x)
    x = x/2;
    vprev = 1;                                      % evaluate p(L) by chebyshev recurrence
    vcurr = x;
    yp = cc(1)*vprev + cc(2)*vcurr;
    for j = 3:length(cc)
        vnext = 2*x*vcurr - vprev;
        yp = yp + cc(j)*vnext;
        vprev = vcurr;
        vcurr = vnext;
    end
end