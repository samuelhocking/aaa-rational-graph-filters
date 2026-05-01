function yp = applypLx(x,L,cheb_poly_coeffs)
    arguments
        x
        L
        cheb_poly_coeffs
    end

    vprev = x;
    vcurr = L*x;
    yp = cheb_poly_coeffs(1)*vprev + cheb_poly_coeffs(2)*vcurr;
    for j = 3:length(cheb_poly_coeffs)
        vnext = 2*L*vcurr - vprev;
        yp = yp + cheb_poly_coeffs(j)*vnext;
        vprev = vcurr;
        vcurr = vnext;
    end
end