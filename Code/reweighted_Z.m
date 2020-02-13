function Z = Reweighted_Z(alpha,beta,P,Z,C,D,nu,rho)
    O = beta - alpha;
    itr = 0;
    max_itr = 100;
    stp = 1e-3;
    while(itr < max_itr)
        grad =  0.5*(Z*P-O)*P' + nu*(Z+D-C);
        %grad =  0.5*(Z*P-O)*P' + nu*(Z-D-C);
        itr = itr + 1;
        Z = Z - stp*grad;
        stp = stp/2;
        Z(Z<0) = 0;
        Z(Z>1) = 1;
        Z = (Z+Z')/2;
    end
end
