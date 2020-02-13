function C = gl_solve_for_C(Z,D,rho,nu,mu)
    %mu = 0.01;
    W = rho*sign(Z-D).*max(abs(Z-D)-mu/rho,0);
    %W = rho*sign(Z-D).*max(abs(Z-D)-nu/rho,0);
    c_norm = sqrt(sum(abs(W).^2,1));
    V = max((c_norm-nu)./(rho*c_norm),0);
    dims = size(W);
    V1 = repmat(V,dims(1),1);
    C = W.*V1;
end
