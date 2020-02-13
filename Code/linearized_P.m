function Pcont = Linearized_P(X,Y,alpha,beta,P,Z,U,V,rho,tau)
    M = alpha+U;
    N = beta+V;
    O = beta-alpha;

    F = 0.5*rho*norm(Z*P-O,'fro')^2 + 0.5*rho*norm(M-P*X,'fro')^2 + 0.5*rho*norm(N-Y*P,'fro')^2;
    Fhat = F+1;
    k = 0;
    while (Fhat > F) && (k<10)
        k = k+1;
        R1 = P + rho*tau*(M-P*X)*X/k;
        R2 = P + rho*tau*Y*(N-Y*P)/k;
        R3 = P + rho*tau*Z*(O-Z*P)/k;
        Pcont = project_EYE((R1 + R2 + R3)/3);
        Fhat = 0.5*rho*norm(Z*Pcont-O,'fro')^2 + 0.5*rho*norm(M-Pcont*X,'fro')^2 + 0.5*rho*norm(N-Y*Pcont,'fro')^2;
    end
end
