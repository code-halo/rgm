function G = project_EYE(P)
    %we implement the BBS algorithm using euclidean projection rather than the KL-divergence projection
    tol=1e-6;
    n=size(P,1);
    err = 1;
    itr = 0;
    max_itr = 25;
    G = P;
    one = ones(n,n);
    I = eye(n,n);
    while (err > tol) && (itr < max_itr)
        itr = itr + 1;
        G_old = G;
        G = G + ((I - G + (one*G)/n)*one)/n - (one*G)/n;
        G = max(0,G);
        err = norm(G-G_old,1);
    end
end
