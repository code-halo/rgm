function [P, Z, obj1, obj2] = kdd_run_ADMM(X,Y,P,Z,params)
    forig = @(X,Y,P,Z,nu,mu) sum(sum(sqrt((P*X).^2 + (Y*P).^2))) + 0.5*norm(P*X + Z*P -Y*P,'fro')^2 + nu*max(sum(abs(Z),1)) + mu*sum(sum(sqrt(Z)));
    itr = 0;
    rho = params.rho;
    nu = params.nu;
    tau = params.tau;
    mu = params.mu;
    fileID = params.fileID;
    obj1 = 1;
    obj2 = 1;
    n = size(P,1);
    P_h = zeros(n); Z_h = zeros(n); U = zeros(n); V = zeros(n); D = zeros(n);
    
    while ((obj1 > params.tol1) && (obj2 > params.tol1) && (itr < params.max_itr))
        itr = itr + 1;
        [alpha,beta] = admom_sub_solve_for_alpha_beta(X,Y,P,U,V,rho);
        C = gl_solve_for_C(Z,D,rho,nu,mu);
        P = linearized_P(X,Y,alpha,beta,P,Z,U,V,rho,tau);
        Z = reweighted_Z(alpha,beta,P,Z,C,D,nu,rho);
        U = U + (alpha-P*X);
        V = V + (beta-Y*P);
        D = D + (C-Z);
        if (mod(itr,3) == 0)
            corr = lapjv(-P,0.01);
            Pdisc = perm2mat(corr);
            obj1 = norm(P*X+Z*P-Y*P,'fro'); %matching error for P continous
            obj2 = norm(Pdisc*X+Z*Pdisc-Y*Pdisc,'fro'); %matching error for P discrete
            fori = forig(X,Y,Pdisc,Z,nu,mu);   %Original Objective
            fprintf(fileID,'Iteration:%d Matching Error(cont):%1.6f   Matching Error(disc):%1.6f  ObjFunc:%1.6f  P-diff:%1.6f   Z-diff:%1.6f\n', itr,obj1,obj2,fori,norm(P-P_h,'fro'),norm(Z-Z_h,'fro'));
            P_h = P;
            Z_h = Z;
        end
    end
end
