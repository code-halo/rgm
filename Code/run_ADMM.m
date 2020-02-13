function [P,Pp]=get_P(X,Y,P,Z,params)
	N = 2;
	narginchk(5,5);
	if nargin<5
	    fprintf(1,'not enough arguments. Exiting');
	    return
	end
	if isfield(params,'rho')
	    rho = params.rho;
	else % set default value
	    rho = 1;
	end
	if isfield(params,'tau')
	    tau = params.tau;
	else % set default value
	    tau = 1;
	end
	if isfield(params,'mu')
	    mu = params.mu;
	else % set default value
	    mu = 1e-2;
	end
	if isfield(params,'tol')
	    tol = params.tol;
	else % set default value
	    tol=2e-2;
	end
	if isfield(params,'fileID')
	    fileID = params.fileID;
	else % set default value
	    fileID=1;
	end
	if isfield(params,'maxIter')
	    maxIter = params.maxIter;
	else % set default value
	    maxIter=50;
	end
	
	n=size(X,1);
	tol2=9e-9;
	forig = @(P) sum(sum(sqrt((X*P).^2 + (P*Y).^2)));
	
	%P=ones(n)/n;
	U=zeros(n);
	V=zeros(n);
	iter=0;
	va=1;
	corr=lapjv(-P,0.01);
	P1=perm2mat(corr);
	normAPPB=norm(P1*X+Z*P1-Y*P1,'fro');
	normAPPB1=normAPPB;
	fori=1;
	
	Pmin=P;
	fmin=1e10;
	
	% Main loop
	while (normAPPB>tol) && (normAPPB1>tol) && (va > tol2) && (iter < maxIter)
	    iter=iter+1;
	    normAPPB_o=normAPPB;
	    Pold=P;
	    fori_o=fori;
	    
	    [alpha,beta] = admom_sub_solve_for_alpha_beta(X,Y,P,U,V,rho);
	    param.tole = 1e-4;
	    P = linearized_P(X,Y,alpha,beta,P,Z,U,V,rho,tau);
	    U = U + alpha - P*X;
	    V = V + beta - Y*P;
	    fori = forig(P);
	  
	    % Each N iterations, compute the closest permutation matrix and check
	    % if it gives a perfect match
	    if (mod(iter,N)==0)
	        corr=lapjv(-P,0.01);
	        P1=perm2mat(corr);
	        normAPPB1=norm(P1*X+Z*P1-Y*P1,'fro');
	        fau = forig(P1);
	        if fau < fmin
	            fmin=fau;
	            Pmin=P1;
	        end
	    end
	    
	    forigi(iter)=fori;
	
	    % if the optimization gets stuck, reinitialize with the identity matrix
	    if (iter > 5)
	        normAPPB=norm(P1*X+Z*P1-Y*P1,'fro');
	        va = abs(fori - fori_o) + abs(normAPPB- normAPPB_o);
	        va = va/fori;
	        
	        if (iter==6) && (va <1e-7)
	            fprintf('Switching to P = EYE \n');
	            var=1;
	            P=eye(size(P,1));
	        end
	    end;
	    
	    fprintf(fileID, 'Iteration :%i  normAPPB: %1.7f  normAPPB1: %1.7f   forig: %1.7f  va: %1.10f  P-Pold: %1.8f \n',iter,normAPPB,normAPPB1,fori,va,norm(P-Pold,'fro'));
    end;
	    
	corr=lapjv(-P,0.01);
	Pp=perm2mat(corr);
	fau = forig(Pp);
	if fmin < fau
	    fau=fmin;
	    Pp=Pmin;
	end
end
