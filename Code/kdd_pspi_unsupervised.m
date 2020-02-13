function [] = kdd_pspi_unsupervised(inputset,iterations,nu,mu)
    warning('off');
    params.rho = 1; %
    params.nu = nu; %Regularizer on Z
    params.mu = mu;
    params.max_itr = 150;       %ADMM maximum iterations
    params.tol1 = 1e-6;
    params.tol2 = 1e-10;
    
    Fdir = ['results/',inputset];
    mkdir(Fdir);
    outputfile = sprintf([Fdir,'/pspi_results_Exp_%d_nu_%f_mu_%f.txt'],iterations,nu,mu);
    fileID = fopen(outputfile,'w');
    [X,Y,Zorig,~] = get_data(inputset,iterations);
    Zorig(Zorig <0) = 0;
    n = size(X,1);
    %Z = Y;
    Z = Y-X;
    Z(Z<0) = 0;
    P = ones(n)/n;
    tauA = max(eigs(X*X'));
    tauB = max(eigs(Y*Y'));
    tau = max(tauA,tauB);
    params.tau = 0.8/tau;
    params.fileID = fileID;
    
    [P, Z, ~, Obj2] = kdd_run_ADMM(X,Y,P,Z,params);
    corr = lapjv(-P,0.01);
    Pdisc = perm2mat(corr);
    Z_rmse = sqrt(sum(sum((Zorig - Z).^2))/(n*n));
    
    
    Z(Z<1)=0;
    Z(Z>=1) = 1;
    sumZ = sum(Z(:));
    if (regexp(inputset,'Exhaustive'))
        min_size = 5;
    else
        min_size = 10;
    end
    fprintf(fileID,'Num of non-zero entries in Z: %d\n',sumZ);
    fprintf(fileID,'Matching Error: %f\n',Obj2);
    fprintf(fileID,'Z RMSE Error: %f\n',Z_rmse);
    if sumZ == 0
        return
    end
    outP_id = sprintf([Fdir,'/P_Exp_%d_nu_%f_mu_%f.txt'],iterations,nu,mu);
    csvwrite(outP_id,Pdisc);
    outZ_id = sprintf([Fdir,'/Z_Exp_%d_nu_%f_mu_%f.txt'],iterations,nu,mu);
    csvwrite(outZ_id,Z);
    
    [Zbins,Fmax_I,opt_k,~] = get_connected_component_info(Z,min_size);
    [ZObins,Fmax_O,true_k,~] = get_connected_component_info(Zorig,min_size);
    fprintf(fileID,'Fmax Inferred: %d for Opt k: %d\n',Fmax_I,opt_k);
    fprintf(fileID,'Fmax Groundtruth: %d for True k: %d\n',Fmax_O,true_k);
    [C1,~,ic1] = unique(Zbins);
    [C2,~,ic2] = unique(ZObins);
    [~,id1] = sort(Zbins);
    [~,id2] = sort(ZObins);
    
    th = min_size;
    cmp1 = C1(accumarray(ic1,1)>=th); %connected components number which appears more than th times
    fprintf(fileID,'No of Connected Components in Infered Z with size %d is : %d\n',th, length(cmp1));
    cmp2 = C2(accumarray(ic2,1)>=th); %connected components number in the original graph which appears more than the times
    fprintf(fileID,'No of Connected Components in Original Z with size %d is : %d\n',th, length(cmp2));
    if isempty(cmp2) || isempty(cmp1)
        fclose(fileID);
        return
    end
    Zc = zeros(size(Z,1),size(Z,2));
    for c=1:length(cmp1)
        c1 = cmp1(c);
        indices = find(Zbins==c1);
        Zc(indices,indices) = Z(indices,indices);
    end
    %imshow(Zorig(id2,id2));
    %imshow(Zc(id1,id1));
    
    indxo = find(Zorig>0);
    indxi = find(Zc>0);
    tpz = intersect(indxo,indxi);
    prec = length(tpz)/length(indxi);
    rec = length(tpz)/length(indxo);
    fprintf(fileID,'Connected component with size %d. Precision : %f Recall: %f \n',th, prec, rec);
    fclose(fileID);
    
end

