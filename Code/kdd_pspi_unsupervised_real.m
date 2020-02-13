function [] = kdd_pspi_unsupervised_real(inputset,nu,mu)
    warning('off');
    params.rho = 1; %
    params.nu = nu; %Regularizer on Z
    params.mu = mu;
    params.max_itr = 150;       %ADMM maximum iterations
    params.tol1 = 1e-6;
    params.tol2 = 1e-10;
    
    Fdir = ['results/',inputset];
    mkdir(Fdir);
    outputfile = sprintf([Fdir,'/pspi_results_nu_%f_mu_%f.txt'],nu,mu);
    fileID = fopen(outputfile,'w');
    [X,Y,Zorig,~] = get_real_data(inputset);
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
    min_size = 10;
    fprintf(fileID,'Num of non-zero entries in Z: %d\n',sumZ);
    fprintf(fileID,'Matching Error: %f\n',Obj2);
    fprintf(fileID,'Z RMSE Error: %f\n',Z_rmse);
    if sumZ == 0
        return
    end
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

% %% Define axes
% xticks = unique([1:25:51,51:50:202]);
% xticklabels = unique([1:25:51,51:50:202]);
% yticks = unique([1:25:51,51:50:202]);
% yticklabels = unique([1:25:51,51:50:202]);
% 
% %% Additional Analysis for CEL network
% X_true = load('../Datasets/Real/CEL/true_X.txt','-ascii');
% h = figure;
% subplot(2,3,1);
% imagesc(X_true);
% set(gca,'XTick',xticks,'XTickLabel',xticklabels);
% set(gca,'YTick',yticks,'YTickLabel',yticklabels);
% title('Original X');
% 
% subplot(2,3,2);
% imagesc(Y);
% set(gca,'XTick',xticks,'XTickLabel',xticklabels);
% set(gca,'YTick',yticks,'YTickLabel',yticklabels);
% title('Original Y');
% 
% subplot(2,3,3);
% imagesc(Z);
% set(gca,'XTick',xticks,'XTickLabel',xticklabels);
% set(gca,'YTick',yticks,'YTickLabel',yticklabels);
% title('Original Z');
% 
% subplot(2,3,4)
% imagesc(X);
% xticklabels = unique([i(1:25:51);i(51:50:202)]);
% yticklabels = unique([i(1:25:51);i(51:50:202)]);
% set(gca,'XTick',xticks,'XTickLabel',xticklabels);
% set(gca,'YTick',yticks,'YTickLabel',yticklabels);
% title('Permuted X');
% 
% subplot(2,3,5);
% imagesc(Zorig(id2,id2));
% xticklabels = unique([id2(1:25:51);id2(51:50:202)]);
% yticklabels = unique([id2(1:25:51);id2(51:50:202)]);
% set(gca,'XTick',xticks,'XTickLabel',xticklabels);
% set(gca,'YTick',yticks,'YTickLabel',yticklabels);
% title('Original Structured Z');
% 
% subplot(2,3,6);
% imagesc(Zc(id1,id1));
% xticklabels = unique([id1(1:25:51);id1(51:50:202)]);
% yticklabels = unique([id1(1:25:51);id1(51:50:202)]);
% set(gca,'XTick',xticks,'XTickLabel',xticklabels);
% set(gca,'YTick',yticks,'YTickLabel',yticklabels);
% title('Inferred Structured Z');
% 
% 
% %% Additional Analysis for TCGA network
% h = figure;
% subplot(2,3,1);
% imshow(X);
% title('Original X');
% subplot(2,3,2);
% imshow(Y);
% title('Original Y');
% subplot(2,3,3);
% imshow(Zorig);
% title('Original Z');
% subplot(2,3,4);
% imshow(Zorig(id2,id2));
% title('Original Structured Z');
% subplot(2,3,5);
% imshow(Zc(id1,id1));
% title('Inferred Structured Z');
% 
% node_cluster_info = [find(Zbins==cmp1(1)),ones(length(find(Zbins==cmp1(1))),1);find(Zbins==cmp1(2)),2*ones(length(find(Zbins==cmp1(2))),1);find(Zbins==cmp1(3)),3*ones(length(find(Zbins==cmp1(3))),1)];
% csvwrite('results/TCGA/output_node_cluster.csv',node_cluster_info);


