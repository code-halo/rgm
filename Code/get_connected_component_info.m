function [clustermembership,F_max,numclu_value,opt_threshold] = get_connected_component_info(Z,min_size)

%Here Z is the original matrix , C is the component ids, ic is component
%membership

%% Perform the self-tuned method
    mink = 1;
    threshold = 0.1:0.1:1;
    breaks = length(threshold);
    clusterinfo = [];
    hierarchy = 0;
    entropyfit = zeros(breaks,1);
    balancefit = zeros(breaks,1);
    numclu_cos = zeros(breaks,1);
    F_measure = zeros(breaks,1);
    cosinedistance  = pdist2(Z,Z,'jaccard');
    cosinedistance(isnan(cosinedistance))=1;
    clusterinfo_cell = cell(breaks,3);
    for i=1:length(threshold)
        str=['Working for threshold ',num2str(threshold(i))];
        disp(str);
        temptarget = cosinedistance;
        %% Find indexes which have more similarity
        pertemptarget = sum(temptarget<threshold(i),2);
        hierarchy = hierarchy+1;
        count=1;
        entropy = 0;

        %% Recursively scan the affinity matrix for clusters
        while ~isempty(temptarget)
            [~,index] = max(pertemptarget);
            t = find(temptarget(index,:)<threshold(i));
            if (isempty(t))
                break;
            end
            clusterinfo = [clusterinfo;[count length(t) hierarchy]];
            clusterinfo_cell{i,1} = [clusterinfo_cell{i,1},count];
            clusterinfo_cell{i,2}{count,1} = t;
            clusterinfo_cell{i,3} = hierarchy;
            count = count+1;
            temptarget(t,:) = [];
            temptarget(:,t) = [];
            pertemptarget=sum(temptarget<threshold(i),2);
        end;
        tempclusterinfo = clusterinfo(clusterinfo(:,3)==hierarchy,:);

        %% Remove clusters of small size or outliers
        tempclusterinfo(tempclusterinfo(:,2)<min_size,:) = [];
        clusterinfo(clusterinfo(:,2)<min_size,:) = [];
        vsize = size(Z,1);

        %% Using entropy to identify the ideal number of clusters
        %% Convert count to probability
        numclu_cos(i) = length(find(clusterinfo(:,3)==i));
        if (numclu_cos(i)>=mink)
            %% Calculate Expected Balance
            for j=1:size(tempclusterinfo,1)
                balancefit(i,1) = balancefit(i,1) + (1.0*tempclusterinfo(j,2)/vsize);
            end;
            tempclusterinfo(:,2) = tempclusterinfo(:,2);
            %% Obtain entropy for each clustering
            for j=1:size(tempclusterinfo,1)
                entropy =  entropy  + (-1)*(tempclusterinfo(j,2)/vsize)*log((tempclusterinfo(j,2)/vsize));
            end;
            entropyfit(i,1) = entropy;
            balancefit(i,1) = balancefit(i,1)/size(tempclusterinfo,1);
            F_measure(i,1) = 2*(entropyfit(i,1)*balancefit(i,1))/(balancefit(i,1)+entropyfit(i,1));
        else
            F_measure(i,1) = 0.0;
        end;
        clear pertemptarget;
    end;
    [F_max,index] = max(F_measure);
    opt_threshold = threshold(index);
    fprintf('Max F-measure = %f and threshold %f\n',F_max,opt_threshold);
    numclu_value = length(find(clusterinfo(:,3)==index));

    %% Plot for the F-measure
%     figure;
%     plot(threshold,F_measure,'r*');
%     xlabel('Threshold t');
%     ylabel('F-measure');
%     title('Plot of threshold vs F-measure value for Self-Tuned Version');

    %% Get the info w.r.t to Z matrix the connected components
    clustermembership_info = clusterinfo_cell{index,2};
    clustermembership = zeros(size(Z,1),1);
    count=1;
    for i=1:length(clustermembership_info)
        component_info=clustermembership_info{i};
        indices = find(clustermembership==0);
        if(length(component_info)>=min_size)
            clustermembership(indices(component_info))=count;
            count=count+1;
        else
            for j=1:length(component_info)
                clustermembership(indices(component_info(j)))=count;
                count=count+1;
            end
        end
    end
    zero_indices = find(clustermembership==0);
    for i=1:length(zero_indices)
        clustermembership(zero_indices(i))=count;
        count=count+1;
    end
end

