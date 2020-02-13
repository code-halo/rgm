function [X,Y,Z,P] = get_data(fname,iteration)
    if (regexp(fname,'Exhaustive'))
        filepath = ['../Datasets/DREAM4/',fname];
    else
        filepath = ['../Datasets/',fname];
    end
    X_file_name = ['Edgelist_A_',num2str(iteration),'.csv'];
    Y_file_name = ['Edgelist_B_',num2str(iteration),'.csv'];
    Z_file_name = [filepath,'/Z_',num2str(iteration),'.csv'];
    if (regexp(fname,'Exhaustive'))
        P_file_name = [filepath,'/P_',num2str(iteration),'.csv'];
    else
        P_file_name = [filepath,'/Permutation_',num2str(iteration),'.csv'];
    end
    X_edgelist = load([filepath,'/',X_file_name],'-ascii');
    Y_edgelist = load([filepath,'/',Y_file_name],'-ascii');
    Z = load(Z_file_name,'-ascii');
    P = load(P_file_name,'-ascii');

    net_X = [X_edgelist(:,1) X_edgelist(:,2) X_edgelist(:,3)];
    net_X = [net_X; [X_edgelist(:,2) X_edgelist(:,1) X_edgelist(:,3)]];
    net_X = unique(net_X,'rows');
    X = full(spconvert(net_X));

    net_Y = [Y_edgelist(:,1) Y_edgelist(:,2) Y_edgelist(:,3)];
    net_Y = [net_Y; [Y_edgelist(:,2) Y_edgelist(:,1) Y_edgelist(:,3)]];
    net_Y = unique(net_Y,'rows');
    Y = full(spconvert(net_Y));
end

