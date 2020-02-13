%% To run of a015_N_500_m1_50_m2_100 datasets (Scale Free Experiments)
datasets = {'a015_N_500_m1_50_m2_100_noise_50','a015_N_500_m1_50_m2_100_noise_50_both',...
            'a015_N_500_m1_50_m2_100_noise_30','a015_N_500_m1_50_m2_100_noise_30_both'};
for j=1:4
    data = datasets{j};
    for i=1:5
        experiment(data,i)
    end
end

%Results are in results/a015* / folder

%% To run of r015_N_500_m1_50_m2_100 datasets (Random-Geometric Experiments)
datasets = {'r015_N_500_m1_50_m2_100_noise_50','r015_N_500_m1_50_m2_100_noise_50_both',...
            'r015_N_500_m1_50_m2_100_noise_30','r015_N_500_m1_50_m2_100_noise_30_both'};
for j=1:4
    data = datasets{j};
    for i=1:5
        experiment(data,i)
    end
end

%Results are in results/r015* / folder

%% To run of DREAM4 Challenge datasets
for j=1:5
    data = sprintf('Exp_%d_Exhaustive',j);
    for i=1:3
        experiment(data,i)
    end
end

%Results are in results/Exp_j_Exhaustive/ folder

%% To run experiment on C. Elegans network
data = sprintf('CEL');
experiment_real(data)

%Results are in results/CEL/ folder

%% To run experiments on TCGA network
data = sprintf('TCGA');
experiment_real(data)

% Results are in results/TCGA/ folder
