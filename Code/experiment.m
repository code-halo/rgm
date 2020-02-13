function [] = experiment(data,iter)
    numItermax = 50;
    for i=-5:0
        nu = 2^i;
        for j=-5:0
            mu = 2^j;
            kdd_pspi_unsupervised(data,iter,nu,mu);
        end
    end
end


