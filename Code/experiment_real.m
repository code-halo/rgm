function [] = experiment_real(data)
    numItermax = 50;
    for i=-5:0
        nu = 2^i;
        for j=-5:0
            mu = 2^j;
            kdd_pspi_unsupervised_real(data,nu,mu);
        end
    end
end

