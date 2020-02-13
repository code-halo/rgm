function [X,Y,Z,P] = get_real_data(dataset)
    filepath = ['../Datasets/Real/',dataset];

    X = load([filepath,'/','X.txt'],'-ascii');
    Y = load([filepath,'/','Y.txt'],'-ascii');
    Z = load([filepath,'/','Z.txt'],'-ascii');
    P = load([filepath,'/','P.txt'],'-ascii');
end
