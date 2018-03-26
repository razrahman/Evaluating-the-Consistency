clear; close all; clc;
nFeats = 50;
nGen = 50;
nSurv = 20;
nPop = 100;
k = 5;
mu = 0.02;
nE = 10;
T = 25;

%% Creatue Synthetic Dataset
nsamp = 100;
ncov = 1000;
X=normrnd(0,1,[nsamp,ncov]);
    xxcov1=[ones(nsamp,1),X];
    nstrong=50;
    nweak=50;
    n0=ncov-nstrong-nweak;
    betaint=1.4; % response offset
    betastrong=unifrnd(1.5,3,[nstrong,1]);
    betaweak=unifrnd(0,.5,[nweak,1]);
    beta0=zeros(n0,1);
    betafull= [betastrong;betaweak;beta0];
    betafull = betafull(randperm(length(betafull)));
    [w, Top] = sort(betafull, 'descend');
    betafull = [betaint; betafull];
    Y=xxcov1*betafull;

[ranked,weights] = relieff(X,Y,10);
relieffP = ranked(1:nFeats);

tic
[picked, err] = GAfs( X, Y, nFeats, nGen, nSurv, nPop, k, mu, nE, T);
toc
cvInd = crossvalind('Kfold', nsamp, k);
foldErrFF = RFcvMAE(X(:,ranked(1:nFeats)), Y, cvInd, T);
foldErrGA = RFcvMAE(X(:,picked(1,:)), Y, cvInd, T);