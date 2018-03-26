function [Is, err] = GAfs( X, y, nFeats, nGen, nSurv, nPop, k, mu, nE, T)
%GAFS Genetic Algorithm based feature selection. Will iteratively try to
% find the top nFeats that minimize the k-fold cross validation error
% for the Random Forest model.
%
%Input: X - MxN matrix, M=Number of samples, N=Number of features.
%       Y - Mx1 vector of responses. Currently Regression only.
%       nFeats - Number of features we are picking, nFeats < N.
%       nGen - Number of generations (iterations) the GA will run.
%       nSurv - Number of chromosomes who survive every generation should 
%           be > nE.
%       nPop - Population of each generation from crossover and mutation. 
%       mu - Mutation Rate, typical be very small, around 0.02.
%       co - Crossover Rate, typically much bigger than mutation.
%       nE - Number of top performing members to keep each generation.
%       T - Number of trees used in RF for generating CV error.
    [M,N] = size(X);
    
    % Start with randomly picked features for each sample.
    [~, out] = sort(rand(nPop,N),2);
    I = unique(sort(out(:,1:nFeats),2), 'rows');
    
    cvErr = zeros(nPop,1);
    cvInd = crossvalind('Kfold', M, k);
    
    % Estimate CV error of our starting set.
    parfor ne = 1:nPop
        cvErr(ne) = RFcvMAE(X(:,I(ne,:)), y, cvInd, T);
    end
    
    % Pick chromosomes with lowest error as our elite
    [err, ind] = sort(cvErr);
    Is = I(ind(1:min(nSurv,size(I,1))),:);
    
    Elite = Is(1:nE,:);
    h = waitbar(0, 'Running GA Feature Selection');
    
    % Main GA loop
    for n1 = 1:nGen
        waitbar(n1/nGen,h);
        I = unique(sort(GACrossFS(Is, nPop, mu, N, Elite),2), 'rows');
        cvErr = zeros(size(I,1),1);
         % Estimate CV error of our starting set.
        parfor ne = 1:size(I,1)
            cvErr(ne) = RFcvMAE(X(:,I(ne,:)), y, cvInd, T);
        end

        % Pick chromosomes with lowest error as our elite
        [err, ind] = sort(cvErr);
        Is = I(ind(1:min(nSurv,size(I,1))),:);

        Elite = Is(1:nE,:);
    end
    
    close(h);
end


function cvErr = RFcvMAE(X, y, cvInd, T)
    k = max(cvInd);
    M = length(y);
    yPred = zeros(M,1);
    
    for n1 = 1:k
        trI = false(M,1);
        tsI = false(M,1);
        tsI(cvInd == n1) = 1;
        trI(cvInd ~= n1) = 1;
        
        Xtr = X(trI,:);
        ytr = y(trI);
        
        RF = TreeBagger(T, Xtr, ytr, 'method', 'regression');
        
        Xts = X(tsI,:);
        
        yPred(tsI) = predict(RF, Xts);
    end
    
    % Calculate MAE
    cvErr = mean(abs((yPred - y)));
    


end

function O = GACrossFS(I, nPop, mu, nF, Elite)
%GENETICALGORITHMTARGETS Performs a single genetic algorithm iteration on
% binary input data.
%
%Input: I - MxN matrix of chromosomes from previous generation.
%                      M=Number of chromosomes, N=Length of chromosome.
%       nPop - Number of offspring in each generation. Should be > M
%       mu - Mutation Rate, should be very small, around 0.02.
%       cRate - Crossover Rate, typically much bigger than mutation.
%       nF - Total number of features.
%       Elite - Top performing chromosomes to keep.

%% Save Elite population
    [M, N] = size(I);
    nElite = size(Elite,1);
    
    O = zeros(nPop, N);
    O(1:nElite,:) = Elite;
    nCross = nPop - nElite;
    
    %% Perform Crossover and mutation on input chromosomes
    % For each child randomly pick a mother and a father
    fatherInd = randi(M, nCross,1);
    motherInd = randi(M, nCross,1);
    
    F = 1:nF;
    
    for n1 = 1:nCross
        % Take the union of the mother and father chromosome
        U =  union(I(fatherInd(n1),:), I(motherInd(n1),:));
   
        % Randomly pick from that union.
        O(nElite+n1,:) = U(randperm(length(U),N));
        
        %% Mutate random chromosome. 
        % Randomly pick gene to mutate.
        Imu = rand(1,N) < mu;
        if any(Imu)
            % Get features contained in mother or father.
            D = setdiff(F, U);
            % Shuffle for randomness
            D = D(randperm(length(D)));
            O(nElite+n1,Imu) = D(1:sum(Imu));
        end
    end
    
    
end

