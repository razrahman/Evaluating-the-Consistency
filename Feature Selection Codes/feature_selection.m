function FF =feature_selection(CCLE_GENE,CCLE_AUC,GDSC_GENE,GDSC_AUC,NN)
% feature_selection function performs RELIEFF filter feature selection,
% Genetic Algorithm wrapper feature selection approach and LASSO Embedded
% feature selection. Be careful with the number of samples as, it need to
% same for all data provided (CCLE_GENE,CCLE_AUC,GDSC_GENE,GDSC_AUC). Also,
% for Genetic Algorithm wrapper feature selection, it takes a lot of time
% to calculate, so unioin of top 500 features of CCLE and GDSC is taken for 
% performing feature selection.
%     Input:
%         CCLE_GENE: Gene Expression or any other genomic characterization of CCLE database
%         CCLE_AUC: Drug sensitivity measure (AUC, IC_50 etc) of CCLE database
%         GDSC_GENE: Gene Expression or any other genomic characterization of GDSC database
%         GDSC_AUC: Drug sensitivity measure (AUC, IC_50 etc) of GDSC database
%         NN: Number of top features among which common features of CCLE and GDSC is looked upon
%     Output:
%         FF: a vector of 3 which gives number of common features of CCLE
%             and GDSC databases using RELIEFF, Genetic Algorithm and LASSO feature selection.
%
% %   Example:
%        CCLE_GENE = randn(100,1000);
%        CCLE_AUC  = .25.*rand(100,1);
%        GDSC_GENE = randn(100,1000);
%        GDSC_AUC  = .25.*rand(100,1);
%        NN=100;
%        FF =feature_selection(CCLE_GENE,CCLE_AUC,GDSC_GENE,GDSC_AUC,NN);
%% Relieff
[RANK_Relieff_CCLE,~] = relieff(CCLE_GENE,CCLE_AUC,10);
[RANK_Relieff_GDSC,~] = relieff(GDSC_GENE,GDSC_AUC,10);
CR = intersect(RANK_Relieff_CCLE(1:NN),RANK_Relieff_GDSC(1:NN));
%% LASSO
B1 = lasso(CCLE_GENE,CCLE_AUC,'CV',10);
[~,Ind_CCLE_LASSO] = sort(B1(:,1),'descend');
B2 = lasso(GDSC_GENE,GDSC_AUC,'CV',10);
[~,Ind_GDSC_LASSO] = sort(B2(:,1),'descend');
CL = intersect(Ind_CCLE_LASSO(1:NN),Ind_GDSC_LASSO(1:NN));
%% Genetic Algorithm wrapper feature selection
NN2=500;
RU=union(RANK_Relieff_CCLE(1:NN2),RANK_Relieff_GDSC(1:NN2));
CCLE_GENE2=CCLE_GENE(:,RU);
GDSC_GENE2=GDSC_GENE(:,RU);
nFeats = NN;     % Number of expected features
nGen = 100;      % Number of Iterations for GA, more is better
nSurv = 50;  
nPop = 500; 
k = 5;           %fold
mu = 0.02;       % mutation rate
nE = 20;         % number of best features
T = 25;          % number of trees

picked_CCLE_GA = GAfs( CCLE_GENE2, CCLE_AUC, nFeats, nGen, nSurv, nPop, k, mu, nE, T);
picked_GDSC_GA = GAfs( GDSC_GENE2, GDSC_AUC, nFeats, nGen, nSurv, nPop, k, mu, nE, T);
CGA = intersect(picked_CCLE_GA(1,:),picked_GDSC_GA(1,:));

FF=[length(CR) length(CL) length(CGA)];