% This code will generate figure 10 of the paper. Before running it,
% download the following 3 files and keep it in the same folder with main code
% (1) CCLE_Drug_Sensitivity_Raziur_Extended.xlsx from
%      https://drive.google.com/file/d/0B51mihVyZFnVNzdvRWR4SVFodFk/view?usp=sharing
% (2) GDSC_v6_v17_fitted_dose_response.xlsx from 
%      https://drive.google.com/file/d/0B51mihVyZFnVTXd3cm1GZlBoLWc/view?usp=sharing
% (3) GDSC_v6_Cell_Lines_Details from
%      https://drive.google.com/file/d/0B51mihVyZFnVdFFOVzdYSE1vNVE/view?usp=sharing

DD=1:15;
In=nchoosek(DD,2);
NN=25;
for II=1:size(In,1)
    Drug=[In(II,1) In(II,2)];
    [CCLE_Cell_lines,CS_Drugs_AUC]=CCLE_find_Common(Drug);
    [GDSC_Cell_lines,GS_Drugs_AUC]=GDSC_find_Common(Drug);
    
    ecop_DrugPair_CCLE{II}=FindOtherCopula(CS_Drugs_AUC,NN,'Gaussian');
    ecop_DrugPair_GDSC{II}=FindOtherCopula(GS_Drugs_AUC,NN,'Gaussian');
    Fro_drugPair1(II,1)=norm(ecop_DrugPair_CCLE{II}-ecop_DrugPair_GDSC{II}, 'fro');
    
    CS_Drugs_AUC_FF=[CS_Drugs_AUC(:,1) CS_Drugs_AUC(randperm(size(CS_Drugs_AUC,1)),2)];
    ecop_DrugPair_CCLE_rand{II}=FindOtherCopula(CS_Drugs_AUC_FF,NN,'Gaussian');
    GS_Drugs_AUC_FF=[GS_Drugs_AUC(randperm(size(GS_Drugs_AUC,1)),1) GS_Drugs_AUC(:,2)];
    ecop_DrugPair_GDSC_rand{II}=FindOtherCopula(GS_Drugs_AUC_FF,NN,'Gaussian');
    
    Fro_drugPair1(II,2)=norm(ecop_DrugPair_CCLE_rand{II}-ecop_DrugPair_GDSC{II}, 'fro');
    Fro_drugPair1(II,3)=norm(ecop_DrugPair_CCLE{II}-ecop_DrugPair_GDSC_rand{II}, 'fro');
end

[N_order,edges_order] = histcounts(Fro_drugPair1(:,1),11);
[N_disorder,edges_disorder] = histcounts([Fro_drugPair1(:,2); Fro_drugPair1(:,3)],21);
figure
b3 = bar(edges_order(1:end-1),(N_order/max(N_order))*max(N_disorder),'g');
b3.FaceAlpha = 0.5;
hold on
b4 = bar(edges_disorder(1:end-1),N_disorder, 'FaceColor',[0,0.7,0.7]);    
b4.FaceAlpha= 0.5;
axis([-0.2 2.5 0 50])
xlabel('Frobenius norm difference of copulas')
ylabel('Distribution of Frobenius norm difference')
legend('Difference between copulas of drugpairs of CCLE and GDSC with ordered cells',...
        'Difference between copulas of drugpairs of CCLE and GDSC with disordered cells')