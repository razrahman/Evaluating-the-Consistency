% This code will generate figure 4 of the paper. Before running it,
% download the following 3 files and keep it in the same folder with main code
% (1) CCLE_Drug_Sensitivity_Raziur_Extended.xlsx from
%      https://drive.google.com/file/d/0B51mihVyZFnVNzdvRWR4SVFodFk/view?usp=sharing
% (2) GDSC_v6_v17_fitted_dose_response.xlsx from 
%      https://drive.google.com/file/d/0B51mihVyZFnVTXd3cm1GZlBoLWc/view?usp=sharing
% (3) GDSC_v6_Cell_Lines_Details from
%      https://drive.google.com/file/d/0B51mihVyZFnVdFFOVzdYSE1vNVE/view?usp=sharing

DD=1:15;
In=nchoosek(DD,2);
NN1=100; % number of bootstrap sampling
for II=1:size(In,1)
    Drug=[In(II,1) In(II,2)];
    [CCLE_Cell_lines,CS_Drugs_AUC]=CCLE_find_Common(Drug);
    [GDSC_Cell_lines,GS_Drugs_AUC]=GDSC_find_Common(Drug);
    
    Corr_Pearson(II,:) = [corr(CS_Drugs_AUC(:,1),CS_Drugs_AUC(:,2),'type','Pearson') corr(GS_Drugs_AUC(:,1),GS_Drugs_AUC(:,2),'type','Pearson')];
    
    [~,bootsam] = bootstrp(NN1,[],CS_Drugs_AUC);
    for PP=1:NN1
        CS_Drugs_AUC_boot=CS_Drugs_AUC(bootsam(:,PP),:);
        Corr_CCLE_bootsam(PP,II)=corr(CS_Drugs_AUC_boot(:,1),CS_Drugs_AUC_boot(:,2),'type','Pearson');
    end
    PDF = fitdist(Corr_CCLE_bootsam(:,II),'Normal');
    CI = paramci(PDF,'Alpha',0.05);
    Corr_CCLE_bootsam_lower(II,:)=CI(1,1);
    Corr_CCLE_bootsam_upper(II,:)=CI(2,1);    
    if isempty(find(Corr_Pearson(II,2)>min(Corr_CCLE_bootsam(:,II)) & Corr_Pearson(II,2)<max(Corr_CCLE_bootsam(:,II)), 1 ))
        KK(II)=0;
    else
        KK(II)=1;
    end
end
fprintf('\n Out of 105 case, in %d times GDSC response pairs correlation coefficient is inside the box.\n',sum(KK));

figure
h1=boxplot(Corr_CCLE_bootsam,'Notch','on', 'colors','b');
hold on
scatter(1:size(Corr_CCLE_bootsam,2),Corr_Pearson(:,2),'g*');
xlabel('105(=15C2) combinations of drug pairs')
ylabel('Pearson Correlation Coefficient for Drug pairs')
