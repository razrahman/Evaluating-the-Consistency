% This code will generate figure 3 of the paper. Before running it,
% download the following 5 files and keep it in the same folder with main code
% (1) CCLE_Drug_Sensitivity_Raziur_Extended.xlsx from
%      https://drive.google.com/file/d/0B51mihVyZFnVNzdvRWR4SVFodFk/view?usp=sharing
% (2) GDSC_v6_v17_fitted_dose_response.xlsx from 
%      https://drive.google.com/file/d/0B51mihVyZFnVTXd3cm1GZlBoLWc/view?usp=sharing
% (3) GDSC_v6_Cell_Lines_Details from
%      https://drive.google.com/file/d/0B51mihVyZFnVdFFOVzdYSE1vNVE/view?usp=sharing

DD=1:15;
In=nchoosek(DD,2);
Corr_Pearson = zeros(size(In,1),2);
for II=1:size(In,1)
    fprintf('\nCalculating information for Drug pair no %d out of 105 drug-pairs\n', II);
    Drug=[In(II,1) In(II,2)];
    [CCLE_Cell_lines,CS_Drugs_AUC]=CCLE_find_Common(Drug);
    [GDSC_Cell_lines,GS_Drugs_AUC]=GDSC_find_Common(Drug);
    
    Corr_Pearson(II,:) = [corr(CS_Drugs_AUC(:,1),CS_Drugs_AUC(:,2),'type','Pearson') corr(GS_Drugs_AUC(:,1),GS_Drugs_AUC(:,2),'type','Pearson')];
end
figure
scatter(Corr_Pearson(:,1),Corr_Pearson(:,2),'r*'); hold on
refline(1,0)
xlabel('Pearson correlation coefficient of CCLE drug pairs')
ylabel('Pearson correlation coefficient of GDSC drug pairs')
