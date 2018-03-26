function [Cell_lines_sensitivity,Drugs_AUC]=GDSC_find_Common(Drug)
GDSC_drugs=[1026 38 1062 1 119 1013 1047 11 1060 1054 37 6 1036 30 35];
[GS_num,~,~]=xlsread('GDSC_v6_v17_fitted_dose_response');
[cell_num,cell_txt,~]=xlsread('C:\Users\Raziur_Rahman\Google Drive\Data CCLE and GDSC\GDSC\Version 6\GDSC_v6_Cell_Lines_Details');

Ind_D1=find(GDSC_drugs(Drug(1))==GS_num(:,4));
D1_cosmic=GS_num(Ind_D1,3);
kk1=0;
for ii=1:length(D1_cosmic)
    if ~isempty(cell_txt(find(D1_cosmic(ii)==cell_num(:,2))+1,1))
        kk1=kk1+1;
        GS_Drug1_cell_lines(kk1,:)=cell_txt(find(D1_cosmic(ii)==cell_num(:,2))+1,1);
        GS_Drug1_AUC(kk1,:)=1-GS_num(Ind_D1(ii),7);
    end
end
kk2=0;
Ind_D2=find(GDSC_drugs(Drug(2))==GS_num(:,4));
D2_cosmic=GS_num(Ind_D2,3);
for jj=1:length(D2_cosmic)
    if ~isempty(cell_txt(find(D2_cosmic(jj)==cell_num(:,2))+1,1))
        kk2=kk2+1;
        GS_Drug2_cell_lines(kk2,:)=cell_txt(find(D2_cosmic(jj)==cell_num(:,2))+1,1);
        GS_Drug2_AUC(kk2,:)=1-GS_num(Ind_D2(jj),7);
    end
end

[Cell_lines_sensitivity,Ind_drug1,Ind_drug2]=intersect(GS_Drug1_cell_lines,GS_Drug2_cell_lines);
Drugs_AUC=[GS_Drug1_AUC(Ind_drug1) GS_Drug2_AUC(Ind_drug2)];
Cell_lines_sensitivity=regexprep(Cell_lines_sensitivity(1:end,:), '-', '');
Cell_lines_sensitivity=regexprep(Cell_lines_sensitivity(1:end), ' ', '');
