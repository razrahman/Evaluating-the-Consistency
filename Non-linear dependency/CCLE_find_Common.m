function [Cell_lines_sensitivity,Drugs_AUC]=CCLE_find_Common(Drug)

[CS_num1,CS_txt1]=xlsread('CCLE_Drug_Sensitivity_Raziur_Extended.xlsx',Drug(1));
[CS_num2,CS_txt2]=xlsread('CCLE_Drug_Sensitivity_Raziur_Extended.xlsx',Drug(2));

CS_Drug1_cell_lines=CS_txt1(2:end,2);
CS_Drug1_AUC=CS_num1(:,6)/8;

CS_Drug2_cell_lines=CS_txt2(2:end,2);
CS_Drug2_AUC=CS_num2(:,6)/8;
[Cell_lines_sensitivity,Ind_drug1,Ind_drug2]=intersect(CS_Drug1_cell_lines,CS_Drug2_cell_lines);
Drugs_AUC=[CS_Drug1_AUC(Ind_drug1) CS_Drug2_AUC(Ind_drug2)];
Cell_lines_sensitivity=regexprep(Cell_lines_sensitivity(1:end,:), '-', '');
Cell_lines_sensitivity=regexprep(Cell_lines_sensitivity(1:end), ' ', '');
