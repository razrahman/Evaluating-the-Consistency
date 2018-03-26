% This code will generate figure 7 and 8 of the paper. Before running it,
% download the following 2 files and keep it in the same folder with main code
% (1) GDSC_Data_v6_(Common_Genes_&_Cell-lines)_SRD_Mar_01.xlsx from
%      https://drive.google.com/file/d/0B51mihVyZFnVYUJFS3lPbGtzUVk/view?usp=sharing
% (2) CCLE_Data_(Common_Genes_&_Cell-lines)_SRD_Mar_01.xlsx from
%      https://drive.google.com/file/d/0B51mihVyZFnVSTJiTFVWWHpwazg/view?usp=sharing
[GG_num,GG_txt]=xlsread('GDSC_Data_v6_(Common_Genes_&_Cell-lines)_SRD_Mar_01.xlsx');
[CG_num,CG_txt]=xlsread('CCLE_Data_(Common_Genes_&_Cell-lines)_SRD_Mar_01.xlsx');
%% Figure 7: Difference between copulas of identical cell lines of CCLE and GDSC with ordered and disordered genes
NN=25;
parfor RR=1:size(CG_num,1)
    Cop_cell{RR}    = FindOtherCopula([CG_num(RR,:)' GG_num(RR,:)'],NN,'Gaussian');
    Cop_cell_rC{RR} = FindOtherCopula([CG_num(RR,randperm(size(CG_num,2)))' GG_num(RR,:)'],NN,'Gaussian');
    Cop_cell_rG{RR} = FindOtherCopula([CG_num(RR,:)' GG_num(RR,randperm(size(CG_num,2)))'],NN,'Gaussian');
end
for RR1=1:size(CG_num,1)
    for RR2=1:size(CG_num,1)
        Fro_cell(RR1,RR2)=norm(Cop_cell{RR1}-Cop_cell{RR2}, 'fro');
    end
    Fro_cell_max(RR1,1)=max(Fro_cell(RR1,:));
    Fro_cell_max(RR1,2)=norm(Cop_cell{RR1}-Cop_cell_rC{RR1}, 'fro');
    Fro_cell_max(RR1,3)=norm(Cop_cell{RR1}-Cop_cell_rG{RR1}, 'fro');
end
Order_gene = Fro_cell(~eye(size(Fro_cell)));
Disorder_gene=[Fro_cell_max(:,2); Fro_cell_max(:,3)];

[N_order,edges_order] = histcounts(Order_gene,101);
[N_disorder,edges_disorder] = histcounts(Disorder_gene,101);

figure
bar(edges_order(1:end-1),(N_order/max(N_order))*max(N_disorder),'g')
hold on
bar(edges_disorder(1:end-1),N_disorder)
axis([0 2.5 0 110])
xlabel('Frobenius norm difference of copulas')
ylabel('Distribution of Frobenius norm difference')
legend('Difference between copulas of identical cell lines of CCLE and GDSC with ordered genes',...
        'Difference between copulas of identical cell lines of CCLE and GDSC with disordered genes')
%% Figure 8: Difference between copulas of identical genes of CCLE and GDSC with ordered and disordered cell lines
for RR=1:size(CG_num,2)
    Cop_cell2{RR}    = FindOtherCopula([CG_num(:,RR) GG_num(:,RR)],NN,'Gaussian');
    Cop_cell2_rC{RR} = FindOtherCopula([CG_num(randperm(size(CG_num,1)),RR) GG_num(:,RR)],NN,'Gaussian');
    Cop_cell2_rG{RR} = FindOtherCopula([CG_num(:,RR) GG_num(randperm(size(CG_num,1)),RR)],NN,'Gaussian');
end
for RR1=1:size(CG_num,2)
    for RR2=1:size(CG_num,2)
        Fro2_cell(RR1,RR2)=norm(Cop_cell2{RR1}-Cop_cell2{RR2}, 'fro');
    end
    Fro_cell_max2(RR1,1)=max(Fro2_cell(RR1,:));
    Fro_cell_max2(RR1,2)=norm(Cop_cell2{RR1}-Cop_cell2_rC{RR1}, 'fro');
    Fro_cell_max2(RR1,3)=norm(Cop_cell2{RR1}-Cop_cell2_rG{RR1}, 'fro');
end

Order_cell = Fro2_cell(~eye(size(Fro2_cell)));
Disorder_cell=[Fro_cell_max2(:,2); Fro_cell_max2(:,3)];

[N_order2,edges_order2] = histcounts(Order_cell,101);
[N_disorder2,edges_disorder2] = histcounts(Disorder_cell,101);

figure
b1 = bar(edges_order2(1:end-1),(N_order2/max(N_order2))*max(N_disorder2),'g');
b1.FaceAlpha = 0.5;
hold on
b2 = bar(edges_disorder2(1:end-1),N_disorder2, 'FaceColor',[0,0.7,0.7]);  
b2.FaceAlpha= 0.5;
axis([0 3 0 800])
xlabel('Frobenius norm difference of copulas')
ylabel('Distribution of Frobenius norm difference')
legend('Difference between copulas of identical genes of CCLE and GDSC with ordered cell lines',...
        'Difference between copulas of genes with disordered cell lines of CCLE and GDSC')