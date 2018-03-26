% This code will generate figure 9 of the paper. Before running it,
% download the following 2 files and keep it in the same folder with main code
% (1) CCLE_Drug_Sensitivity_Raziur_Extended.xlsx from
%      https://drive.google.com/file/d/0B51mihVyZFnVNzdvRWR4SVFodFk/view?usp=sharing
% (2) GDSC_drug_common from
%      https://drive.google.com/file/d/1MWr4sY3DQAH7XZYNgOQ3ZRCETXQ-EeRz/view?usp=sharing
DD=1:15;
NN=25;
for Drug=1:length(DD)
    GDSC_drug_sensitivity=[];    GDSC_cell_lines=[];
    [GS_num, GS_txt]=xlsread('GDSC_drug_common.xlsx');
    GDSC_drug_sensitivity1=1-GS_num(8:end,4+Drug);   
    T=[];
    T=isnan(sum(GDSC_drug_sensitivity1,2));
    CellLineIndex=[];
    CellLineIndex=find(T==0);
    GDSC_drug_sensitivity=GDSC_drug_sensitivity1(CellLineIndex,:);
    GDSC_cell_lines=GS_txt(7:end,1);
    GDSC_cell_lines=GDSC_cell_lines(CellLineIndex,:);
    
    CCLE_drug_sensitivity=[];     CCLE_cell_lines=[];
    [CS_num, CS_txt]=xlsread('CCLE_Drug_Sensitivity_Raziur_Extended.xlsx',Drug);
    CCLE_drug_sensitivity=CS_num(1:end,10)/8;
    CCLE_cell_lines=CS_txt(2:end,2);
    g(Drug)=0;
    cell_line_index_GDSC=[];
    cell_line_index_CCLE=[];
    h=[];
    for i=1:length(GDSC_cell_lines)
        for j=1:length(CCLE_cell_lines)
            h(i,j)=strcmp(CCLE_cell_lines(j),GDSC_cell_lines(i));
            
        end
        if any(h(i,:))==1
            g(Drug)=g(Drug)+1;
            cell_line_index_GDSC(g(Drug),:)=i;
            cell_line_index_CCLE(g(Drug),:)=find(h(i,:)==1);
        end
    end
    finalY_GDSC{Drug}=GDSC_drug_sensitivity(cell_line_index_GDSC,:);
    finalY_CCLE{Drug}=CCLE_drug_sensitivity(cell_line_index_CCLE(1:end,:),:);
    %%    
    Y1=[finalY_CCLE{Drug} finalY_GDSC{Drug}];
    ecop_Drug{Drug}=FindOtherCopula(Y1,NN,'Gaussian');
    
    Y1_CC=[Y1(randperm(size(Y1,1)),1) Y1(:,2)];
    ecop_CCLE_rand{Drug}=FindOtherCopula(Y1_CC,NN,'Gaussian');
    
    Y1_GD=[Y1(:,1) Y1(randperm(size(Y1,1)),2)];
    ecop_GDSC_rand{Drug}=FindOtherCopula(Y1_GD,NN,'Gaussian');
end

for RR1=1:15
    for RR2=1:15
        Fro_drug(RR1,RR2)=norm(ecop_Drug{RR1}-ecop_Drug{RR2}, 'fro');
    end
    Fro_drug_max(RR1,1)=max(Fro_drug(RR1,:));
    Fro_drug_max(RR1,2)=norm(ecop_Drug{RR1}-ecop_CCLE_rand{RR1}, 'fro');
    Fro_drug_max(RR1,3)=norm(ecop_Drug{RR1}-ecop_GDSC_rand{RR1}, 'fro');
end
%%
Fro_drug2 = Fro_drug(~eye(size(Fro_drug)));
[N_order,edges_order] = histcounts(Fro_drug2,11);
[N_disorder,edges_disorder] = histcounts([Fro_drug_max(:,2); Fro_drug_max(:,3)],11);

figure
b1 = bar(edges_order(1:end-1),(N_order/max(N_order))*max(N_disorder),'g');
b1.FaceAlpha = 0.5;
hold on
b2 = bar(edges_disorder(1:end-1),N_disorder, 'FaceColor',[0,0.7,0.7]);    
b2.FaceAlpha= 0.5;
axis([-0.1 1.4 0 7])
xlabel('Frobenius norm difference of copulas')
ylabel('Distribution of Frobenius norm difference')
legend('Difference between copulas of identical drugs of CCLE and GDSC with ordered cells',...
        'Difference between copulas of identical drugs of CCLE and GDSC with disordered cells')
