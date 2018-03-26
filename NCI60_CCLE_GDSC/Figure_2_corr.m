% This code will generate figure 4 of the paper. Before running it,
% run the following 3 codes which will generate the required files of this code
% (1) CCLE_GDSC_CORR.R
% (2) CCLE_NCI60_CORR.R
% (3) NCI60_GDSC_CORR.R
%% CCLE GDSC 
[CCLE_GDSC_num, CCLE_GDSC_txt, CCLE_GDSC_raw] = xlsread('corr_data.xlsx');
figure
subplot(331)
bar(CCLE_GDSC_num(:,3:4)); hold on
hline = refline([0 .5]);
hline.Color = 'r';
set(gca,'XTick',1:16)
set(gca,'XTickLabel',CCLE_GDSC_txt(4:19,1))
set(gca,'XTickLabelRotation',90)
axis([0 17 -0.2 1])
legend('Pearson Correlation','Spearman Correlation')
title('Direct Correlation')
ylabel('CCLE & GDSC')
hold off

subplot(332)
bar(CCLE_GDSC_num(:,6:7)); hold on
hline = refline([0 .5]);
hline.Color = 'r';
set(gca,'XTick',1:16)
set(gca,'XTickLabel',CCLE_GDSC_txt(4:19,1))
set(gca,'XTickLabelRotation',90)
axis([0 17 -0.2 1])
legend('Pearson Correlation','Spearman Correlation')
title('Range Adjusted Correlation')
hold off

subplot(333)
bar(CCLE_GDSC_num(:,10:11)); hold on
hline = refline([0 .5]);
hline.Color = 'r';
set(gca,'XTick',1:16)
set(gca,'XTickLabel',CCLE_GDSC_txt(4:19,1))
set(gca,'XTickLabelRotation',90)
axis([0 17 -0.2 1])
legend('Pearson Correlation','Spearman Correlation')
title('Log Converted Correlation')
ylabel('Correlation Coefficient')
hold off

%% CCLE NCI60 
[CCLE_NCI60_num, CCLE_NCI60_txt, CCLE_NCI60_raw] = xlsread('corr_data.xlsx',2);
% figure
% title('Correlation Coefficient of CCLE and NCI60 IC_{50} values')
subplot(334)
bar(CCLE_NCI60_num(:,3:4)); hold on
set(gca,'XTick',1:10)
set(gca,'XTickLabel',CCLE_NCI60_txt(2:11,1))
set(gca,'XTickLabelRotation',90)
axis([0 11 -0.4 0.4])
legend('Pearson Correlation','Spearman Correlation')
title('Direct Correlation')
ylabel('CCLE & NCI60')
hold off

subplot(335)
bar(CCLE_NCI60_num(:,6:7)); hold on
set(gca,'XTick',1:10)
set(gca,'XTickLabel',CCLE_NCI60_txt(2:11,1))
set(gca,'XTickLabelRotation',90)
axis([0 11 -0.4 0.4])
legend('Pearson Correlation','Spearman Correlation')
title('Range Adjusted Correlation')
hold off

subplot(336)
bar(CCLE_NCI60_num(:,9:10)); hold on
set(gca,'XTick',1:10)
set(gca,'XTickLabel',CCLE_NCI60_txt(2:11,1))
set(gca,'XTickLabelRotation',90)
axis([0 11 -0.4 0.4])
legend('Pearson Correlation','Spearman Correlation')
title('Log Converted Correlation')
ylabel('Correlation Coefficient')
hold off

%% NCI60 GDSC 
[NCI60_GDSC_num, NCI60_GDSC_txt, NCI60_GDSC_raw] = xlsread('corr_data.xlsx',3);
% title('Correlation Coefficient of NCI60 and GDSC IC_{50} values')
subplot(337)
hist(NCI60_GDSC_num(:,3:4),10); hold on
h1= histfit(NCI60_GDSC_num(:,3),10); 
set(h1(1),'facecolor','b'); set(h1(2),'color','r','linestyle','--')
delete(h1(1))
h2 =histfit(NCI60_GDSC_num(:,4),10);
set(h2(1),'facecolor','y'); set(h2(2),'color','g','linestyle','--')
delete(h2(1))
axis([-1 1 0 50])
% legend('Pearson Correlation','Spearman Correlation')
title('Direct Correlation')
ylabel('NCI60 & GDSC')
xlabel('Correlation Coefficient')
hold off

subplot(338)
hist(NCI60_GDSC_num(:,8:9),10); hold on
h1= histfit(NCI60_GDSC_num(:,8),10); hold on
set(h1(1),'facecolor','b'); set(h1(2),'color','r','linestyle','--')
delete(h1(1))
h2 =histfit(NCI60_GDSC_num(:,9),10);
set(h2(1),'facecolor','y'); set(h2(2),'color','g','linestyle','--')
delete(h2(1))
axis([-1 1 0 50])
% legend('Pearson Correlation','Spearman Correlation','Distribution of drugs against Pearson correlation','Distribution of drugs against Spearman correlation')
title('Range Adjusted Correlation')
xlabel('Correlation Coefficient')
% ylabel('Number of drugs with that correlation coefficient range')
hold off

subplot(339)
hist(NCI60_GDSC_num(:,12:13),10); hold on
h1= histfit(NCI60_GDSC_num(:,12),10); 
set(h1(1),'facecolor','b'); set(h1(2),'color','r','linestyle','--')
delete(h1(1))
h2 =histfit(NCI60_GDSC_num(:,13),10);
set(h2(1),'facecolor','y'); set(h2(2),'color','g','linestyle','--')
delete(h2(1))
axis([-1 1 0 50])
% legend('Pearson Correlation','Spearman Correlation')
title('Log Converted Correlation')
xlabel('Correlation Coefficient')
ylabel('Number of drugs')
hold off