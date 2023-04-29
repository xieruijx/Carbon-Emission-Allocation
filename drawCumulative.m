clear
close all

load Result_CEF_ES_Feedback_Loss.mat
CCostAccum_CEF = CCostAccum;
load Result_Lyapunov_ES_Feedback_Loss_V1_k0.05.mat
CCostAccum_AS = CCostAccum;

figure;
plot(0:Num_T,CCostAccum_AS(:,1)*1e5,0:Num_T,CCostAccum_CEF(:,1)*1e5,'LineWidth',2);
legend('Proposed','CEF','location','northwest');
xlabel('Time (h)');
ylabel('Cumulative emission of ES 15 (kgCO_2)');
set(gca,'FontName','Times New Roman','FontSize',14);

figure;
set(gcf,'position',[100,100,560,270])
subplot(1,2,1);
plot(0:Num_T,CCostAccum_AS(:,1)*1e5,'LineWidth',2);
legend('Proposed','location','southwest');
xlabel('Time (h)');
ylabel('Cumulative emission of ES 15 (kgCO_2)');
set(gca,'FontName','Times New Roman','FontSize',10);

subplot(1,2,2);
plot(0:Num_T,CCostAccum_CEF(:,1)*1e5,'r','LineWidth',2);
legend('CEF','location','southeast');
xlabel('Time (h)');
ylabel('Cumulative emission of ES 15 (kgCO_2)');
set(gca,'FontName','Times New Roman','FontSize',10);