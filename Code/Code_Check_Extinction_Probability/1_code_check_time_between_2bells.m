close all 
clear all

VE0 = 95;
R0 = 2.5;
sR0 = size(R0);
l = 1000; aa = [0.25,0.5,0.75,0.9]; aasize = size(aa);
N = 67529761;
I = 1; % initial infectious individuals
N = N-I;
L = 0;
A = 0;
R = 0;
t1 = 0:1e-1:3000;
st = size(t1);
%eS = 0.953; eI = 0.9;
%eS = 0.546; eI = 0.770;
eS = 0.552; eI = 0.032;
% tzero = NaN(l,sR0(2));
% tcount = 0;
data_Ncount= NaN(11,100);
data_edges= NaN(11,100);

filename = 'UniformBatch1month1R03.5VC0.9REP1000.mat';
    d1 = load(filename,'tzero');
data1 =  table2array(struct2table(d1));
   
filename = 'UniformBatch2month1R03.5VC0.9REP1000.mat';
 d2 = load(filename,'tzero');
data2 =  table2array(struct2table(d2));

data = [data1;data2];
    
for i=1:aasize(2)
    S = N - (aa(i)*N);
    SV = eS*aa(i)*N;
    SP = (1-eS)*aa(i)*N;

    data_tzero = data(:,i);
    data_tzero_total(:,i) = data_tzero;

data_Pext(i) = 100 - (100.*sum(isnan(data_tzero))/l);
data_tzeromax(i) = max(data_tzero);
if i ==1
s(i,:) = size(find(data_tzero_total(:,i)< 250));
data_Pext_cri(i,1) = (s(i,1)/l);
data_Pout_cri(i,1) = 1 - (s(i,1)/l);

elseif i>=2 && i<5
    
s(i,:) = size(find(data_tzero_total(:,i)< 300));
data_Pext_cri(i,1) = (s(i,1)/l);
data_Pout_cri(i,1) = 1 - (s(i,1)/l);
else 
    s(i,:) = size(find(data_tzero_total(:,i)< 600));
data_Pext_cri(i,1) = (s(i,1)/l);
data_Pout_cri(i,1) = 1 - (s(i,1)/l);
end
figure(i)
hist(data_tzero(:,1),500);
% legend('Simulation','Theory','C');
set(gca,'FontSize',14,'FontName','monospaced');
xlabel('Days','FontName','Helvetica','FontSize',16); xlim([0 500]);
ylabel('Frequency','FontName','Helvetica','FontSize',16);
% title('R_0 = 1.2 dt = 0.1','FontName','Helvetica','FontSize',12);
name = ['dt01','R0',num2str(R0),'VC',num2str(aa(i)) '.fig'];
%savefig(figure(i),name);

[Ncount,edeges,bins] = histcounts(data_tzero(:,1),500);
Ncount = Ncount/sum(Ncount); 

for n = 1:500
    bin_means(:,n) = mean(data_tzero(bins==n,1)); %average Y in each bin/subinterval
%     bin_mdpt 
end

figure(13)
semilogy(bin_means,Ncount);
hold on

end
% figure(13)
% legend('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%');
% set(gca,'FontSize',12,'FontName','monospaced');
% xlabel('Days','FontName','Helvetica','FontSize',14); xlim([0 800]);
% ylabel('Probability of Extinction','FontName','Helvetica','FontSize',14);
% title('R_0 = 2.5 dt = 0.1','FontName','Helvetica','FontSize',12);
% name = ['histfitdt01','R0',num2str(R0),'VC',num2str(aa(i)) '.fig'];
% % savefig(figure(13),name);
% 
% VC = 0:10:100;
% figure(12)
% plot(VC,data_Pext_cri,'LineWidth',1.2);
% set(gca,'FontSize',12,'FontName','monospaced'); ylim([0 1]);
% xlabel('% Vaccinated Population','FontName','Helvetica','FontSize',14);
% ylabel('Extinction Probability','FontName','Helvetica','FontSize',14);
% title('R_0 = 2.5 dt = 0.1','FontName','Helvetica','FontSize',12);
% savefig(figure(12),'SummaryGraph');

% save('R025Hetdt01.mat');