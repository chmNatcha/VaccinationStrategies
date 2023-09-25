close all 
clear all

VC = [0.25,0.5,0.75,0.9]; % vaccine coverage

R0 = 2.5;
R0size = size(R0);
l = 1000; aa = VC; aasize = size(aa);
N = 67529761;
I = 1; % initial infectious individuals
N = N-I;
L = 0;
A = 0;
R = 0;
eS = 0.552; eI = 0.032;
FA = 0.4;
FI = 1-FA;
qA =0.58;
% tzero = NaN(l,sR0(2));
% tcount = 0;

filename = 'UniformBatch1month1R03.5VC0.9REP1000.mat';
    d1 = load(filename,'tzero');
data1 =  table2array(struct2table(d1));
   
filename = 'UniformBatch2month1R03.5VC0.9REP1000.mat';
 d2 = load(filename,'tzero');
data2 =  table2array(struct2table(d2));

 filename = 'UniformBatch3month1R03.5VC0.9REP1000.mat';
    d3 = load(filename,'tzero');
data3 =  table2array(struct2table(d3));
   

%filename = 'UniformBatch4boosterR02.5VC1REP1000.mat';
% d4 = load(filename,'tzero');
% data4 =  table2array(struct2table(d4));
% 
% filename = 'UniformBatch5boosterR02.5VC1REP1000.mat';
% d5 = load(filename,'tzero');
% data5 =  table2array(struct2table(d5));


for i=1:aasize(2)
% 
% if i <= 4
batch1(i,:) = size(find(data1(:,i)<100));
data_Pout_cri(i,1) = 1- (batch1(i,1)/1000);
data_Pext_cri(i,1) =  (batch1(i,1)/1000);

batch2(i,:) = size(find(data2(:,i)<100));
data_Pout_cri(i,2) = 1-(batch2(i,1)/1000);
data_Pext_cri(i,2) =  (batch2(i,1)/1000);

batch3(i,:) = size(find(data3(:,i)<100));
data_Pout_cri(i,3) = 1- (batch3(i,1)/1000);
data_Pext_cri(i,3) =  (batch3(i,1)/1000);

%batch4(i,:) = size(find(data4(:,i)<150));
%data_Pout_cri(i,4) = 1-(batch4(i,1)/1000);
%data_Pext_cri(i,4) =  (batch4(i,1)/1000);

%batch5(i,:) = size(find(data5(:,i)<150));
%data_Pout_cri(i,5) = 1- (batch5(i,1)/1000);
%data_Pext_cri(i,5) =  (batch5(i,1)/1000);


data_Pout_error(i,1) = std(data_Pout_cri(i,1:3))/sqrt(3);
data_Pout_mean(i,1) = mean(data_Pout_cri(i,1:3));

data_Pext_error(i,1) = std(data_Pext_cri(i,1:3))/sqrt(3);
data_Pext_mean(i,1) = mean(data_Pext_cri(i,1:3));

end

for j = 1:aasize(2)
    for k = 1:R0size(2)
    S(j,k) = 1 - aa(j);
    SP(j,k) = (1-eS)*aa(j);
   Rt(j,k) = ((R0(k).*S(j,k)) + ((1 - eI).*R0(k).*SP(j,k)));
%Rt(j,k) = R0(k).*(FI + qA*FA).*(1-eI*aa(j));

    Pext_theory(j,k) = (1./Rt(j,k));
    Pout_theory(j,k) = 1-Pext_theory(j,k);

      if Pext_theory(j,k) > 1
        Pext_theory(j,k) = 1;
        Pout_theory(j,k) = 0;
      end
      
    end
end
save('Uniformmonth1_Pext.mat');


%figure(1)
%plot(VC,data_Pout_cri,'o','LineWidth',1.4);
%hold on
%plot(VC,Pout_theory,'--','LineWidth',1.5);
%set(gca,'FontSize',12,'FontName','monospaced');ylim([0 1]);
%xlabel('% Vaccinated Population','FontName','Helvetica','FontSize',14);
%ylabel('Probability of Suscessful Outbreak','FontName','Helvetica','FontSize',14);
%title('R_0 = 3.5 dt = 0.1 1000 runs','FontName','Helvetica','FontSize',12);
% savefig(figure(1),'R025dt01_Pout');

figure(2)
plot(VC,data_Pext_cri,'o','LineWidth',1.4);
hold on
plot(VC,Pext_theory,'--','LineWidth',1.5);
set(gca,'FontSize',12,'FontName','monospaced');ylim([0 1]);
xlabel('% Vaccinated Population','FontName','Helvetica','FontSize',14);
ylabel('Probability of Extinction','FontName','Helvetica','FontSize',14);
title('R_0 = 3.5 dt = 0.1 1000 runs','FontName','Helvetica','FontSize',12);
% savefig(figure(2),'R025dt01_Pext');
