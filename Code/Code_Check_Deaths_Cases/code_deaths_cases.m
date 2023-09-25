clc;
clear all;




 filename1 = ['UniformBatch1month1R03.5VC0.9REP1000.mat'];
 filename2 = ['UniformBatch2month1R03.5VC0.9REP1000.mat'];
 filename3 = ['UniformBatch3month1R03.5VC0.9REP1000.mat'];


 tzero1 =  table2array(struct2table(load(filename1,'tzero')));
 tzero2 =  table2array(struct2table(load(filename2,'tzero')));
 tzero3 =  table2array(struct2table(load(filename3,'tzero')));
 tzero1=[tzero1;tzero2;tzero3];



 Czero1 =  table2array(struct2table(load(filename1,'Czero')));
 Czero2 =  table2array(struct2table(load(filename2,'Czero')));
 Czero3 =  table2array(struct2table(load(filename3,'Czero')));
 Czero1 = [Czero1;Czero2;Czero3];
%Czero1 = [Czero1;Czero2];

 CDzero1 =  table2array(struct2table(load(filename1,'CDzero')));
 CDzero2 =  table2array(struct2table(load(filename2,'CDzero')));
 CDzero3 =  table2array(struct2table(load(filename3,'CDzero')));
 CDzero1 = [CDzero1;CDzero2;CDzero3];




N =  table2array(struct2table(load(filename1,'N')));
i=5;
CDzero1(CDzero1<=100) = 0;
% CDzero2(CDzero2<=10) = 0;

fraction_death1 = CDzero1./sum(N) ;
% fraction_death2 = CDzero2./Czero2;

death1 = fraction_death1.*1e5;
% death2 = fraction_death2.*1e5;

death1(isnan(death1)) = 0;
% death2(isnan(death2)) = 0;

death1(isinf(death1)) = 0;
% death2(isinf(death2)) = 0;

countzero1 = zeros(3000,i);
% countzero2 = zeros(1000,11);

countzero1(death1(:)>0)=1;
% countzero2(death2(:)>0)=1;

count1 = sum(countzero1,1);
% count2 = sum(countzero2,1);



findD1 = find(death1(:,1)>0);
findD2 = find(death1(:,2)>0);
findD3 = find(death1(:,3)>0);
findD4 = find(death1(:,4)>0);
% findD5 = find(death1(:,5)>0);
% findD6 = find(death1(:,6)>0);
% findD7 = find(death1(:,7)>0);
% findD8 = find(death1(:,8)>0);
% findD9 = find(death1(:,9)>0);
% findD10 = find(death1(:,10)>0);
% findD11 = find(death1(:,11)>0);

t1= size(findD1);
t2= size(findD2);
t3= size(findD3);
t4= size(findD4);
% t5= size(findD5);
% t6= size(findD6);
% t7= size(findD7);
% t8= size(findD8);
% t9= size(findD9);
% t10= size(findD10);
% t11= size(findD11);


for j1 = 1:t1(1)
    
    D1(j1) = death1(findD1(j1),1);

end

for j2 = 1:t2(1)
    
    D2(j2) = death1(findD2(j2),2);

end

for j3 = 1:t3(1)
    
    D3(j3) = death1(findD3(j3),3);

end

for j4 = 1:t4(1)
    
    D4(j4) = death1(findD4(j4),4);

end

% for j5 = 1:t5(1)
% 
%     D5(j5) = death1(findD5(j5),5);
% 
% end

meandeath1(1) = mean(D1);
meandeath1(2) = mean(D2);
meandeath1(3) = mean(D3);
meandeath1(4) = mean(D4);
%meandeath1(5) = mean(D5);

meandeath1(isnan(meandeath1)) = 0;
meandeath1(isinf(meandeath1)) = 0;
Derror(1) = std(D1)/sqrt(t1(1));
Derror(2) = std(D2)/sqrt(t2(1));
Derror(3) = std(D3)/sqrt(t3(1));
Derror(4) = std(D4)/sqrt(t4(1));
%Derror(5) = std(D5)/sqrt(t5(1));

Czero1(Czero1<=100) = 0;

fraction_C1 = Czero1./sum(N);

C1 = fraction_C1.*1e5;
% death2 = fraction_death2.*1e5;

C1(isnan(C1)) = 0;
% death2(isnan(death2)) = 0;

C1(isinf(C1)) = 0;
% death2(isinf(death2)) = 0;

countzeroC1 = zeros(3000,i);
% countzero2 = zeros(1000,11);

countzeroC1(C1(:)>0)=1;
% countzero2(death2(:)>0)=1;

countC1 = sum(countzeroC1,1);
% count2 = sum(countzero2,1);

% meanC1 = sum(C1)./countC1;
% meanC1(isnan(meanC1)) = 0;
% meanC1(isinf(meanC1)) = 0;

findC1 = find(C1(:,1)>0);
findC2 = find(C1(:,2)>0);
findC3 = find(C1(:,3)>0);
findC4 = find(C1(:,4)>0);
% findC5 = find(C1(:,5)>0);
% findC6 = find(C1(:,6)>0);
% findC7 = find(C1(:,7)>0);
% findC8 = find(C1(:,8)>0);
% findC9 = find(C1(:,9)>0);
% findC10 = find(C1(:,10)>0);
% findC11 = find(C1(:,11)>0);
t1= size(findC1);
t2= size(findC2);
t3= size(findC3);
t4= size(findC4);
% t5= size(findC5);
% t6= size(findC6);
% t7= size(findC7);
% t8= size(findC8);
% t9= size(findC9);
% t10= size(findC10);
% t11= size(findC11);


for j1 = 1:t1(1)
    
    C0(j1) = C1(findC1(j1),1);

end

for j2 = 1:t2(1)
    
    C2(j2) = C1(findC2(j2),2);

end

for j3 = 1:t3(1)
    
    C3(j3) = C1(findC3(j3),3);

end

for j4 = 1:t4(1)
    
    C4(j4) = C1(findC4(j4),4);

end

% for j5 = 1:t5(1)
% 
%     C5(j5) = C1(findC5(j5),5);
% 
% end

% 
% for j6 = 1:t6(1)
%     
%     C6(j6) = C1(findC6(j6),6);
% 
% end
% 
% for j7 = 1:t7(1)
%     
%     C7(j7) = C1(findC7(j7),7);
% 
% end
% 
% for j8 = 1:t8(1)
%     
%     C8(j8) = C1(findC8(j8),8);
% 
% end
% 
% for j9 = 1:t9(1)
%     
%     C9(j9) = C1(findC9(j9),9);
% 
% end
% 
% for j10 = 1:t10(1)
%     
%     C10(j10) = C1(findC10(j10),10);
% 
% end
% 
% for j11 = 1:t11(1)
%     
%     C11(j11) = C1(findC11(j11),11);
% 
% end

meanC1(1) = mean(C0);
meanC1(2) = mean(C2);
meanC1(3) = mean(C3);
meanC1(4) = mean(C4);
% meanC1(5) = mean(C5);
% meanC1(6) = mean(C6);
% meanC1(7) = mean(C7);
% meanC1(8) = mean(C8);
% meanC1(9) = mean(C9);
% meanC1(10) = mean(C10);
% meanC1(11) = mean(C11);
% meanC1(isnan(meanC1)) = 0;
% meanC1(isinf(meanC1)) = 0;

Cerror(1) = std(C0)/sqrt(t1(1));
Cerror(2) = std(C2)/sqrt(t2(1));
Cerror(3) = std(C3)/sqrt(t3(1));
Cerror(4) = std(C4)/sqrt(t4(1));
% Cerror(5) = std(C5)/sqrt(t5(1));
% Cerror(6) = std(C6)/sqrt(t6(1));
% Cerror(7) = std(C7)/sqrt(t7(1));
% Cerror(8) = std(C8)/sqrt(t8(1));
% Cerror(9) = std(C9)/sqrt(t9(1));
% Cerror(10) = std(C10)/sqrt(t10(1));
% Cerror(11) = std(C11)/sqrt(t11(1));

a = [0.25,0.5,0.75,0.9];

figure(1)
errorbar(a,meandeath1,Derror,'--o','LineWidth',1.8,'MarkerSize',8);
colororder({'#fab73d','#e9723d','#259086','#c6d763','#6e93d6','#84817d'});
% legend('batch1','batch2','FontSize',12);
set(gca,'FontSize',14,'FontName','Helvetica'); %ylim([0 1.1]); %xlim([0 10]);
xlabel('Vaccinated Population (%)','FontName','Helvetica','FontSize',16);
ylabel('Deaths/100,000 cases','FontName','Helvetica','FontSize',16);
% title('Delta Variant','FontName','Helvetica','FontSize',12);
savefig(figure(1),'allR0Graph_D');

figure(2)
errorbar(a,meanC1,Cerror,'--o','LineWidth',1.8,'MarkerSize',8);
colororder({'#fab73d','#e9723d','#259086','#c6d763','#6e93d6','#84817d'});
% legend('batch1','batch2','FontSize',12);
set(gca,'FontSize',14,'FontName','Helvetica'); %ylim([0 1.1]); %xlim([0 10]);
xlabel('Vaccinated Population (%)','FontName','Helvetica','FontSize',16);
ylabel('Cases/100,000 population','FontName','Helvetica','FontSize',16);
% title('Delta Variant','FontName','Helvetica','FontSize',12);
savefig(figure(2),'allR0Graph_C');

save('Uniformmonth1_Death.mat');