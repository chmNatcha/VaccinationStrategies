close all
clear all
clc
rng('shuffle');
%% Vaccination - Herd Immunity Model : Using stochastic method
% COVID-19 // UK 2020 // Well-mixed model
%% Set Time
dt = 0.1; %time step
t_final = 1000; % final time
total_step = round(t_final/dt);
Repetition = 1000;  % Round
N_save = total_step;
save_step = total_step;%% Constant parameters
count_r = 1;
gamma = 1/2.5; %1/infectous period
sigma = 1/3.0; %1/Latend period
FA = 0.4;
FI = 1 - FA;
qA = 0.58; % the asymptomatic infectiousness
%% Contact Matirx UK
num1 = readtable('UK_Data.xlsx','sheet','CM 16g','Range','A2:P17');
 c1 = table2array(num1); % Contact Matrix
%  c1 = ones(16);
%% %changing R0
R0_data = 2.5; % Change R0 here
R0 = R0_data;
e = max(eig(c1));
k1 = (R0*gamma)/(e*(FI+qA*FA));
k2 = (R0*gamma)/(e*(FI+qA*FA));
beta1 = k1.*c1;
beta2 = k2.*c1;
%%  Set Populations
group = 16; %Number of age groups
num2 = readtable('UK_Data.xlsx','sheet','Pop 16g','Range','A2:A17');
M1 = table2array(num2);
N = M1;
Ntot = zeros(total_step,Repetition);
%% Vaccine Efficacy Parameters
  % Change vaccine effectiveness here
        eS = 0.552; %efficacy against infection
        eI = 0.032; %efficacy against transmission
        eD = 0.487; %efficacy against diesase
        eDE = 0.95;  %effetiveness agaisnt death
      FAP = 1 - ((1-FA)*(1-eD)/(1-eS));% fraction of symptomatic vaccinated individuals
num3 = readtable('DataOmicron.xlsx','sheet','IFR','Range','I3:I18');
fD = table2array(num3);

C = (1-eDE)/(1-eD);
kIPDP = gamma.*fD.*C;
kIPRP = (gamma.*ones(16,1)) - kIPDP;
kID = gamma.*fD;
kIR = (gamma.*ones(16,1))- kID;

aadata1 =    readtable('OrderOmicron.xlsx','sheet','Fraction I 2.5','Range','G2:G17');  %load order for number of vaccine doses
aadata2 =    readtable('OrderOmicron.xlsx','sheet','Fraction I 2.5','Range','F2:F17');  %load order for vaccination 
order_vaccaine = table2array(aadata1); 
order_group = table2array(aadata2); 

aa = [0.25,0.5,0.75,0.9]; %Vaccine coverage %change vc here
ss = size(aa);
tzero = NaN(Repetition,ss(2));
Czero = zeros(Repetition,ss(2));
CDzero = zeros(Repetition,ss(2));
%% loop
tic
for mm = 1:ss(2)
Infected_save = zeros(Repetition,N_save); % L+I+IP+LP
% R0 = R0_data(mm);
% beta1 = R0*gamma_S;

S_save = zeros(Repetition,N_save);
L_save = zeros(Repetition,N_save);
I_save = zeros(Repetition,N_save);
A_save = zeros(Repetition,N_save);
R_save = zeros(Repetition,N_save);
D_save = zeros(Repetition,N_save);
C_save = zeros(Repetition,N_save);
CD_save = zeros(Repetition,N_save);
t_save = zeros(Repetition,N_save);


SP_save = zeros(Repetition,N_save);
LP_save = zeros(Repetition,N_save);
IP_save = zeros(Repetition,N_save);
AP_save = zeros(Repetition,N_save);
RP_save = zeros(Repetition,N_save);
DP_save = zeros(Repetition,N_save);
CP_save = zeros(Repetition,N_save);

S_data = zeros(group,N_save);
L_data = zeros(group,N_save);
I_data = zeros(group,N_save);
A_data = zeros(group,N_save);
R_data = zeros(group,N_save);
D_data = zeros(group,N_save);
C_data = zeros(group,N_save);
CD_data = zeros(group,N_save);

SP_data = zeros(group,N_save);
LP_data = zeros(group,N_save);
IP_data = zeros(group,N_save);
AP_data = zeros(group,N_save);
RP_data = zeros(group,N_save);
DP_data = zeros(group,N_save);
Infected_data = zeros(group,N_save);

for l = 1:Repetition
    % Saved Data Parameters

fA = zeros(group,1);
fAP = zeros(group,1);
fA2 = zeros(group,1);
fAP2 = zeros(group,1);
fI = zeros(group,1);
fIP = zeros(group,1);
fI2 = zeros(group,1);
fIP2 = zeros(group,1);

t_zero = zeros(ss(2),1);
%---------set time-------------
% dt = 0.001; %time step
t = 0; % initial time
% t_final = 100; % final time
step = 1;
TimeStochastic= zeros(total_step,1);
CheckTime = 0;
%----- Initial Conditions --------%
I = zeros(group,1); % initial infectious individuals
A = zeros(group,1);
r1 = rand;
if r1 <= FA
   r = randi([1 sum(M1)],1,1);
    if r <= M1(1)
    A(1) = 1;    
elseif r> M1(1) &&  r<= sum(M1(1:2))
    A(2) = 1;  
elseif r> sum(M1(1:2)) &&  r<= sum(M1(1:3))
    A(3) = 1;
elseif r> sum(M1(1:3)) &&  r<= sum(M1(1:4))  
    A(4) = 1;
elseif r> sum(M1(1:4)) &&  r<= sum(M1(1:5))
    A(5) = 1;
elseif r> sum(M1(1:5)) &&  r<= sum(M1(1:6))  
    A(6) = 1;
elseif r> sum(M1(1:6)) &&  r<= sum(M1(1:7))   
    A(7) = 1;
elseif r> sum(M1(1:7)) &&  r<= sum(M1(1:8))   
    A(8) = 1;
elseif r> sum(M1(1:8)) &&  r<= sum(M1(1:9))
    A(9) = 1;
elseif r> sum(M1(1:9)) &&  r<= sum(M1(1:10))    
    A(10) = 1;
elseif r> sum(M1(1:10)) &&  r<= sum(M1(1:11)) 
    A(11) = 1;
elseif r> sum(M1(1:11)) &&  r<= sum(M1(1:12))   
    A(12) = 1;
elseif r> sum(M1(1:12)) &&  r<= sum(M1(1:13))    
    A(13) = 1;
elseif r> sum(M1(1:13)) &&  r<= sum(M1(1:14))    
    A(14) = 1;
elseif r> sum(M1(1:14)) &&  r<= sum(M1(1:15))    
    A(15) = 1;
else %r> sum(M1(1:15)) &&  r<= sum(M1(1:16))   
    A(16) = 1;
    end  
else
   r = randi([1 sum(M1)],1,1);
if r <= M1(1)
    I(1) = 1;    
elseif r> M1(1) &&  r<= sum(M1(1:2))
    I(2) = 1;  
elseif r> sum(M1(1:2)) &&  r<= sum(M1(1:3))
    I(3) = 1;
elseif r> sum(M1(1:3)) &&  r<= sum(M1(1:4))  
    I(4) = 1;
elseif r> sum(M1(1:4)) &&  r<= sum(M1(1:5))
    I(5) = 1;
elseif r> sum(M1(1:5)) &&  r<= sum(M1(1:6))  
    I(6) = 1;
elseif r> sum(M1(1:6)) &&  r<= sum(M1(1:7))   
    I(7) = 1;
elseif r> sum(M1(1:7)) &&  r<= sum(M1(1:8))   
    I(8) = 1;
elseif r> sum(M1(1:8)) &&  r<= sum(M1(1:9))
    I(9) = 1;
elseif r> sum(M1(1:9)) &&  r<= sum(M1(1:10))    
    I(10) = 1;
elseif r> sum(M1(1:10)) &&  r<= sum(M1(1:11)) 
    I(11) = 1;
elseif r> sum(M1(1:11)) &&  r<= sum(M1(1:12))   
    I(12) = 1;
elseif r> sum(M1(1:12)) &&  r<= sum(M1(1:13))    
    I(13) = 1;
elseif r> sum(M1(1:13)) &&  r<= sum(M1(1:14))    
    I(14) = 1;
elseif r> sum(M1(1:14)) &&  r<= sum(M1(1:15))    
    I(15) = 1;
else %r> sum(M1(1:15)) &&  r<= sum(M1(1:16))   
    I(16) = 1;
end 

end



    
N = N-I-A;
S = N;
L = zeros(group,1);
R = zeros(group,1);
D = zeros(group,1);
SP = zeros(group,1);
SV = zeros(group,1);
NV = aa(mm).*sum(N);

 if NV <= order_vaccaine(1) %used up
     S(order_group(1)) = N(order_group(1)) - NV;
     SV(order_group(1)) = eS*NV;
     SP(order_group(1)) = (1-eS)*NV;
     remain_vaccine = 0; %used up
 else
     S(order_group(1)) = 0;
     SV(order_group(1)) = eS*N(order_group(1));
     SP(order_group(1)) = (1-eS)*N(order_group(1));
     remain_vaccine = NV - N(order_group(1)); %remaining vaccine doses
     count_order =2;
 end

%----------------------------------------------------------------
while remain_vaccine > 0
    if remain_vaccine <= order_vaccaine(count_order) %used up
     S(order_group(count_order)) = N(order_group(count_order)) - remain_vaccine;
     SV(order_group(count_order)) = eS*remain_vaccine;
     SP(order_group(count_order)) = (1-eS)*remain_vaccine;
     remain_vaccine = 0; %used up
 else
     S(order_group(count_order)) = 0;
     SV(order_group(count_order)) = eS*N(order_group(count_order));
     SP(order_group(count_order)) = (1-eS)*N(order_group(count_order));
     remain_vaccine = remain_vaccine - N(order_group(count_order)); %remaining vaccine doses
     count_order = count_order+1;
    end    
end
    

% S = N - (aa(mm).*N);
% SV = eS.*(aa(mm).*N);
% SP = (1-eS).*(aa(mm).*N);
LP = zeros(group,1);
IP = zeros(group,1);
RP = zeros(group,1);
DP = zeros(group,1);
AP = zeros(group,1);

Infected_save(:,1) = sum(I)+sum(A);
C_save(:,1) = sum(I);


pop = zeros(group,12);
pop(:,1) = S;
pop(:,2) = L;
pop(:,3) = I;
pop(:,4) = A;
pop(:,5) = R;
pop(:,6) = SP;
pop(:,7) = LP;
pop(:,8) = IP;
pop(:,9) = AP;
pop(:,10) = RP;
pop(:,11) = D;
pop(:,12) = DP;

S_save(:,1) = sum(S);
L_save(:,1) = sum(L);
I_save(:,1) = sum(I);
A_save(:,1) = sum(A);
R_save(:,1) = sum(R);
D_save(:,1) = sum(D);
C_save(:,1) = sum(I)+sum(A);

SP_save(:,1) = sum(SP);
LP_save(:,1) = sum(LP);
IP_save(:,1) = sum(IP);
AP_save(:,1) = sum(AP);
RP_save(:,1) = sum(RP);
DP_save(:,1) = sum(DP);
CD_save(:,1) = sum(D)+sum(DP);

S_data(:,1) = S;
L_data(:,1) = L;
I_data(:,1) = I;
A_data(:,1) = A;
R_data(:,1) = R;
D_data(:,1) = D;
C_data(:,1) = I+A+IP+AP;


SP_data(:,1) = SP;
LP_data(:,1) = LP;
IP_data(:,1) = IP;
AP_data(:,1) = AP;
RP_data(:,1) = RP;
DP_data(:,1) = DP;
CD_data(:,1) = D+DP;

Infected_data(:,1) = I+A;


while CheckTime < t_final
       
       %% forces of infection
    for j = 1:group 
        
        fI(j) = sum((beta1(j,:).'.*S(j)).*((I)./(N)));
         % S VS A
        fA(j) = sum(qA.*(beta1(j,:).'.*S(j).*(A))./(N));
        % S VS VI
        fIP(j) = sum((1 - eI).*((beta1(j,:).').*S(j).*IP./(N))); 
         % S VS VA
        fAP(j) = sum((1 - eI).*qA.*((beta1(j,:).').*S(j).*(AP./(N))));
        % SP VS I      
        fI2(j) = sum((beta1(j,:).').*SP(j).*(I)./(N));  
          % SP VS A      
        fA2(j) = sum(qA.*(beta1(j,:).').*SP(j).*(A)./(N));
        % SP VS IP
        fIP2(j) = sum((1 - eI).*((beta1(j,:).').*SP(j).*IP./(N)));
        % SP VS AP       
        fAP2(j) = sum((1 - eI).*qA.*((beta1(j,:).').*SP(j).*(AP./(N))));
    end
 
 % ---------- Rate of Each Event---------%
 % [S L I A R SP LP IP AP RP D DP]
% -------------- spread disease form I to S ----------%
Rate1(1) = fI(1);  Change1(1,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate1(2) = fI(2);  Change1(2,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g2
Rate1(3) = fI(3);  Change1(3,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g3
Rate1(4) = fI(4);  Change1(4,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g4
Rate1(5) = fI(5);  Change1(5,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g5
Rate1(6) = fI(6);  Change1(6,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g6
Rate1(7) = fI(7);  Change1(7,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g7
Rate1(8) = fI(8);  Change1(8,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g8
Rate1(9) = fI(9);  Change1(9,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g9
Rate1(10) = fI(10);  Change1(10,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g10
Rate1(11) = fI(11);  Change1(11,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g11
Rate1(12) = fI(12);  Change1(12,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g12
Rate1(13) = fI(13);  Change1(13,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g13
Rate1(14) = fI(14);  Change1(14,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g14
Rate1(15) = fI(15);  Change1(15,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g15
Rate1(16) = fI(16);  Change1(16,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g16
% -------------- spread disease form IP to S ----------%
Rate2(1) = fIP(1);  Change2(1,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate2(2) = fIP(2);  Change2(2,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g2
Rate2(3) = fIP(3);  Change2(3,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g3
Rate2(4) = fIP(4);  Change2(4,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g4
Rate2(5) = fIP(5);  Change2(5,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g5
Rate2(6) = fIP(6);  Change2(6,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g6
Rate2(7) = fIP(7);  Change2(7,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g7
Rate2(8) = fIP(8);  Change2(8,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g8
Rate2(9) = fIP(9);  Change2(9,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g9
Rate2(10) = fIP(10);  Change2(10,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g10
Rate2(11) = fIP(11);  Change2(11,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g11
Rate2(12) = fIP(12);  Change2(12,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g12
Rate2(13) = fIP(13);  Change2(13,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g13
Rate2(14) = fIP(14);  Change2(14,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g14
Rate2(15) = fIP(15);  Change2(15,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g15
Rate2(16) = fIP(16);  Change2(16,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g16
% -------------- spread disease form I to SP ----------%
Rate3(1) = fI2(1);  Change3(1,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate3(2) = fI2(2);  Change3(2,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0];%g2
Rate3(3) = fI2(3);  Change3(3,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g3
Rate3(4) = fI2(4);  Change3(4,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g4
Rate3(5) = fI2(5);  Change3(5,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g5
Rate3(6) = fI2(6);  Change3(6,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g6
Rate3(7) = fI2(7);  Change3(7,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g7
Rate3(8) = fI2(8);  Change3(8,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g8
Rate3(9) = fI2(9);  Change3(9,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g9
Rate3(10) = fI2(10);  Change3(10,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g10
Rate3(11) = fI2(11);  Change3(11,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g11
Rate3(12) = fI2(12);  Change3(12,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g12
Rate3(13) = fI2(13);  Change3(13,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g13
Rate3(14) = fI2(14);  Change3(14,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g14
Rate3(15) = fI2(15);  Change3(15,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g15
Rate3(16) = fI2(16);  Change3(16,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g16
% -------------- spread disease form IP to SP ----------%
Rate4(1) = fIP2(1);  Change4(1,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate4(2) = fIP2(2);  Change4(2,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g2
Rate4(3) = fIP2(3);  Change4(3,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g3
Rate4(4) = fIP2(4);  Change4(4,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g4
Rate4(5) = fIP2(5);  Change4(5,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g5
Rate4(6) = fIP2(6);  Change4(6,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g6
Rate4(7) = fIP2(7);  Change4(7,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g7
Rate4(8) = fIP2(8);  Change4(8,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0 ]; %g8
Rate4(9) = fIP2(9);  Change4(9,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g9
Rate4(10) = fIP2(10);  Change4(10,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g10
Rate4(11) = fIP2(11);  Change4(11,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g11
Rate4(12) = fIP2(12);  Change4(12,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g12
Rate4(13) = fIP2(13);  Change4(13,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g13
Rate4(14) = fIP2(14);  Change4(14,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g14
Rate4(15) = fIP2(15);  Change4(15,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g15
Rate4(16) = fIP2(16);  Change4(16,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g16
% -------------- L ----> I ----------%
Rate5(1) = FI.*sigma*L(1);  Change5(1,:)= [0 -1 +1 0 0 0 0 0 0 0 0 0]; %g1
Rate5(2) = FI.*sigma*L(2);  Change5(2,:)= [0 -1 +1 0 0 0 0 0 0 0 0 0]; %g2
Rate5(3) = FI.*sigma*L(3);  Change5(3,:)= [0 -1 +1 0 0 0 0 0 0 0 0 0]; %g3
Rate5(4) = FI.*sigma*L(4);  Change5(4,:)= [0 -1 +1 0 0 0 0 0 0 0 0 0]; %g4
Rate5(5) = FI.*sigma*L(5);  Change5(5,:)= [0 -1 +1 0 0 0 0 0 0 0 0 0]; %g5
Rate5(6) = FI.*sigma*L(6);  Change5(6,:)= [0 -1 +1 0 0 0 0 0 0 0 0 0]; %g6
Rate5(7) = FI.*sigma*L(7);  Change5(7,:)= [0 -1 +1 0 0 0 0 0 0 0 0 0]; %g7
Rate5(8) = FI.*sigma*L(8);  Change5(8,:)= [0 -1 +1 0 0 0 0 0 0 0 0 0]; %g8
Rate5(9) = FI.*sigma*L(9);  Change5(9,:)= [0 -1 +1 0 0 0 0 0 0 0 0 0]; %g9
Rate5(10) = FI.*sigma*L(10);  Change5(10,:)= [0 -1 +1 0 0 0 0 0 0 0 0 0]; %g10
Rate5(11) = FI.*sigma*L(11);  Change5(11,:)= [0 -1 +1 0 0 0 0 0 0 0 0 0]; %g11
Rate5(12) = FI.*sigma*L(12);  Change5(12,:)= [0 -1 +1 0 0 0 0 0 0 0 0 0]; %g12
Rate5(13) = FI.*sigma*L(13);  Change5(13,:)= [0 -1 +1 0 0 0 0 0 0 0 0 0]; %g13
Rate5(14) = FI.*sigma*L(14);  Change5(14,:)= [0 -1 +1 0 0 0 0 0 0 0 0 0]; %g14
Rate5(15) = FI.*sigma*L(15);  Change5(15,:)= [0 -1 +1 0 0 0 0 0 0 0 0 0]; %g15
Rate5(16) = FI.*sigma*L(16);  Change5(16,:)= [0 -1 +1 0 0 0 0 0 0 0 0 0]; %g16
% -------------- LP ----> IP ----------%
Rate6(1) = (1-FAP).*sigma*LP(1);  Change6(1,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g1
Rate6(2) = (1-FAP).*sigma*LP(2);  Change6(2,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g2
Rate6(3) = (1-FAP).*sigma*LP(3);  Change6(3,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g3
Rate6(4) = (1-FAP).*sigma*LP(4);  Change6(4,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g4
Rate6(5) = (1-FAP).*sigma*LP(5);  Change6(5,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g5
Rate6(6) = (1-FAP).*sigma*LP(6);  Change6(6,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g6
Rate6(7) = (1-FAP).*sigma*LP(7);  Change6(7,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g7
Rate6(8) = (1-FAP).*sigma*LP(8);  Change6(8,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g8
Rate6(9) = (1-FAP).*sigma*LP(9);  Change6(9,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g9
Rate6(10) = (1-FAP).*sigma*LP(10);  Change6(10,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g10
Rate6(11) = (1-FAP).*sigma*LP(11);  Change6(11,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g11
Rate6(12) = (1-FAP).*sigma*LP(12);  Change6(12,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g12
Rate6(13) = (1-FAP).*sigma*LP(13);  Change6(13,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g13
Rate6(14) = (1-FAP).*sigma*LP(14);  Change6(14,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g14
Rate6(15) = (1-FAP).*sigma*LP(15);  Change6(15,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g15
Rate6(16) = (1-FAP).*sigma*LP(16);  Change6(16,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g16
% -------------- I ----> R ----------%
Rate7(1) = kIR(1)*I(1);  Change7(1,:)= [0 0 -1 0 +1 0 0 0 0 0 0 0]; %g1
Rate7(2) = kIR(2)*I(2);  Change7(2,:)= [0 0 -1 0 +1 0 0 0 0 0 0 0]; %g2
Rate7(3) = kIR(3)*I(3);  Change7(3,:)= [0 0 -1 0 +1 0 0 0 0 0 0 0]; %g3
Rate7(4) = kIR(4)*I(4);  Change7(4,:)= [0 0 -1 0 +1 0 0 0 0 0 0 0]; %g4
Rate7(5) = kIR(5)*I(5);  Change7(5,:)= [0 0 -1 0 +1 0 0 0 0 0 0 0]; %g5
Rate7(6) = kIR(6)*I(6);  Change7(6,:)= [0 0 -1 0 +1 0 0 0 0 0 0 0]; %g6
Rate7(7) = kIR(7)*I(7);  Change7(7,:)= [0 0 -1 0 +1 0 0 0 0 0 0 0]; %g7
Rate7(8) = kIR(8)*I(8);  Change7(8,:)= [0 0 -1 0 +1 0 0 0 0 0 0 0]; %g8
Rate7(9) = kIR(9)*I(9);  Change7(9,:)= [0 0 -1 0 +1 0 0 0 0 0 0 0]; %g9
Rate7(10) = kIR(10)*I(10);  Change7(10,:)= [0 0 -1 0 +1 0 0 0 0 0 0 0]; %g10
Rate7(11) = kIR(11)*I(11);  Change7(11,:)= [0 0 -1 0 +1 0 0 0 0 0 0 0]; %g11
Rate7(12) = kIR(12)*I(12);  Change7(12,:)= [0 0 -1 0 +1 0 0 0 0 0 0 0]; %g12
Rate7(13) = kIR(13)*I(13);  Change7(13,:)= [0 0 -1 0 +1 0 0 0 0 0 0 0]; %g13
Rate7(14) = kIR(14)*I(14);  Change7(14,:)= [0 0 -1 0 +1 0 0 0 0 0 0 0]; %g14
Rate7(15) = kIR(15)*I(15);  Change7(15,:)= [0 0 -1 0 +1 0 0 0 0 0 0 0]; %g15
Rate7(16) = kIR(16)*I(16);  Change7(16,:)= [0 0 -1 0 +1 0 0 0 0 0 0 0]; %g16
% -------------- IP ----> RP ----------%
Rate8(1) = kIPRP(1)*IP(1);  Change8(1,:)= [0 0 0 0 0 0 0 -1 0 +1 0 0]; %g1
Rate8(2) =  kIPRP(2)*IP(2);  Change8(2,:)= [0 0 0 0 0 0 0 -1 0 +1 0 0]; %g2
Rate8(3) =  kIPRP(3)*IP(3);  Change8(3,:)= [0 0 0 0 0 0 0 -1 0 +1 0 0]; %g3
Rate8(4) =  kIPRP(4)*IP(4);  Change8(4,:)= [0 0 0 0 0 0 0 -1 0 +1 0 0]; %g4
Rate8(5) =  kIPRP(5)*IP(5);  Change8(5,:)= [0 0 0 0 0 0 0 -1 0 +1 0 0]; %g5
Rate8(6) =  kIPRP(6)*IP(6);  Change8(6,:)= [0 0 0 0 0 0 0 -1 0 +1 0 0]; %g6
Rate8(7) =  kIPRP(7)*IP(7);  Change8(7,:)= [0 0 0 0 0 0 0 -1 0 +1 0 0]; %g7
Rate8(8) =  kIPRP(8)*IP(8);  Change8(8,:)= [0 0 0 0 0 0 0 -1 0 +1 0 0]; %g8
Rate8(9) =  kIPRP(9)*IP(9);  Change8(9,:)= [0 0 0 0 0 0 0 -1 0 +1 0 0]; %g9
Rate8(10) =  kIPRP(10)*IP(10);  Change8(10,:)= [0 0 0 0 0 0 0 -1 0 +1 0 0]; %g10
Rate8(11) =  kIPRP(11)*IP(11);  Change8(11,:)= [0 0 0 0 0 0 0 -1 0 +1 0 0]; %g11
Rate8(12) =  kIPRP(12)*IP(12);  Change8(12,:)= [0 0 0 0 0 0 0 -1 0 +1 0 0]; %g12
Rate8(13) =  kIPRP(13)*IP(13);  Change8(13,:)= [0 0 0 0 0 0 0 -1 0 +1 0 0]; %g13
Rate8(14) =  kIPRP(14)*IP(14);  Change8(14,:)= [0 0 0 0 0 0 0 -1 0 +1 0 0]; %g14
Rate8(15) =  kIPRP(15)*IP(15);  Change8(15,:)= [0 0 0 0 0 0 0 -1 0 +1 0 0]; %g15
Rate8(16) =  kIPRP(16)*IP(16);  Change8(16,:)= [0 0 0 0 0 0 0 -1 0 +1 0 0]; %g16
% -------------- spread disease form A to S ----------%
Rate9(1) = fA(1);  Change9(1,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate9(2) = fA(2);  Change9(2,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate9(3) = fA(3);  Change9(3,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate9(4) = fA(4);  Change9(4,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate9(5) = fA(5);  Change9(5,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate9(6) = fA(6);  Change9(6,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate9(7) = fA(7);  Change9(7,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate9(8) = fA(8);  Change9(8,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate9(9) = fA(9);  Change9(9,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate9(10) = fA(10);  Change9(10,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate9(11) = fA(11);  Change9(11,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate9(12) = fA(12);  Change9(12,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate9(13) = fA(13);  Change9(13,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate9(14) = fA(14);  Change9(14,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate9(15) = fA(15);  Change9(15,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate9(16) = fA(16);  Change9(16,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
% -------------- spread disease form AP to S ----------%
Rate10(1) = fAP(1);  Change10(1,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate10(2) = fAP(2);  Change10(2,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate10(3) = fAP(3);  Change10(3,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate10(4) = fAP(4);  Change10(4,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate10(5) = fAP(5);  Change10(5,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate10(6) = fAP(6);  Change10(6,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate10(7) = fAP(7);  Change10(7,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate10(8) = fAP(8);  Change10(8,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate10(9) = fAP(9);  Change10(9,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate10(10) = fAP(10);  Change10(10,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate10(11) = fAP(11);  Change10(11,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate10(12) = fAP(12);  Change10(12,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate10(13) = fAP(13);  Change10(13,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate10(14) = fAP(14);  Change10(14,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate10(15) = fAP(15);  Change10(15,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
Rate10(16) = fAP(16);  Change10(16,:)= [-1 +1 0 0 0 0 0 0 0 0 0 0]; %g1
% -------------- spread disease form A to SP ----------%
Rate11(1) = fA2(1);  Change11(1,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate11(2) = fA2(2);  Change11(2,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate11(3) = fA2(3);  Change11(3,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate11(4) = fA2(4);  Change11(4,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate11(5) = fA2(5);  Change11(5,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate11(6) = fA2(6);  Change11(6,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate11(7) = fA2(7);  Change11(7,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate11(8) = fA2(8);  Change11(8,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate11(9) = fA2(9);  Change11(9,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate11(10) = fA2(10);  Change11(10,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate11(11) = fA2(11);  Change11(11,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate11(12) = fA2(12);  Change11(12,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate11(13) = fA2(13);  Change11(13,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate11(14) = fA2(14);  Change11(14,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate11(15) = fA2(15);  Change11(15,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate11(16) = fA2(16);  Change11(16,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1

% -------------- spread disease form AP to SP ----------%
Rate12(1) = fAP2(1);  Change12(1,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate12(2) = fAP2(2);  Change12(2,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate12(3) = fAP2(3);  Change12(3,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate12(4) = fAP2(4);  Change12(4,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate12(5) = fAP2(5);  Change12(5,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate12(6) = fAP2(6);  Change12(6,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate12(7) = fAP2(7);  Change12(7,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate12(8) = fAP2(8);  Change12(8,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate12(9) = fAP2(9);  Change12(9,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate12(10) = fAP2(10);  Change12(10,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate12(11) = fAP2(11);  Change12(11,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate12(12) = fAP2(12);  Change12(12,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate12(13) = fAP2(13);  Change12(13,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate12(14) = fAP2(14);  Change12(14,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate12(15) = fAP2(15);  Change12(15,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
Rate12(16) = fAP2(16);  Change12(16,:)= [0 0 0 0 0 -1 +1 0 0 0 0 0]; %g1
% -------------- L ----> A ----------%
Rate13(1) = FA.*sigma*L(1);  Change13(1,:)= [0 -1 0 +1 0 0 0 0 0 0 0 0]; %g1
Rate13(2) = FA.*sigma*L(2);  Change13(2,:)= [0 -1 0 +1 0 0 0 0 0 0 0 0]; %g1
Rate13(3) = FA.*sigma*L(3);  Change13(3,:)= [0 -1 0 +1 0 0 0 0 0 0 0 0]; %g1
Rate13(4) = FA.*sigma*L(4);  Change13(4,:)= [0 -1 0 +1 0 0 0 0 0 0 0 0]; %g1
Rate13(5) = FA.*sigma*L(5);  Change13(5,:)= [0 -1 0 +1 0 0 0 0 0 0 0 0]; %g1
Rate13(6) = FA.*sigma*L(6);  Change13(6,:)= [0 -1 0 +1 0 0 0 0 0 0 0 0]; %g1
Rate13(7) = FA.*sigma*L(7);  Change13(7,:)= [0 -1 0 +1 0 0 0 0 0 0 0 0]; %g1
Rate13(8) = FA.*sigma*L(8);  Change13(8,:)= [0 -1 0 +1 0 0 0 0 0 0 0 0]; %g1
Rate13(9) = FA.*sigma*L(9);  Change13(9,:)= [0 -1 0 +1 0 0 0 0 0 0 0 0]; %g1
Rate13(10) = FA.*sigma*L(10);  Change13(10,:)= [0 -1 0 +1 0 0 0 0 0 0 0 0]; %g1
Rate13(11) = FA.*sigma*L(11);  Change13(11,:)= [0 -1 0 +1 0 0 0 0 0 0 0 0]; %g1
Rate13(12) = FA.*sigma*L(12);  Change13(12,:)= [0 -1 0 +1 0 0 0 0 0 0 0 0]; %g1
Rate13(13) = FA.*sigma*L(13);  Change13(13,:)= [0 -1 0 +1 0 0 0 0 0 0 0 0]; %g1
Rate13(14) = FA.*sigma*L(14);  Change13(14,:)= [0 -1 0 +1 0 0 0 0 0 0 0 0]; %g1
Rate13(15) = FA.*sigma*L(15);  Change13(15,:)= [0 -1 0 +1 0 0 0 0 0 0 0 0]; %g1
Rate13(16) = FA.*sigma*L(16);  Change13(16,:)= [0 -1 0 +1 0 0 0 0 0 0 0 0]; %g1

% -------------- LP ----> AP ----------%
Rate14(1) = FAP.*sigma*LP(1);  Change14(1,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g1
Rate14(2) = FAP.*sigma*LP(2);  Change14(2,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g1
Rate14(3) = FAP.*sigma*LP(3);  Change14(3,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g1
Rate14(4) = FAP.*sigma*LP(4);  Change14(4,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g1
Rate14(5) = FAP.*sigma*LP(5);  Change14(5,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g1
Rate14(6) = FAP.*sigma*LP(6);  Change14(6,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g1
Rate14(7) = FAP.*sigma*LP(7);  Change14(7,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g1
Rate14(8) = FAP.*sigma*LP(8);  Change14(8,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g1
Rate14(9) = FAP.*sigma*LP(9);  Change14(9,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g1
Rate14(10) = FAP.*sigma*LP(10);  Change14(10,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g1
Rate14(11) = FAP.*sigma*LP(11);  Change14(11,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g1
Rate14(12) = FAP.*sigma*LP(12);  Change14(12,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g1
Rate14(13) = FAP.*sigma*LP(13);  Change14(13,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g1
Rate14(14) = FAP.*sigma*LP(14);  Change14(14,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g1
Rate14(15) = FAP.*sigma*LP(15);  Change14(15,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g1
Rate14(16) = FAP.*sigma*LP(16);  Change14(16,:)= [0 0 0 0 0 0 -1 +1 0 0 0 0]; %g1
% -------------- A ----> R ----------%
Rate15(1) = gamma*A(1);  Change15(1,:)= [0 0 0 -1 +1 0 0 0 0 0 0 0]; %g1
Rate15(2) = gamma*A(2);  Change15(2,:)= [0 0 0 -1 +1 0 0 0 0 0 0 0]; %g1
Rate15(3) = gamma*A(3);  Change15(3,:)= [0 0 0 -1 +1 0 0 0 0 0 0 0]; %g1
Rate15(4) = gamma*A(4);  Change15(4,:)= [0 0 0 -1 +1 0 0 0 0 0 0 0]; %g1
Rate15(5) = gamma*A(5);  Change15(5,:)= [0 0 0 -1 +1 0 0 0 0 0 0 0]; %g1
Rate15(6) = gamma*A(6);  Change15(6,:)= [0 0 0 -1 +1 0 0 0 0 0 0 0]; %g1
Rate15(7) = gamma*A(7);  Change15(7,:)= [0 0 0 -1 +1 0 0 0 0 0 0 0]; %g1
Rate15(8) = gamma*A(8);  Change15(8,:)= [0 0 0 -1 +1 0 0 0 0 0 0 0]; %g1
Rate15(9) = gamma*A(9);  Change15(9,:)= [0 0 0 -1 +1 0 0 0 0 0 0 0]; %g1
Rate15(10) = gamma*A(10);  Change15(10,:)= [0 0 0 -1 +1 0 0 0 0 0 0 0]; %g1
Rate15(11) = gamma*A(11);  Change15(11,:)= [0 0 0 -1 +1 0 0 0 0 0 0 0]; %g1
Rate15(12) = gamma*A(12);  Change15(12,:)= [0 0 0 -1 +1 0 0 0 0 0 0 0]; %g1
Rate15(13) = gamma*A(13);  Change15(13,:)= [0 0 0 -1 +1 0 0 0 0 0 0 0]; %g1
Rate15(14) = gamma*A(14);  Change15(14,:)= [0 0 0 -1 +1 0 0 0 0 0 0 0]; %g1
Rate15(15) = gamma*A(15);  Change15(15,:)= [0 0 0 -1 +1 0 0 0 0 0 0 0]; %g1
Rate15(16) = gamma*A(16);  Change15(16,:)= [0 0 0 -1 +1 0 0 0 0 0 0 0]; %g1
% -------------- AP ----> RP ----------%
Rate16(1) = gamma*AP(1);  Change16(1,:)= [0 0 0 0 0 0 0 0 -1 +1 0 0]; %g1
Rate16(2) = gamma*AP(2);  Change16(2,:)= [0 0 0 0 0 0 0 0 -1 +1 0 0]; %g1
Rate16(3) = gamma*AP(3);  Change16(3,:)= [0 0 0 0 0 0 0 0 -1 +1 0 0]; %g1
Rate16(4) = gamma*AP(4);  Change16(4,:)= [0 0 0 0 0 0 0 0 -1 +1 0 0]; %g1
Rate16(5) = gamma*AP(5);  Change16(5,:)= [0 0 0 0 0 0 0 0 -1 +1 0 0]; %g1
Rate16(6) = gamma*AP(6);  Change16(6,:)= [0 0 0 0 0 0 0 0 -1 +1 0 0]; %g1
Rate16(7) = gamma*AP(7);  Change16(7,:)= [0 0 0 0 0 0 0 0 -1 +1 0 0]; %g1
Rate16(8) = gamma*AP(8);  Change16(8,:)= [0 0 0 0 0 0 0 0 -1 +1 0 0]; %g1
Rate16(9) = gamma*AP(9);  Change16(9,:)= [0 0 0 0 0 0 0 0 -1 +1 0 0]; %g1
Rate16(10) = gamma*AP(10);  Change16(10,:)= [0 0 0 0 0 0 0 0 -1 +1 0 0]; %g1
Rate16(11) = gamma*AP(11);  Change16(11,:)= [0 0 0 0 0 0 0 0 -1 +1 0 0]; %g1
Rate16(12) = gamma*AP(12);  Change16(12,:)= [0 0 0 0 0 0 0 0 -1 +1 0 0]; %g1
Rate16(13) = gamma*AP(13);  Change16(13,:)= [0 0 0 0 0 0 0 0 -1 +1 0 0]; %g1
Rate16(14) = gamma*AP(14);  Change16(14,:)= [0 0 0 0 0 0 0 0 -1 +1 0 0]; %g1
Rate16(15) = gamma*AP(15);  Change16(15,:)= [0 0 0 0 0 0 0 0 -1 +1 0 0]; %g1
Rate16(16) = gamma*AP(16);  Change16(16,:)= [0 0 0 0 0 0 0 0 -1 +1 0 0]; %g1
 % [S L I A R SP LP IP AP RP D DP]
% -------------- Death form I to D ----------%
Rate17(1) = kID(1)*I(1);  Change17(1,:)= [0 0 -1 0 0 0 0 0 0 0 +1 0]; %g1
Rate17(2) = kID(2)*I(2);  Change17(2,:)= [0 0 -1 0 0 0 0 0 0 0 +1 0]; %g2
Rate17(3) = kID(3)*I(3);  Change17(3,:)= [0 0 -1 0 0 0 0 0 0 0 +1 0]; %g3
Rate17(4) = kID(4)*I(4);  Change17(4,:)= [0 0 -1 0 0 0 0 0 0 0 +1 0]; %g4
Rate17(5) = kID(5)*I(5);  Change17(5,:)= [0 0 -1 0 0 0 0 0 0 0 +1 0]; %g5
Rate17(6) = kID(6)*I(6);  Change17(6,:)= [0 0 -1 0 0 0 0 0 0 0 +1 0]; %g6
Rate17(7) = kID(7)*I(7);  Change17(7,:)= [0 0 -1 0 0 0 0 0 0 0 +1 0]; %g7
Rate17(8) = kID(8)*I(8);  Change17(8,:)= [0 0 -1 0 0 0 0 0 0 0 +1 0]; %g8
Rate17(9) = kID(9)*I(9);  Change17(9,:)= [0 0 -1 0 0 0 0 0 0 0 +1 0]; %g9
Rate17(10) = kID(10)*I(10);  Change17(10,:)= [0 0 -1 0 0 0 0 0 0 0 +1 0]; %g10
Rate17(11) = kID(11)*I(11);  Change17(11,:)= [0 0 -1 0 0 0 0 0 0 0 +1 0]; %g11
Rate17(12) = kID(12)*I(12);  Change17(12,:)= [0 0 -1 0 0 0 0 0 0 0 +1 0]; %g12
Rate17(13) = kID(13)*I(13);  Change17(13,:)= [0 0 -1 0 0 0 0 0 0 0 +1 0]; %g13
Rate17(14) = kID(14)*I(14);  Change17(14,:)= [0 0 -1 0 0 0 0 0 0 0 +1 0]; %g14
Rate17(15) = kID(15)*I(15);  Change17(15,:)= [0 0 -1 0 0 0 0 0 0 0 +1 0]; %g15
Rate17(16) = kID(16)*I(16);  Change17(16,:)= [0 0 -1 0 0 0 0 0 0 0 +1 0]; %g16

% -------------- IP ----> DP ----------%
Rate18(1) = kIPDP(1)*IP(1);  Change18(1,:)= [0 0 0 0 0 0 0 -1 0 0 0 +1]; %g1
Rate18(2) =  kIPDP(2)*IP(2);  Change18(2,:)= [0 0 0 0 0 0 0 -1 0 0 0 +1]; %g2
Rate18(3) =  kIPDP(3)*IP(3);  Change18(3,:)= [0 0 0 0 0 0 0 -1 0 0 0 +1]; %g3
Rate18(4) =  kIPDP(4)*IP(4);  Change18(4,:)= [0 0 0 0 0 0 0 -1 0 0 0 +1]; %g4
Rate18(5) =  kIPDP(5)*IP(5);  Change18(5,:)= [0 0 0 0 0 0 0 -1 0 0 0 +1]; %g5
Rate18(6) =  kIPDP(6)*IP(6);  Change18(6,:)= [0 0 0 0 0 0 0 -1 0 0 0 +1]; %g6
Rate18(7) =  kIPDP(7)*IP(7);  Change18(7,:)= [0 0 0 0 0 0 0 -1 0 0 0 +1]; %g7
Rate18(8) =  kIPDP(8)*IP(8);  Change18(8,:)= [0 0 0 0 0 0 0 -1 0 0 0 +1]; %g8
Rate18(9) =  kIPDP(9)*IP(9);  Change18(9,:)= [0 0 0 0 0 0 0 -1 0 0 0 +1]; %g9
Rate18(10) =  kIPDP(10)*IP(10);  Change18(10,:)= [0 0 0 0 0 0 0 -1 0 0 0 +1]; %g10
Rate18(11) =  kIPDP(11)*IP(11);  Change18(11,:)= [0 0 0 0 0 0 0 -1 0 0 0 +1]; %g11
Rate18(12) =  kIPDP(12)*IP(12);  Change18(12,:)= [0 0 0 0 0 0 0 -1 0 0 0 +1]; %g12
Rate18(13) =  kIPDP(13)*IP(13);  Change18(13,:)= [0 0 0 0 0 0 0 -1 0 0 0 +1];%g13
Rate18(14) =  kIPDP(14)*IP(14);  Change18(14,:)= [0 0 0 0 0 0 0 -1 0 0 0 +1]; %g14
Rate18(15) =  kIPDP(15)*IP(15);  Change18(15,:)= [0 0 0 0 0 0 0 -1 0 0 0 +1]; %g15
Rate18(16) =  kIPDP(16)*IP(16);  Change18(16,:)= [0 0 0 0 0 0 0 -1 0 0 0 +1]; %g16
 

for kk = 1:16
parameter = [S(kk) L(kk) I(kk) A(kk) R(kk) SP(kk) LP(kk) IP(kk) AP(kk) RP(kk) D(kk) DP(kk) fI(kk) fI2(kk)...
              fIP(kk) fIP2(kk) fA(kk) fA2(kk) fAP(kk) fAP2(kk) sigma gamma kID(kk) kIR(kk) kIPDP(kk) kIPRP(kk) dt];
i = kk;
          

    Num1 = poissrnd(Rate1(i)*dt);
    Num2 = poissrnd(Rate2(i)*dt);
    Num3 = poissrnd(Rate3(i)*dt);
    Num4 = poissrnd(Rate4(i)*dt);
    Num5 = poissrnd(Rate5(i)*dt);
    Num6 = poissrnd(Rate6(i)*dt);
    Num7 = poissrnd(Rate7(i)*dt);
    Num8 = poissrnd(Rate8(i)*dt);
    Num9 = poissrnd(Rate9(i)*dt);
    Num10 = poissrnd(Rate10(i)*dt);
    Num11 = poissrnd(Rate11(i)*dt);
    Num12 = poissrnd(Rate12(i)*dt);
    Num13 = poissrnd(Rate13(i)*dt);
    Num14 = poissrnd(Rate14(i)*dt);
    Num15 = poissrnd(Rate15(i)*dt);
    Num16 = poissrnd(Rate16(i)*dt);
    Num17 = poissrnd(Rate17(i)*dt);
    Num18 = poissrnd(Rate18(i)*dt);    
    
    
    % Make sure things don't go negative
    Use1 = min([Num1 parameter(find(Change1(i,:)<0))]);
    Use2 = min([Num2 parameter(find(Change2(i,:)<0))]);
    Use3 = min([Num3 parameter(find(Change3(i,:)<0))]);
    Use4 = min([Num4 parameter(find(Change4(i,:)<0))]);
    Use5 = min([Num5 parameter(find(Change5(i,:)<0))]);
    Use6 = min([Num6 parameter(find(Change6(i,:)<0))]);
    Use7 = min([Num7 parameter(find(Change7(i,:)<0))]);
    Use8 = min([Num8 parameter(find(Change8(i,:)<0))]);
    Use9 = min([Num9 parameter(find(Change9(i,:)<0))]);
    Use10 = min([Num10 parameter(find(Change10(i,:)<0))]);
    Use11 = min([Num11 parameter(find(Change11(i,:)<0))]);
    Use12 = min([Num12 parameter(find(Change12(i,:)<0))]);
    Use13 = min([Num13 parameter(find(Change13(i,:)<0))]);
    Use14 = min([Num14 parameter(find(Change14(i,:)<0))]);
    Use15 = min([Num15 parameter(find(Change15(i,:)<0))]);
    Use16 = min([Num16 parameter(find(Change16(i,:)<0))]);
    Use17 = min([Num17 parameter(find(Change17(i,:)<0))]);
    Use18 = min([Num18 parameter(find(Change18(i,:)<0))]);
    
   pop(i,:) = pop(i,:) + Change1(i,:)*Use1 + Change2(i,:)*Use2 + Change3(i,:)*Use3...
          + Change4(i,:)*Use4 + Change5(i,:)*Use5 + Change6(i,:)*Use6 ...
          + Change7(i,:)*Use7 + Change8(i,:)*Use8  + Change9(i,:)*Use9...
          + Change10(i,:)*Use10  + Change11(i,:)*Use11  + Change12(i,:)*Use12...
          + Change13(i,:)*Use13  + Change14(i,:)*Use14  + Change15(i,:)*Use15...
          + Change16(i,:)*Use16 + Change17(i,:)*Use17 + Change18(i,:)*Use18;


end
% Update Parameters
S_save(l,step+1) = sum(pop(:,1));
L_save(l,step+1) = sum(pop(:,2));
I_save(l,step+1) = sum(pop(:,3));
A_save(l,step+1) = sum(pop(:,4));
R_save(l,step+1) = sum(pop(:,5));
SP_save(l,step+1) = sum(pop(:,6));
LP_save(l,step+1) = sum(pop(:,7));
IP_save(l,step+1) = sum(pop(:,8));
AP_save(l,step+1) = sum(pop(:,9));
RP_save(l,step+1) = sum(pop(:,10));
D_save(l,step+1) = sum(pop(:,11));
DP_save(l,step+1) = sum(pop(:,12));

Infected_save(l,step+1) = sum(pop(:,2)) + sum(pop(:,3)) + sum(pop(:,4)) + sum(pop(:,7)) + sum(pop(:,8))+ sum(pop(:,9));
CD_save(l,step+1) = sum(pop(:,11)) + sum(pop(:,12));
C_save(l,step+1) = Infected_save(l,step+1) + sum(pop(:,5)) + sum(pop(:,10)) + CD_save(l,step+1);

% 
S_data(:,step+1) = S(:);
L_data(:,step+1) = L(:);
I_data(:,step+1) = I(:);
A_data(:,step+1) = A(:);
R_data(:,step+1) = R(:);
D_data(:,step+1) = D(:);
        
SP_data(:,step+1) = SP(:);
LP_data(:,step+1) = LP(:);        
IP_data(:,step+1) = IP(:); 
AP_data(:,step+1) = AP(:); 
RP_data(:,step+1) = RP(:);   
DP_data(:,step+1) = DP(:); 
%         
Infected_data(:,step+1) = L(:)+LP(:)+I(:)+IP(:)+A(:)+AP(:);
C_data(:,step+1) = L(:)+LP(:)+I(:)+IP(:)+A(:)+AP(:)+R(:)+RP(:)+D(:)+DP(:);
CD_data(:,step+1) = D(:)+DP(:);
% 
% N_data(l,step+1)= sum(sum(pop)); 
CheckTime = CheckTime + dt;
if Infected_save(l,step+1) == 0
    tzero(l,mm) = CheckTime;
    Czero(l,mm) = C_save(l,step+1);
    CDzero(l,mm) = CD_save(l,step+1);
%     tcount = tcount+1;
    break;
end

step = step + 1;

% TimeStochastic(stepp+step) = CheckTime;
S = pop(:,1);
L = pop(:,2);
I = pop(:,3);
A = pop(:,4);
R = pop(:,5);
SP = pop(:,6);
LP = pop(:,7);
IP = pop(:,8);
AP = pop(:,9);
RP = pop(:,10);
D = pop(:,11);
DP = pop(:,12);

end

%% Save
%  filename = ['R0',num2str(R0),'VC',num2str(aa(mm)),'REP',num2str(l) '.mat'];
% save(filename);
l
end
mm
%% Save
% Infected_mean = mean(Infected_save);
C_mean = mean(C_save);
S_mean = mean(S_save);
L_mean = mean(L_save);
I_mean = mean(I_save);
A_mean = mean(A_save);
R_mean = mean(R_save);
D_mean = mean(D_save);

SP_mean = mean(SP_save);
LP_mean = mean(LP_save);
IP_mean = mean(IP_save);
AP_mean = mean(AP_save);
RP_mean = mean(RP_save);
DP_mean = mean(DP_save);
CD_mean = mean(CD_save);

S_error = std(S_save)/sqrt(Repetition);
L_error = std(L_save)/sqrt(Repetition);
I_error = std(I_save)/sqrt(Repetition);
A_error = std(A_save)/sqrt(Repetition);
R_error = std(R_save)/sqrt(Repetition);
D_error = std(D_save)/sqrt(Repetition);
SP_error = std(SP_save)/sqrt(Repetition);
LP_error = std(LP_save)/sqrt(Repetition);
IP_error = std(IP_save)/sqrt(Repetition);
AP_error = std(AP_save)/sqrt(Repetition);
RP_error = std(RP_save)/sqrt(Repetition);
DP_error = std(DP_save)/sqrt(Repetition);
C_error = std(C_save)/sqrt(Repetition);
CD_error = std(CD_save)/sqrt(Repetition);
% Infected_error = std(Infected_save)/sqrt(Repetition);

filename = ['IorderBatch1month1','R0',num2str(R0),'VC',num2str(aa(mm)),'REP',num2str(l) '.mat'];
save(filename);
end
toc