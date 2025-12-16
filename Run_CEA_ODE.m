%Script for running and plotting tau-leap SIR open model with CEA (with
%uncertainty in parameters and costs)

%Clear all the previous variables (important if you change iteration number
%and re-run). Also close previous figures
clear all
close all
%% 

%Specify plotting defaults
set(0,'defaultaxesfontsize',22)
set(0,'defaultlinelinewidth',2)
load('Posterior.mat')
%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Key health economic parameters

%Set up time horizon
NumYears = 20;
maxtime = 365*NumYears + 1;


%Set up number of iterations and timestep and population
iterations = 1000; 
%On my machine total run time is 0.4*NumStrat*iterations
timestep = 1;
N = 1e5;

%Set up discounting and DALY parameters
dw_HI = 0.15;
dw_MI = 0.01;
yll = 0;%Undiscounted average years of life lost per death
r=0.03;%0.035;
yll_disc = sum((1/(1+r)).^([1:yll]-1));
disc = (1/(1+r)).^([1:NumYears]-1);

%Pre-compute costs of strategies
Cost_c = gamrnd(4.5*ones(iterations,1),1/3*ones(iterations,1)); 
Cost_a = gamrnd(45*ones(iterations,1),2/30*ones(iterations,1)); 

Cost_ad = gamrnd(4.5*ones(iterations,1),1/9*ones(iterations,1)); %Exp $0.5
Cost_wt = gamrnd(18*ones(iterations,1),1/9*ones(iterations,1));  %Exp $2

g_c = 0.85;
g_a = 0.2;
h = 0.95;

%Gives start point (index) of each year
t_Yr=[1:365:NumYears*365+1];

%Set up the number of strategies
NumStrats = 6;

%Preallocate matrices for storing DALYs and Costs
CostMat = zeros(iterations,NumStrats);
CostMat_disc = zeros(iterations,NumStrats);
DALYMat = zeros(iterations,NumStrats);
DALYMat_disc = zeros(iterations,NumStrats);

CostMat_avg = zeros(1,NumStrats);
CostMat_disc_avg = zeros(1,NumStrats);
DALYMat_avg = zeros(1,NumStrats);
DALYMat_disc_avg = zeros(1,NumStrats);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the transmission model including running for 5 years to find a
% range of ICs for simulating the strategies
%% Initial Search for Simple Equillibrium

%Given original model parameters
para = struct('beta_c',2/365,'beta_a',1/365,'n_c',0.3,'n_a',0.7,'sigma',1/365,'rho',2/3,'mu',5/365,'k',0.7,'z',0.93,'R0',3);
para0=para;
%Find End Eq
%Compute endemic equilibrium Do Heuristically
M_c_star = 11.5;
M_a_star = 5.6;
l_star = 5.6;

ICs0 = struct('M_c',M_c_star,'M_a',M_a_star,'l',l_star);

%Run ODE model for a lil bit to fine tune equilibrium for given params
[Classes] = ODE_model(para,ICs0,0,5*365);
ICs0 = struct('M_c',Classes.M_c(end),'M_a',Classes.M_a(end),'l',Classes.l(end));
%Define initial conditions as a structure

%% Compute The IC Matrix

%Run the no intervention model for 5 years to generate variation in
%ICs for time 0
for r = 1:iterations 
    %Generate a parameter set using data from the posterior
    para = struct('beta_c',2/365,'beta_a',1/365,'n_c',0.3,'n_a',0.7,'sigma',1/365,'rho',2/3,'mu',5/365,'k',Posterior.k(r),'z',0.93,'R0',Posterior.R_0(r));
    
    %Initial Run to find Endemic Equilibrium. Assume close to that found
    %before.
    [Classes] = ODE_model(para,ICs0,0,365*NumYears);    
    ICs(r) = struct('M_c',Classes.M_c(end),'M_a',Classes.M_a(end),'l',Classes.l(end));
end



%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the ODE model for each startegy for the number of replicates and store the key
% metrics

for Strat = 1:NumStrats
    %Strat
    %Update parameters as relevent
    if Strat==1
        g_c = 0;
        g_a = 0;
        h  = 0;
        
    elseif Strat==2
        g_c = 0.85;
        g_a = 0;
        h = 0.95;
        
        e_c = (1-0.85*0.95);
        e_a = 1;
        
    elseif Strat==3
        g_c = 0.85*0.85;
        g_a = 0;
        h = 0.95;
        
        e_c = (1-0.85*0.95)^2;
        e_a = 1;
     elseif Strat==4
        g_c = 0.85;
        g_a = 0.2;
        h = 0.95;
        
        e_c = (1-0.85*0.95);
        e_a = (1-0.2*0.95);
    elseif Strat==5
        para=para0;
        para.beta_c = (3/5)*para.beta_c;
        para.beta_a = (3/5)*para.beta_a;
    elseif Strat==6
        para=para0;
        para.beta_c = (3/5)*para.beta_c;
    end
    
    
    %Loop through each iteration and re-run model with different IC
    for r = 1:iterations
        [Strat, r] %Helps seeing progress when running.
        %If strat 1, just run the model.
        if Strat == 1
            para = para0;
            para.R0 = Posterior.R_0(r);
            para.k = Posterior.k(r);
            [Classes] = ODE_model(para,ICs(r),0,maxtime);
            M_c_Mat{Strat}(r,:) = Classes.M_c;
            M_a_Mat{Strat}(r,:) = Classes.M_a;
        elseif Strat <=4
            %Set up parameters etc
            para = para0;
            para.R0 = Posterior.R_0(r);
            para.k = Posterior.k(r);
            
            IC = ICs(r);
            M_cVec = [];
            M_aVec = [];
            %Resrat ODE solver every year when Ms are changed
            for Y = 1:NumYears
                %Y
                [Classes1] = ODE_model(para,IC,365*(Y-1),365*Y +1 );
                IC = struct('M_c',e_c*Classes1.M_c(end),'M_a', e_a*Classes1.M_a(end),'l', Classes1.l(end));
                M_cVec = [M_cVec;Classes1.M_c];
                M_aVec = [M_aVec;Classes1.M_a];
            end
            
            M_c_Mat{Strat}(r,:) = M_cVec;
            M_a_Mat{Strat}(r,:) = M_aVec;
        elseif Strat >=4
            %Run ODE solver with altered betas.
            para.R0 = Posterior.R_0(r);
            para.k = Posterior.k(r);
            [Classes] = ODE_model(para,ICs(r),0,maxtime);
            M_c_Mat{Strat}(r,:) = Classes.M_c;
            M_a_Mat{Strat}(r,:) = Classes.M_a;
        end
        
        %Convert Mean Worm Burdens to Infection Dynamics.
        
        %Those in Small Class
        Low_I_c{Strat}(r,:) = my_nbin(0,14,M_c_Mat{Strat}(r,:),para.k,true)*N*para.n_c;
        Low_I_a{Strat}(r,:) = my_nbin(0,14,M_a_Mat{Strat}(r,:),para.k,true)*N*para.n_a;
        %Those in Medium Class
        Mid_I_c{Strat}(r,:) = my_nbin(15,29,M_c_Mat{Strat}(r,:),para.k,true)*N*para.n_c;
        Mid_I_a{Strat}(r,:) = my_nbin(15,29,M_a_Mat{Strat}(r,:),para.k,true)*N*para.n_a;
        %Those in High Class
        High_I_c{Strat}(r,:) = N*para.n_c -  Mid_I_c{Strat}(r,:) - Low_I_c{Strat}(r,:);
        High_I_a{Strat}(r,:) = N*para.n_a -  Mid_I_a{Strat}(r,:) - Low_I_a{Strat}(r,:);
        
        
        

        %Compute treatments, person-years and deaths each year
        for i=1:length(t_Yr)-1
                %Person years infected
                Low_PersonYears_annual(i) = (sum(Low_I_c{Strat}(r,[t_Yr(i): t_Yr(i+1)])) +sum(Low_I_a{Strat}(r,[t_Yr(i): t_Yr(i+1)]))) /365; %divide by 365 to go from person days to person years
                Mid_PersonYears_annual(i) = (sum(Mid_I_c{Strat}(r,[t_Yr(i): t_Yr(i+1)])) +sum(Mid_I_a{Strat}(r,[t_Yr(i): t_Yr(i+1)]))) /365;
                High_PersonYears_annual(i) = (sum(High_I_c{Strat}(r,[t_Yr(i): t_Yr(i+1)])) +sum(High_I_a{Strat}(r,[t_Yr(i): t_Yr(i+1)]))) /365;
                %Treated
                if Strat <=4
                    Treatments_annual_c(i) = (g_c*N*para.n_c);
                    Treatments_annual_a(i) = (g_a*N*para.n_a);
                elseif Strat ==5
                    Treatments_annual_c(i) = ((2/3)*N*para.n_c);
                    Treatments_annual_a(i) = ((2/3)*N*para.n_a);
                elseif Strat ==6
                    Treatments_annual_c(i) = (0.9*N*para.n_c);
                    Treatments_annual_a(i) = (0*N*para.n_a);
                end
        end
    
        %Compute DALYs and store
        DALYs = dw_MI*Mid_PersonYears_annual(i) + dw_MI*High_PersonYears_annual;
        
        DALYMat(r,Strat) = sum(DALYs,2);
        DALYMat_disc(r,Strat) = sum(DALYs.*disc,2);
        
        DALYs_annual{Strat}(r,:)=DALYs;
        
        
        %Compute Costs and store
        Cost=zeros(1,NumYears);
        if Strat == 1
            Cost = 0;
        end
        if Strat==2
            Cost = Treatments_annual_c.*Cost_c(r);
        end
        if Strat==3
            Cost = Treatments_annual_c.*Cost_c(r);
        end
        if Strat==4
            Cost = Treatments_annual_c.*Cost_c(r) + Treatments_annual_a.*Cost_a(r);
        end
        if Strat==5
            Cost = Treatments_annual_c.*Cost_ad(r) + Treatments_annual_a.*Cost_ad(r);
        end
        if Strat==6
            Cost = Treatments_annual_c.*Cost_wt(r);
        end
        
        CostMat(r,Strat) = sum(Cost,2);
        CostMat_disc(r,Strat) = sum(Cost.*disc,2);
        Costs_annual{Strat}(r,:)=Cost;

    end
    CostMat_avg(Strat) = mean(CostMat(:,Strat)); %zeros(1,NumStrats);
    CostMat_disc_avg(Strat) = mean(CostMat_disc(:,Strat)); %zeros(1,NumStrats);
    DALYMat_avg(Strat) = mean(DALYMat(:,Strat));%zeros(1,NumStrats);
    DALYMat_disc_avg(Strat) = mean(DALYMat_disc(:,Strat));%zeros(1,NumStrats);

end

CMap=[228,26,28;55,126,184;77,175,74]./255;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Initial Distribution
low_c = my_nbin(0,14,ICs0.M_c,para0.k,true)*para.n_c*N
mid_c = my_nbin(15,29,ICs0.M_c,para0.k,true)*para.n_c*N
high_c = para.n_c*N - low_c - mid_c

low_a = my_nbin(0,14,ICs0.M_a,para0.k,true)*para.n_a*N
mid_a = my_nbin(15,29,ICs0.M_a,para0.k,true)*para.n_a*N
high_a = para.n_a*N - low_a - mid_a

%bar([low_c mid_c high_c;low_a mid_a high_a])

%Plot distributions
figure(1)
clf
subplot(1,2,1)
hold on
bar([0:50],my_nbin(0,50,ICs0.M_c,para0.k,false));
xline(15,'LineWidth',2,'color','r');
xline(29,'LineWidth',2,'color','r');
hold off
title('Initial Worm Distribution for Children')
xlabel('Number of Worms')
ylabel('Probability')
xlim([0 50])

subplot(1,2,2)
hold on
bar([0:50],my_nbin(0,50,ICs0.M_a,para0.k,false));
xline(15,'LineWidth',2,'color','r');
xline(29,'LineWidth',2,'color','r');
hold off
title('Initial Worm Distribution for Adults')
xlabel('Number of Worms')
ylabel('Probability')
xlim([0 50])

%% Plot Infection Dynamics
figure(2)
clf

subplot(2,1,1)
hold on
plot(mean(Mid_I_c{2}),'-b')
plot(mean(Mid_I_a{2}),'--b')
plot(mean(Mid_I_c{3}),'-r')
plot(mean(Mid_I_a{3}),'--r')
plot(mean(Mid_I_c{4}),'-k')
plot(mean(Mid_I_a{4}),'--k')

%Comment out or in to add extra strategies

% plot(mean(Mid_I_c{5}),'-g')
% plot(mean(Mid_I_a{5}),'--g')
% plot(mean(Mid_I_c{6}),'-y')
% plot(mean(Mid_I_a{6}),'--y')

title('Infection Dynamics for Medium Infected')
xlabel('Time (Days)')
ylabel('Number Infected')
legend('S2 (Children)','S2 (Adults)','S3 (Children)','S3 (Adults)','S4 (Children)','S4 (Adults)')
%legend('S2 (Children)','S2 (Adults)','S3 (Children)','S3 (Adults)','S4 (Children)','S4 (Adults)','S5 (Children)','S5 (Adults)','S6 (Children)','S6 (Adults)')
%xlim([0 inf]);
xlim([0 9000]);


subplot(2,1,2)
hold on
plot(mean(High_I_c{2}),'-b')
plot(mean(High_I_a{2}),'--b')
plot(mean(High_I_c{3}),'-r')
plot(mean(High_I_a{3}),'--r')
plot(mean(High_I_c{4}),'-k')
plot(mean(High_I_a{4}),'--k')

%Comment out or in to add extra strategies

% plot(mean(High_I_c{5}),'-g')
% plot(mean(High_I_a{5}),'--g')
% plot(mean(High_I_c{6}),'-y')
% plot(mean(High_I_a{6}),'--y')
title('Infection Dynamics for Highly Infected')
xlabel('Time (Days)')
ylabel('Number Infected')
legend('S2 (Children)','S2 (Adults)','S3 (Children)','S3 (Adults)','S4 (Children)','S4 (Adults)')
%legend('S2 (Children)','S2 (Adults)','S3 (Children)','S3 (Adults)','S4 (Children)','S4 (Adults)','S5 (Children)','S5 (Adults)','S6 (Children)','S6 (Adults)')
%xlim([0 inf]);
xlim([0 9000]);

%% Plot DALYS and COSTS
figure(3)
subplot(2,1,1)
b=bar(transpose([mean(DALYs_annual{2});mean(DALYs_annual{3});mean(DALYs_annual{4})]))
%b=bar(transpose([mean(DALYs_annual{2});mean(DALYs_annual{3});mean(DALYs_annual{4});mean(DALYs_annual{5});mean(DALYs_annual{6})]))
b(1).FaceColor = ['b'];
b(2).FaceColor = ['r'];
b(3).FaceColor = ['k'];
% b(4).FaceColor = ['g'];
% b(5).FaceColor = ['y'];
title('Yearly DALYs Lost')
xlabel('Time (years)')
ylabel('DALYs Lost')
%legend('S2','S3','S4')
legend('S2','S3','S4','S5','S6')

%Again can change what is commented to add extra strategies


subplot(2,1,2)
b=bar(transpose([mean(Costs_annual{2});mean(Costs_annual{3});mean(Costs_annual{4})]))
%b = bar(transpose([mean(Costs_annual{2});mean(Costs_annual{3});mean(Costs_annual{4});mean(Costs_annual{5});mean(Costs_annual{6})]))
b(1).FaceColor = ['b'];
b(2).FaceColor = ['r'];
b(3).FaceColor = ['k'];
% b(4).FaceColor = ['g'];
% b(5).FaceColor = ['y'];
title('Yearly Costs')
xlabel('Time (years)')
ylabel('Costs ($)')
legend('S2','S3','S4')
%legend('S2','S3','S4','S5','S6')

%Again can change what is commented to add extra strategies

%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute Delta Costs/ Delta DALYs and compute NMBs Q1

%Discounted
%DCostMat_disc = CostMat_disc_avg - repmat(CostMat_disc_avg(:,1),1,NumStrats);
DCost_disc = CostMat_disc_avg - CostMat_disc_avg(1)*ones(1,6);
%DDALYMat_disc = repmat(DALYMat_disc(:,1),1,NumStrats) - DALYMat_disc;
DDALY_disc = DALYMat_disc_avg(1)*ones(1,6) - DALYMat_disc_avg;

%Reorder by Cost
a = DCost_disc(3);
DCost_disc(3) = DCost_disc(2);
DCost_disc(2) = a;

a = DDALY_disc(3);
DDALY_disc(3) = DDALY_disc(2);
DDALY_disc(2) = a;

%ICERs 
ICER(1)=0;
ICER(2)=(DCost_disc(2)-DCost_disc(1))/(DDALY_disc(2)-DDALY_disc(1));
ICER(3)=(DCost_disc(3)-DCost_disc(2))/(DDALY_disc(3)-DDALY_disc(2));
ICER(4)=(DCost_disc(4)-DCost_disc(3))/(DDALY_disc(4)-DDALY_disc(3));

Name = [1 3 2 4];
%Make table

ICERtable = array2table([Name' DCost_disc([1:4])' DDALY_disc([1:4])' ICER']);
ICERtable.Properties.VariableNames={'Strategy','Delta Costs','Delta DALYs','ICER'};
%Is incorrect, Remove dominated stat 3
ICER2([1:2]) = ICER([1:2]);
ICER2(3) = (DCost_disc(4)-DCost_disc(2))/(DDALY_disc(4)-DDALY_disc(2));
DCost_disc2 = DCost_disc([[1:2], 4]);
DDALY_disc2 = DDALY_disc([[1:2], 4]);
Name2 =Name([[1:2], 4]);

%Remake the table
ICERtable2 = array2table([Name2' DCost_disc2' DDALY_disc2' ICER2']);
ICERtable2.Properties.VariableNames={'Strategy','Delta Costs','Delta DALYs','ICER'}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

%Plot CE plane
figure(4)
clf
hold on
%Draw on ICERs
plot(DDALY_disc([1,3,4])./1e3,DCost_disc([1,3,4])./1e6,'-k')

%Plot outcomes as a scatter (with transparency)
for Strat=1:4
   scatter((DALYMat_disc(:,1) - DALYMat_disc(:,Strat))./1e3,-(CostMat_disc(:,1) - CostMat_disc(:,Strat))./1e6)%,[])%,CMap(Strat,:),'filled','markerfacealpha',0.5);
end

%Mark on means
for Strat=1:4
hs(Strat)=scatter(DDALY_disc(Strat)/1e3,DCost_disc(Strat)/1e6);%,[],CMap(Strat,:),'filled');
end

%Added details using other software for ease.
xlabel('DALYs averted (thousands)')
ylabel('Additional costs ($M)')
title('CE Plane')
%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NMB for $1 increments
%N.B. using original order of strategies
Diff_DALY = (DALYMat_disc(:,1) - DALYMat_disc(:,[[1:2],4]));
Diff_Cost = -(CostMat_disc(:,1) - CostMat_disc(:,[[1:2],4]));
w=0:1:1e4; %WTP thresholds
for j=1:length(w)
    for Strat=1:3    
        NMB(:,(Strat-1)*length(w)+j) = w(j)*Diff_DALY(:,Strat) - Diff_Cost(:,Strat);
    end
end

%Check if each strategy is optimal (lowest NBM for each WTP) per replicate
Optimal = zeros(iterations,length(w)*3);
for j=1:length(w)
    for r= 1:iterations
        [im,iy] = max(NMB(r,j:length(w):Strat*length(w)));
        Optimal(r,j+ (iy-1)*length(w)) = 1;
    end
end


%Compute the probability that each strategy is CE
for j=1:length(w)
    ProbCE(:,j) = (sum(Optimal(:,j:length(w):length(w)*3))/iterations)';
end
%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot CEACs
figure(4)
clf
set(gca, 'ColorOrder', CMap, 'NextPlot', 'replacechildren');

%Plot each CEAC 
p=plot(w,ProbCE);
hold on

%Plot CEAF in different segments according to ICERs
% Strategy 1
x=0:round(ICER2(2));
scatter(x,ProbCE(1,x+1),[],CMap(1,:),'filled')
% Strategy 2
x=round(ICER2(2)):round(ICER2(3));
scatter(x,ProbCE(2,x+1),[],CMap(2,:),'filled')
% Strategy 3
x=round(ICER2(3)):max(w);
scatter(x,ProbCE(3,x+1),[],CMap(3,:),'filled')


%Generate an additonal marker in black for the legened
pp = scatter(-1,-1,'filled','k');

axis([ 0 1000 0 1])
ylabel('Probability cost-effective')
xlabel('Willingness to pay ($) per DALY averted')
title('Cost-Effectiveness Acceptability Curves With Frontier')

legend([p; pp],'S1 (Comparator)','S3','S4','CEAF','location','east')

%% 

%Compute Delta Costs/ Delta DALYs and compute NMBs Q2

%Discounted
%DCostMat_disc = CostMat_disc_avg - repmat(CostMat_disc_avg(:,1),1,NumStrats);
DCost_disc = CostMat_disc_avg - CostMat_disc_avg(1)*ones(1,6);
%DDALYMat_disc = repmat(DALYMat_disc(:,1),1,NumStrats) - DALYMat_disc;
DDALY_disc = DALYMat_disc_avg(1)*ones(1,6) - DALYMat_disc_avg;

%Reorder by Cost
C2 = DCost_disc(2);
C3 = DCost_disc(3);
C4 = DCost_disc(4);
C5 = DCost_disc(5);
C6 = DCost_disc(6);

D2 = DDALY_disc(2);
D3 = DDALY_disc(3);
D4 = DDALY_disc(4);
D5 = DDALY_disc(5);
D6 = DDALY_disc(6);

DCost_disc2(2) = C5;
DCost_disc2(3) = C3;
DCost_disc2(4) = C2;
DCost_disc2(5) = C6;
DCost_disc2(6) = C4;

DDALY_disc2(2) = D5;
DDALY_disc2(3) = D3;
DDALY_disc2(4) = D2;
DDALY_disc2(5) = D6;
DDALY_disc2(6) = D4;

%ICERs 
ICER(1)=0;
ICER(2)=(DCost_disc2(2)-DCost_disc2(1))/(DDALY_disc2(2)-DDALY_disc2(1));
ICER(3)=(DCost_disc2(3)-DCost_disc2(2))/(DDALY_disc2(3)-DDALY_disc2(2));
ICER(4)=(DCost_disc2(4)-DCost_disc2(3))/(DDALY_disc2(4)-DDALY_disc2(3));
ICER(5)=(DCost_disc2(5)-DCost_disc2(4))/(DDALY_disc2(5)-DDALY_disc2(4));
ICER(6)=(DCost_disc2(6)-DCost_disc2(5))/(DDALY_disc2(6)-DDALY_disc2(5));
Name = [1 5 3 2 6 4];

ICERtable = array2table([Name' DCost_disc2' DDALY_disc2' ICER']);
ICERtable.Properties.VariableNames={'Strategy','Delta Costs','Delta DALYs','ICER'};

%First table incorrect
%Remove Strat 2 and 6(4th and 5th Element)
DCost_disc2([4:5]) = [];
DDALY_disc2([4:5]) = [];
Name([4:5]) = [];

ICER(1)=0;
ICER(2)=(DCost_disc2(2)-DCost_disc2(1))/(DDALY_disc2(2)-DDALY_disc2(1));
ICER(3)=(DCost_disc2(3)-DCost_disc2(2))/(DDALY_disc2(3)-DDALY_disc2(2));
ICER(4)=(DCost_disc2(4)-DCost_disc2(3))/(DDALY_disc2(4)-DDALY_disc2(3));
ICER([5:6]) = [];
ICERtable = array2table([Name' DCost_disc2' DDALY_disc2' ICER']);
ICERtable.Properties.VariableNames={'Strategy','Delta Costs','Delta DALYs','ICER'};

%Remove 2nd element 5th strat as weakly dominated.
DCost_disc2(2) = [];
DDALY_disc2(2) = [];
Name(2) = [];
ICER(1)=0;
ICER(2)=(DCost_disc2(2)-DCost_disc2(1))/(DDALY_disc2(2)-DDALY_disc2(1));
ICER(3)=(DCost_disc2(3)-DCost_disc2(2))/(DDALY_disc2(3)-DDALY_disc2(2));
ICER(4) = [];
ICERtable = array2table([Name' DCost_disc2' DDALY_disc2' ICER']);
ICERtable.Properties.VariableNames={'Strategy','Delta Costs','Delta DALYs','ICER'}

%% Plot CE Plane
figure(6)
clf
hold on
%Draw on ICERs
plot(DDALY_disc([1,3,4])./1e3,DCost_disc([1,3,4])./1e6,'-k')

%Plot outcomes as a scatter (with transparency)
%for Strat=1:4
%   scatter(DDALYMat_disc(:,Strat)./1e3,DCostMat_disc(:,Strat)./1e6,[],CMap(Strat,:),'filled','markerfacealpha',0.5);
%end

for Strat=1:6
   scatter((DALYMat_disc(:,1) - DALYMat_disc(:,Strat))./1e3,-(CostMat_disc(:,1) - CostMat_disc(:,Strat))./1e6)%,[])%,CMap(Strat,:),'filled','markerfacealpha',0.5);
end

%Mark on means
for Strat=1:6
hs(Strat)=scatter(DDALY_disc(Strat)/1e3,DCost_disc(Strat)/1e6);%,[],CMap(Strat,:),'filled');
end

%legend(hs,'No intervention','Drug 1','Drug 2','location','SouthEast')
xlabel('DALYs averted (thousands)')
ylabel('Additional costs ($M)')
title('CE Plane')

%% 
%NMB for $1 increments
%N.B. using original order of strategies
Diff_DALY = (DALYMat_disc(:,1) - DALYMat_disc(:,[[1:2],4]));
Diff_Cost = -(CostMat_disc(:,1) - CostMat_disc(:,[[1:2],4]));
w=0:1:1e4; %WTP thresholds
for j=1:length(w)
    for Strat=1:3    
        NMB(:,(Strat-1)*length(w)+j) = w(j)*Diff_DALY(:,Strat) - Diff_Cost(:,Strat);
    end
end

%Check if each strategy is optimal (lowest NBM for each WTP) per replicate
Optimal = zeros(iterations,length(w)*3);
for j=1:length(w)
    for r= 1:iterations
        [im,iy] = max(NMB(r,j:length(w):Strat*length(w)));
        Optimal(r,j+ (iy-1)*length(w)) = 1;
    end
end


%Compute the probability that each strategy is CE
for j=1:length(w)
    ProbCE(:,j) = (sum(Optimal(:,j:length(w):length(w)*3))/iterations)';
end
%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot CEACs
figure(7)
clf
set(gca, 'ColorOrder', CMap, 'NextPlot', 'replacechildren');

%Plot each CEAC 
p=plot(w,ProbCE);
hold on

%Plot CEAF in different segments according to ICERs
% Strategy 1
x=0:round(ICER(2));
scatter(x,ProbCE(1,x+1),[],CMap(1,:),'filled')
% Strategy 2
x=round(ICER(2)):round(ICER(3));
scatter(x,ProbCE(2,x+1),[],CMap(2,:),'filled')
% Strategy 3
x=round(ICER(3)):max(w);
scatter(x,ProbCE(3,x+1),[],CMap(3,:),'filled')


%Generate an additonal marker in black for the legened
pp = scatter(-1,-1,'filled','k');

axis([ 0 2e3 0 1])
ylabel('Probability cost-effective')
xlabel('Willingness to pay ($) per DALY averted')
title('Cost-Effectiveness Acceptability Curves With Frontier')

legend([p; pp],'S1 (Comparator)','S3','S4','CEAF','location','east')



