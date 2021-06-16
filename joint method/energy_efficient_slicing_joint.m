clc
clear
close all
cvx_quiet true

%% Initialization
%N=6;%input('number of cloud nodes: '); %no. of cloud nodes in all of our infrastructure
% TYPE=3; %1 is computing, 2 is memory, 3 is storage
% T=6; %input('number of tenants: '); %number of tenants
% K=5;%input('max number of slices for a tenant: '); %max no. of slices for each tenant
% SliceNum=randi([1 K],1,T);%Number of slices for each tenant
% M=7;%input('max nu mber of VMs for each slice: '); %MaxVMs  %max no. of slices for a tenant
% NumReqVMs=randi([1 M],T,K);%Number of requested VMs for each slice (for determined t and k)  randi(MaxVMs,[3,2,3])

% Physical links initialization
% ini_physical_links;
% 
% % Physical Nodes resources
% ini_physical_nodes;
% 
% % %SAVE PHYSICAL INITIALIZATION%
% filename = 'ini_phy_lesser_resources.mat';
% save(filename);
load('ini_phy_lesser_resources.mat');
%% Begin Looping

%MC=1; %Monte carlo runs
TT=16;%Max tenants

overall_cost=zeros(1,TT);
links_cost=zeros(1,TT);
nodes_cost=zeros(1,TT);

used_nodes=zeros(1,TT);
sum_comp_used=zeros(1,TT);
sum_power_used=zeros(1,TT);
sum_rate_used=zeros(1,TT);

sum_bw_consumption=zeros(1,TT);
avg_time=zeros(1,TT);

%%%%
%Acceptance_ratio=zeros(TT);


Acceptance_ratio_Link=zeros(1,TT);
Acceptance_ratio_Link_Tau=zeros(1,TT);
Acceptance_ratio_Link_BW=zeros(1,TT);

Acceptance_ratio_Node=zeros(1,TT);
Acceptance_ratio_Node_Com=zeros(1,TT);
Acceptance_ratio_Node_Mem=zeros(1,TT);
Acceptance_ratio_Node_Sto=zeros(1,TT);


Acceptance_ratio_Overall=zeros(1,TT);
Number_of_Accepted_Slices=zeros(1,TT);
Number_of_Rejected_Slices=zeros(1,TT);


diary;
diary_lesser_resources_27Mordad = 'diary_lesser_resources_27Mordad.txt';
diary diary_lesser_resources_27Mordad;
% load('first_figs1to8_lessresources.mat');
% TT=24;
%for mc=1:MC
    
    for T=1:TT
        tic
%         disp('Run#:');
%         disp(mc);
        
        
        disp('Tenant no. (Slice number=tenants*1):');
        disp(T);
        
        % Slice requests
        %Initialization of VMs
        ini_virtual_machines;
        %Initialization of VLs
        ini_virtual_links;
        %End of initialization of requests
        
        % Solution
        % Admission Control via Elasticization
        admission_control;
        
        % Joint Problem
%         weight_node=1;
%         weight_link=10^-3;%%%%%%%%%%%%%%%%%%%%%%%
%         

        weight_node=1;
        weight_link=9*10^-5;%%%%%%%%%%%%%%%%%%%%%%%

        joint_problem;
        used_nodes(T)=sum(gamma_var)
        sum_comp_used(T)=sum(node_comp_capacity_used)
        
        sum_rate_used(T)=sum(sum(temp_sum_rate))
        sum_bw_consumption(T)=sum(sum(temp_sum_rate.*psi))
        sum_power_used(T)=sum(((P_max-P_idle)./(r_n(1,:))).*node_comp_capacity_used+((gamma_var)'.*P_idle))
        links_cost(T)=weight_link.*sum_bw_consumption(T)
        
        disp("joint problem cost");
        overall_cost(T)=cvx_optval
        nodes_cost(T)=overall_cost(T)-links_cost(T)
        
        toc
        avg_time(T)=toc        
    end
    %toc
    %avg_time(mc)=toc/TT
    
%end


set(gca,'FontSize',20);

figure (1)
% plot(1:14,overall_cost(1:14));
plot(1:TT,overall_cost);%plot(3*(1:TT),smooth(mean(cost_pi)));
title('InP''s Cost')
ylabel('C_{Total}')
xlabel('Number of Slices') %(\sum_{t\in\mathcal{{T}}}\sum_{k\in\mathcal{K_\text{t}}})')
hold on;

figure (2)
plot(1:TT,ceil(used_nodes));
% plot(1:14,ceil(used_nodes(1:14)));
%plot(K*(1:TT),mean(used_nodes));
set(gca, 'YTick', 0:N) %not showing numbers like 2.2
title('Used Cloud Nodes')
ylabel('sum of \gamma')
xlabel('Number of Slices')
hold on;

figure (3)
plot(1:TT,sum_comp_used);
% plot(1:14,sum_comp_used(1:14));
title('Computing Capacity Used')
ylabel('sum of all \phi^{Com} (KHz)')
xlabel('Number of Slices')
hold on;


figure (4)
% plot(1:14,sum_rate_used(1:14));
plot(1:TT,sum_rate_used);
title('Total Used Bandwidth')
ylabel('Sum of all \varpi (Kbps)')
xlabel('Number of Slices')
hold on;


figure (5)
% plot(1:14,sum_bw_consumption(1:14));
plot(1:TT,sum_bw_consumption);
title('Total Bandwidth Consumption Cost')
ylabel('Total BW Consumption Cost ( \psi * \sum \varpi )')
ylabel('Total BW Consumption Cost')
xlabel('Number of Slices')
hold on;

figure (6)
% plot(1:14,avg_time(1:14));
plot(1:TT,avg_time);
title('Average Execution Time')
ylabel('Avg Execution Time (s)')
xlabel('Number of Slices')
hold on;

figure (7)
plot(1:TT,Acceptance_ratio_Overall);
%plot(1:14,Acceptance_ratio(1:14));
title('Average Acceptance Ratio of the Slice Requests')
ylabel('Acceptance Ratio')
xlabel('Number of Slices')
hold on;


figure (8)
plot(1:TT,sum_power_used);
%plot(1:14,sum_power_used(1:14));
title('Power Consumed by All Cloud Nodes')
ylabel('sum of all P_{n} (W)')
xlabel('Number of Slices')
hold on;


figure (9)
% plot(1:16,overall_cost(1:16));
plot(1:TT,links_cost);%plot(3*(1:TT),smooth(mean(cost_pi)));
title('Cost of Links')
ylabel('C_{Links}')
xlabel('Number of Slices') %(\sum_{t\in\mathcal{{T}}}\sum_{k\in\mathcal{K_\text{t}}})')
hold on;

figure (10)
% plot(1:16,overall_cost(1:16));
plot(1:TT,nodes_cost);%plot(3*(1:TT),smooth(mean(cost_pi)));
title('Cost of Nodes')
ylabel('C_{Nodes}')
xlabel('Number of Slices') %(\sum_{t\in\mathcal{{T}}}\sum_{k\in\mathcal{K_\text{t}}})')
hold on;

figure (11)
plot(1:TT,Acceptance_ratio_Node);
%plot(1:16,Acceptance_ratio(1:16));
title('Average Acceptance Ratio of the Slice Requests (Nodes Subproblem)')
ylabel('Acceptance Ratio')
xlabel('Number of Slices')
hold on;

figure (12)
plot(1:TT,Acceptance_ratio_Node_Com);
%plot(1:16,Acceptance_ratio(1:16));
title('Average Acceptance Ratio of the Slice Requests (Nodes Subproblem (Computing Resources))')
ylabel('Acceptance Ratio')
xlabel('Number of Slices')
hold on;

figure (13)
plot(1:TT,Acceptance_ratio_Node_Mem);
%plot(1:16,Acceptance_ratio(1:16));
title('Average Acceptance Ratio of the Slice Requests (Nodes Subproblem (Memory Resources))')
ylabel('Acceptance Ratio')
xlabel('Number of Slices')
hold on;

figure (14)
plot(1:TT,Acceptance_ratio_Node_Sto);
%plot(1:16,Acceptance_ratio(1:16));
title('Average Acceptance Ratio of the Slice Requests (Nodes Subproblem (Storage Resources))')
ylabel('Acceptance Ratio')
xlabel('Number of Slices')
hold on;

figure (15)
plot(1:TT,Acceptance_ratio_Link);
%plot(1:16,Acceptance_ratio(1:16));
title('Average Acceptance Ratio of the Slice Requests (Links Subproblem)')
ylabel('Acceptance Ratio')
xlabel('Number of Slices')
hold on;

figure (16)
plot(1:TT,Acceptance_ratio_Link_Tau);
%plot(1:16,Acceptance_ratio(1:16));
title('Average Acceptance Ratio of the Slice Requests (Links Subproblem (Tau))')
ylabel('Acceptance Ratio')
xlabel('Number of Slices')
hold on;

figure (17)
plot(1:TT,Acceptance_ratio_Link_BW);
%plot(1:16,Acceptance_ratio(1:16));
title('Average Acceptance Ratio of the Slice Requests (Links Subproblem (BW))')
ylabel('Acceptance Ratio')
xlabel('Number of Slices')
hold on;

figure (18)
% plot(1:16,overall_cost(1:16));
plot(1:TT,Number_of_Accepted_Slices);%plot(3*(1:TT),smooth(mean(cost_pi)));
title('Number of Accepted Slices')
ylabel('Accepted Slices')
xlabel('Number of Slices') %(\sum_{t\in\mathcal{{T}}}\sum_{k\in\mathcal{K_\text{t}}})')
hold on;

figure (19)
plot(1:TT,Number_of_Rejected_Slices);%plot(3*(1:TT),smooth(mean(cost_pi)));
title('Number of Rejected Slices')
ylabel('Rejected Slices')
xlabel('Number of Slices') %(\sum_{t\in\mathcal{{T}}}\sum_{k\in\mathcal{K_\text{t}}})')
hold on;


save joint_fixed_resources;
diary off;