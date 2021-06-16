%% Initializing Cloud Nodes
%% Initializing Cloud Nodes Capacities (In terms of Computing, Memory and Storage)
TYPE=3; %1 is computing, 2 is memory, 3 is storage
r_n=zeros(TYPE,N); %Comp,Mem, Storage resources of node n
% 
% PhyComSet=[4000,4200,4400,4600,4800,5000,5200,5400,5600,5800,...
%     6000,6200,6400,6600,6800,7000,7200,7400,7600,7800,...
%     8000,8200,8400,8600,8800,9000,9200,9400,9600,9800]; %A feasible and sense-making set for Servers' computing capacity (in MHz)
% 
% PhyMemSet=[256,320,384,452,512,640,768,896,1024,1280,1536,...
%     1792,2048,2560,3072,3584,4096]; %A feasible and sense-making set for Servers' memory capacity (in GB)
% 
% PhyStoSet=[2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,...
%     8000,9000,10000]; %A feasible and sense-making set for Servers' storage capacity (in GB)

PhyComSet=[7000]; %A feasible and sense-making set for Servers' computing capacity (in MHz)

PhyMemSet=[800]; %A feasible and sense-making set for Servers' memory capacity (in GB)

PhyStoSet=[2000]; %A feasible and sense-making set for Servers' storage capacity (in GB)
%% Assigning Capacities randomly
for n=1:N
    SelectedC_Com=randi(length(PhyComSet));
    r_n(1,n)=PhyComSet(SelectedC_Com);%SelectedC_Com(pos); %Computing capacity of node n (in MHz)
    
    SelectedC_Mem=randi(length(PhyMemSet));
    r_n(2,n)=PhyMemSet(SelectedC_Mem); %Memory capacity of node n (in GB)
    
    SelectedC_Sto=randi(length(PhyStoSet));
    r_n(3,n)=PhyStoSet(SelectedC_Sto); %Storage capacity of node n (in GB)
end

%% Power of each cloud node
P_max=zeros(1,N);
P_max=r_n(1,:).*0.1; %Power of a node is related to its CPU capacity (in Watts)
P_idle=zeros(1,N);
P_idle=P_max./4;