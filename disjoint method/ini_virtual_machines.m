%% Initializing slice requests
%T=3; %input('number of tenants: '); %number of tenants
K=1;%input('max number of slices for a tenant: '); %max no. of slices for each tenant
SliceNum=randi([K K],1,T);%Number of slices for each tenant
M=3;%input('max number of VMs for each slice: '); %MaxVMs  %max no. of slices for a tenant
NumReqVMs=randi([M M],T,K);%Number of requested VMs for each slice (for determined t and k)  randi(MaxVMs,[3,2,3])

%% Virtual Nodes Initialization
for t=1:T
    for k=1:K
        if k>SliceNum(1,t)
            NumReqVMs(t,k)=0;%Making NumReqVMs zero wherever a slice isn't defined for a tenant t
        end
    end
end

%% Requested VM Capacities
phi_vm=zeros([T,K,M,TYPE]);%%Requested VM capacities (1:Com,2:Mem,3:Sto) for each VM of slice t,k

VComSet=[1000]; %A feasible and sense-making set for Servers' computing capacity (in MHz)

VMemSet=[64]; %A feasible and sense-making set for Servers' memory capacity (in GB)

VStoSet=[120]; %A feasible and sense-making set for Servers' storage capacity (in GB)
% 
% VComSet=[400,600,800,1000,1200,1400,1600,1800,2000,2200,...
%     2400,2600,2800,3000,3200,3400,3600,3800]; %A feasible and sense-making set for Servers' computing capacity (in MHz)
% 
% VMemSet=[1,2,4,6,8,10,12,14,16,20,24,28,32,40,48,64,...
%     80,96,112,128,160,192,224,256]; %A feasible and sense-making set for Servers' memory capacity (in GB)
% 
% VStoSet=[4,8,16,32,48,64,96,120,128,250,256,320,500,750,1000,1500,2000]; %A feasible and sense-making set for Servers' storage capacity (in GB)


for t=1:T
    for k=1:K
        for m=1:M
            if m<=NumReqVMs(t,k) && k<=SliceNum(1,t) %Wherever the numbers of VMs are exceeded, we zeroed the requested capacity
                SelectedVC_Com=randi(length(VComSet));
                phi_vm(t,k,m,1)=VComSet(SelectedVC_Com);  %requested Computing Capacity for the requested VM (in MHz)
                
                SelectedVC_Mem=randi(length(VMemSet));
                phi_vm(t,k,m,2)=VMemSet(SelectedVC_Mem); %requested Memory Capacity for the requested VM (in GB)
                
                SelectedVC_Sto=randi(length(VStoSet));
                phi_vm(t,k,m,3)=VStoSet(SelectedVC_Sto); %requested Storage Capacity for the requested VM (in GB)
            end
        end
    end
end