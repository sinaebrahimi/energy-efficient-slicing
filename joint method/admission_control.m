%% Admission Control Heuristic
%s=Tau_tx+Tau_prop+ones(length(K_set),1)*1000*max(Tau); %s_k (1000)
I_max=15; %max allowed iteration of feasibilty problem algorithm
count_run=0;
rejected_slices=zeros(T,K,5);
while true
    count_run=count_run+1 %counter of while loop for convergence of feasibility problem (16)
    %convergence_check_node=sum(eq);
    %% solving s(k)
    disp('elastic variables problem ');
    cvx_begin
    cvx_solver MOSEK
    variable xii_var(T,K,M,N) binary
    variable gamma_var(N) binary
    variable pi_var(N,N,max(max(possible_paths)),T,K,M,M) binary
    variable theta_var(T,K,M,N,M,N) binary %auxiliary variable for xii_var(t,k,m,n)*xii_var(t,k,mm,nn)
    variable elas(N,3) %1:computing, 2:memory, 3:storage
    variable elas_tau(T,K) % Delay
    variable elas_bw(N,N) % BW
    minimize (sum(sum(elas))+sum(sum(elas_tau))+sum(sum(elas_bw)))
    subject to
    
    % Relaxed Constraints
    
    %______________________________________________________
    %         0<=xii_var<=1;
    %         0<=gamma_var<=1;
    %         0<=pi_var<=1;
    %         0<=theta_var<=1;
    %_____________________________________________________________
    
    temp_sum_rate=cvx(zeros(N,N));
    %% Summing up all rates (part of (11))
    for u=1:N
        for uu=1:N
            for t=1:T
                for k=1:K
                    if k<=SliceNum(1,t)
                        for m=1:M
                            if m<=NumReqVMs(t,k)
                                for mm=1:M
                                    if m~=mm
                                        if mm<=NumReqVMs(t,k)
                                            if Vlink_adj(m,mm,t,k)==1
                                                for n=1:N
                                                    for nn=1:N
                                                        for b=best_path_sorted(n,nn,:)
                                                            if b~=0
                                                                Requested_rate_mm=Varpi_vl(m,mm,t,k);
                                                                temp_sum_rate(u,uu)=temp_sum_rate(u,uu)+I_l2p(n,nn,b,u,uu).*pi_var(n,nn,b,t,k,m,mm).*Requested_rate_mm;
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% C1 (2) (Computing, Memory and Storage Capacity of Cloud Nodes)
    for rt=1:TYPE
        for n=1:N
            sum(sum(sum(xii_var(:,:,:,n).*phi_vm(:,:,:,rt)))) <= r_n(rt,n)+elas(n,rt); %const (2) C1
        end
    end
    
    %% C2 and C3(3) (C2: Each VM can only run on one cloud node ,,, C3: The cloud node must be turned-on to place the VM)
    for t=1:T
        for  k=1:K
            for m=1:M
                if m<=NumReqVMs(t,k) && k<=SliceNum(1,t)
                    sum(xii_var(t,k,m,:))==1; %const C2
                    
                    for n=1:N
                        xii_var(t,k,m,n)<=gamma_var(n); %const (3) C3
                    end
                end
            end
        end
    end
    
    
    %% C4 (9) & C5 (10)-->(16) (C4: Each VL is offloaded to one path ,,, C5: Each VL's ingress and egress physical nodes should be the cloud nodes assigned to the corresponding VMs)
    for t=1:T
        for k=1:K
            if k<=SliceNum(1,t)
                for m=1:M
                    if m<=NumReqVMs(t,k)
                        for mm=1:M
                            if m~=mm
                                if mm<=NumReqVMs(t,k)
                                    sum(sum(sum(pi_var(:,:,:,t,k,m,mm))))==1; %C4(9)
                                    for n=1:N
                                        for nn=1:N
                                            %%linearization
                                            for b=1:possible_paths(n,nn)
                                                pi_var(n,nn,b,t,k,m,mm)==theta_var(t,k,m,n,mm,nn); %C5-a (16)
                                                theta_var(t,k,m,n,mm,nn)<=xii_var(t,k,m,n)+1-xii_var(t,k,mm,nn); %C5-b (16)
                                                xii_var(t,k,m,n)<=theta_var(t,k,m,n,mm,nn)+1-xii_var(t,k,mm,nn); %C5-c (16)
                                                theta_var(t,k,m,n,mm,nn)<=xii_var(t,k,mm,nn); %C5-d (16)
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    %% C6 (11) (C6: BW constraint)
    for u=1:N
        for uu=1:N
            temp_sum_rate(u,uu)/2<=BW(u,uu)+elas_bw(u,uu); %C6 (11)
        end
    end
    
    %% C7 (13) (Delay Constraint)
    for n=1:N
        for nn=1:N
            for t=1:T
                for k=1:K
                    if k<=SliceNum(1,t)
                        for m=1:M
                            if m<=NumReqVMs(t,k)
                                for mm=1:M
                                    if mm<=NumReqVMs(t,k)
                                        if m~=mm
                                            if Vlink_adj(m,mm,t,k)==1
                                                for b=best_path_sorted(n,nn,:)
                                                    if b~=0
                                                        propdelay_path(n,nn,b).*pi_var(n,nn,b,t,k,m,mm)<=Tau_max_vl(m,mm,t,k)+elas_tau(t,k); % C7 (13)
                                                        %sum(sum(I_l2p(n,nn,b,:,:).*Tau_prop)).*pi_var(n,nn,b,t,k,m,mm)<=Tau_max_vl(m,mm,t,k)+elas_tau(t,k); % C7 (13)
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    %% Elasticization constraint
    elas >= 0;
    elas_tau>=0;
    elas_bw>=0;
    
    cvx_end
    temp_bw_consumption=sum(sum(temp_sum_rate.*psi));
    disp("AC");
    
    disp(cvx_optval);
    disp(cvx_status);
    
    %% Rejecting the slice bottlenecking the network if there is any
    elas_bw_flag=0;
    elas_tau_flag=0;
    elas_flag=0;
    %if (convergence_check-sum(s))/convergence_check <1e-1|| count_sum>=I_max || convergence_check-sum(s)<1e-2
    if count_run>=I_max || sum(sum(elas_tau))==0 || sum(sum(elas_tau))<1e-3
        elas_tau_flag=1;
    end
    
    if count_run>=I_max || sum(sum(elas_bw))==0 || sum(sum(elas_bw))<1e-3
        elas_bw_flag=1;
    end
    
    if count_run>=I_max || sum(sum(elas))==0 || sum(sum(elas))<1e-3
        elas_flag=1;
    end
    
    
    if elas_tau_flag==1 && elas_flag==1 && elas_bw_flag==1
        break; %There's no need to reject a slice because the problem will be feasible
    else
        %% Find which slice's VL requests is the main bottleneck
        if elas_tau_flag==0
            disp('elas_tau');
            disp(elas_tau);
            
            [tenant_arr, slicenumber_arr]=find(elas_tau==max(max(elas_tau)));
            tenant=tenant_arr(1);
            slicenumber=slicenumber_arr(1);
            %% Rejecting the slice which its VL requests is bottleneck
            
            for rt=1:TYPE
                for m=1:NumReqVMs(tenant,slicenumber)
                    phi_vm(tenant,slicenumber,m,rt)=0;
                end
            end
            NumReqVMs(tenant,slicenumber)=0;
            Vlink_adj(:,:,tenant,slicenumber)=0;
            Varpi_vl(:,:,tenant,slicenumber)=0;
            Tau_max_vl(:,:,tenant,slicenumber)=0;
            rejected_slices(tenant,slicenumber,4)=1;
            
            disp('Rejecting slice due to Tau_VL:');
            disp('Tenant #');
            disp(tenant);
            disp('SliceNumber #');
            disp(slicenumber);
            %break;
            %%
        elseif elas_bw_flag==0
            disp('elas_bw');
            disp(elas_bw);
            %% Rejecting the slice which has more requested data rate for its VLs
            %             [FirstNode,SecondNode]=find(elas_bw==max(max(elas_bw)));
            sum_used_data_rate=zeros(T,K);
            %if FirstNode~=SecondNode
            %                 if temp_sum_rate(FirstNode,SecondNode)/2>BW(FirstNode,SecondNode)
            %if temp_sum_rate(FirstNode,SecondNode)>BW(FirstNode,SecondNode)
            for t=1:T
                for k=1:K
                    if k<=SliceNum(1,t)
                        sum_used_data_rate(t,k)=sum(sum(Varpi_vl(:,:,t,k)));
                    end
                end
            end
            [tenant_arr,slicenumber_arr]=find(sum_used_data_rate==max(max(sum_used_data_rate)));
            
            tenant=tenant_arr(1);
            slicenumber=slicenumber_arr(1);
            %
            %                 else
            %                     break;
            %                 end
            %else
            
            %end
            
            %% Rejecting the slice which has  more requested data rate for its VLs
            
            for rt=1:TYPE
                for m=1:NumReqVMs(tenant,slicenumber)
                    phi_vm(tenant,slicenumber,m,rt)=0;
                end
            end
            NumReqVMs(tenant,slicenumber)=0;
            Vlink_adj(:,:,tenant,slicenumber)=0;
            Tau_max_vl(:,:,tenant,slicenumber)=0;
            Varpi_vl(:,:,tenant,slicenumber)=0;
            rejected_slices(tenant,slicenumber,5)=1;
            
            disp('Rejecting slice due to BW bottlenecks:');
            disp('Tenant #');
            disp(tenant);
            disp('SliceNumber #');
            disp(slicenumber);
            % break;
            %%
        elseif elas_flag==0
            disp('elas');
            disp(elas);
            which_constraint=find(sum(elas)==max(sum(elas)));
            %% Find the slice which has more requested computing capacity
            if which_constraint==1 %Computing
                disp('elas Computing');
                disp(elas);
                % x=[];xx=[];
                for t=1:T
                    for k=1:K
                        x(t,k)=sum(phi_vm(t,k,:,1));
                    end
                end
                [tenant_arr, slicenumber_arr]=find(x==max(max(x)));
                tenant=tenant_arr(1);
                slicenumber=slicenumber_arr(1);
                % Rejecting the slice which has more requested computing capacity
                
                for rt=1:TYPE
                    for m=1:NumReqVMs(tenant,slicenumber)
                        phi_vm(tenant,slicenumber,m,rt)=0;
                    end
                end
                
                NumReqVMs(tenant,slicenumber)=0;
                Vlink_adj(:,:,tenant,slicenumber)=0;
                Tau_max_vl(:,:,tenant,slicenumber)=0;
                Varpi_vl(:,:,tenant,slicenumber)=0;
                rejected_slices(tenant,slicenumber,1)=1;
                
                disp('Rejecting slice due to Computing capacity constraints:');
                disp('Tenant #');
                disp(tenant);
                disp('SliceNumber #');
                disp(slicenumber);
                
                %%  Find the slice which has more requested memory capacity
            elseif which_constraint==2 %Memory
                disp('elas Memory');
                disp(elas);
                for t=1:T
                    for k=1:K
                        xx(t,k)=sum(phi_vm(t,k,:,2));
                    end
                end
                [tenant_arr, slicenumber_arr]=find(xx==max(max(xx)));
                tenant=tenant_arr(1);
                slicenumber=slicenumber_arr(1);
                
                %% Rejecting the slice which has more requested memory capacity
                for rt=1:TYPE
                    for m=1:NumReqVMs(tenant,slicenumber)
                        phi_vm(tenant,slicenumber,m,rt)=0;
                    end
                end
                
                NumReqVMs(tenant,slicenumber)=0;
                Vlink_adj(:,:,tenant,slicenumber)=0;
                Tau_max_vl(:,:,tenant,slicenumber)=0;
                Varpi_vl(:,:,tenant,slicenumber)=0;
                rejected_slices(tenant,slicenumber,2)=1;
                
                disp('Rejecting slice due to Memory capacity constraints:');
                disp('Tenant #');
                disp(tenant);
                disp('SliceNumber #');
                disp(slicenumber);
                
                %% Rejecting the slice which has more requested storage capacity
            elseif which_constraint==3 %Storage
                disp('elas Storage');
                disp(elas);
                for t=1:T
                    for k=1:K
                        xxx(t,k)=sum(phi_vm(t,k,:,3));
                    end
                end
                [tenant_arr, slicenumber_arr]=find(xxx==max(max(xxx)));
                tenant=tenant_arr(1);
                slicenumber=slicenumber_arr(1);
                
                %% Rejecting the slice which has more requested memory capacity
                for rt=1:TYPE
                    for m=1:NumReqVMs(tenant,slicenumber)
                        phi_vm(tenant,slicenumber,m,rt)=0;
                    end
                end
                
                NumReqVMs(tenant,slicenumber)=0;
                Vlink_adj(:,:,tenant,slicenumber)=0;
                Tau_max_vl(:,:,tenant,slicenumber)=0;
                Varpi_vl(:,:,tenant,slicenumber)=0;
                rejected_slices(tenant,slicenumber,3)=1;
                %break;
                disp('Rejecting slice due to Computing capacity constraints:');
                disp('Tenant #');
                disp(tenant);
                disp('SliceNumber #');
                disp(slicenumber);
                
            end %ending if of com.mem,sto (which_constraint)
        end %ending if of tau,bw,or node capacity
    end %ending if for breaks
end %ending while


%% Calculating the acceptance ratio

if sum(sum(sum(rejected_slices)))==0
    Acceptance_ratio_Overall(T)=1    
    
    Number_of_Accepted_Slices(T)=sum(SliceNum)
    Number_of_Rejected_Slices(T)=sum(sum(sum(rejected_slices)))
else
    Acceptance_ratio_Overall(T)=1-(sum(sum(sum(rejected_slices)))/sum(SliceNum))
    
    Number_of_Accepted_Slices(T)=sum(SliceNum)-sum(sum(sum(rejected_slices)))
    Number_of_Rejected_Slices(T)=sum(sum(sum(rejected_slices)))
    %Acceptance_ratio(mc,T)=1-(sum(sum(rejected_slices))/sum(SliceNum))
end

% 
if sum(sum(rejected_slices(:,:,1)))==0
    Acceptance_ratio_Node_Com(T)=1
else
    Acceptance_ratio_Node_Com(T)=1-(sum(sum(rejected_slices(:,:,1)))/sum(SliceNum))
end

% 
if sum(sum(rejected_slices(:,:,2)))==0
    Acceptance_ratio_Node_Mem(T)=1
else
    Acceptance_ratio_Node_Mem(T)=1-(sum(sum(rejected_slices(:,:,2)))/sum(SliceNum))
end

% 
if sum(sum(rejected_slices(:,:,3)))==0
    Acceptance_ratio_Node_Sto(T)=1
else
    Acceptance_ratio_Node_Sto(T)=1-(sum(sum(rejected_slices(:,:,3)))/sum(SliceNum))
end


if sum(sum(rejected_slices(:,:,1)))==0 && sum(sum(rejected_slices(:,:,2)))==0 && sum(sum(rejected_slices(:,:,3)))==0
    Acceptance_ratio_Node_Com(T)=1
else
    Acceptance_ratio_Node_Com(T)=1-((sum(sum(rejected_slices(:,:,1)))+sum(sum(rejected_slices(:,:,2)))+sum(sum(rejected_slices(:,:,3))))/sum(SliceNum))
end


%% Link Acceptance Ratio

if sum(sum(rejected_slices(:,:,4)))==0
    Acceptance_ratio_Link_Tau(T)=1
else
    Acceptance_ratio_Link_Tau(T)=1-(sum(sum(rejected_slices(:,:,4)))/sum(SliceNum))
end


if sum(sum(rejected_slices(:,:,5)))==0
    Acceptance_ratio_Link_BW(T)=1
else
    Acceptance_ratio_Link_BW(T)=1-(sum(sum(rejected_slices(:,:,5)))/sum(SliceNum))
end

if sum(sum(rejected_slices(:,:,4)))==0 && sum(sum(rejected_slices(:,:,5)))==0 
    Acceptance_ratio_Link(T)=1
else
    Acceptance_ratio_Link(T)=1-((sum(sum(rejected_slices(:,:,4)))+sum(sum(rejected_slices(:,:,5))))/sum(SliceNum))
end


% if sum(sum(rejected_slices))==0
%     Acceptance_ratio(mc,T)=1
% else
%     Acceptance_ratio(mc,T)=1-(sum(sum(rejected_slices))/sum(SliceNum))
%     %Acceptance_ratio(mc,T)=1-(sum(sum(rejected_slices))/sum(SliceNum))
% end
