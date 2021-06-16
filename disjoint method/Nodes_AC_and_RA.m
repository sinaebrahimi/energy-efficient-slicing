%% Admission Control Heuristic
%s=Tau_tx+Tau_prop+ones(length(K_set),1)*1000*max(Tau); %s_k (1000)
I_max=15; %max allowed iteration of feasibilty problem algorithm
count_run_Node_AC=0;
rejected_slices_node=zeros(T,K,TYPE);

while true
    count_run_Node_AC=count_run_Node_AC+1 %counter of while loop for convergence of feasibility problem (16)
    %convergence_check_node=sum(eq);
    %% solving s(k)
    disp('elastic variables problem for NODES');
    cvx_begin
    cvx_solver MOSEK
    variable xii_var(T,K,M,N) binary
    variable gamma_var(N) binary
    %variable pi_var(N,N,max(max(possible_paths)),T,K,M,M) binary
    %variable theta_var(T,K,M,N,M,N) binary %auxiliary variable for xii_var(t,k,m,n)*xii_var(t,k,mm,nn)
    variable elas(N,3) %1:computing, 2:memory, 3:storage
    minimize sum(sum(elas))
    subject to
    
    for rt=1:TYPE
        for n=1:N
            sum(sum(sum(xii_var(:,:,:,n).*phi_vm(:,:,:,rt)))) <= r_n(rt,n)+elas(n,rt); %const (2)
        end
    end
    
    for t=1:T
        for  k=1:K
            for m=1:M
                if m<=NumReqVMs(t,k) && k<=SliceNum(1,t)
                    sum(xii_var(t,k,m,:))==1; %const (3)
                end
                for n=1:N
                    if m<=NumReqVMs(t,k) && k<=SliceNum(1,t)
                        xii_var(t,k,m,n)<=gamma_var(n); %const (12)
                    end
                end
            end
        end
    end

    
    %% Elasticization constraint
    elas >= 0;
    
    cvx_end
    
    
    disp("AC Node");
    
    disp(cvx_optval);
    disp(cvx_status);
    
    %% Rejecting the slice bottlenecking the network if there is any

    elas_flag=0;

    if count_run_Node_AC>=I_max || sum(sum(elas))==0 || sum(sum(elas))<1e-3
        elas_flag=1;
    end
    
    
    if elas_flag==1
        break; %There's no need to reject a slice because the problem will be feasible
    else
        %% Find which slice's VL requests is the main bottleneck
        if elas_flag==0
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
                rejected_slices_node(tenant,slicenumber,1)=1;
                
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
                rejected_slices_node(tenant,slicenumber,2)=1;
                
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
                rejected_slices_node(tenant,slicenumber,3)=1;
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
if sum(sum(sum(rejected_slices_node)))==0
    Acceptance_ratio_Node(T)=1
else
    Acceptance_ratio_Node(T)=1-(sum(sum(sum(rejected_slices_node)))/sum(SliceNum))
    %Acceptance_ratio(mc,T)=1-(sum(sum(rejected_slices))/sum(SliceNum))
end

% 
if sum(sum(rejected_slices_node(:,:,1)))==0
    Acceptance_ratio_Node_Com(T)=1
else
    Acceptance_ratio_Node_Com(T)=1-(sum(sum(rejected_slices_node(:,:,1)))/sum(SliceNum))
end

% 
if sum(sum(rejected_slices_node(:,:,2)))==0
    Acceptance_ratio_Node_Mem(T)=1
else
    Acceptance_ratio_Node_Mem(T)=1-(sum(sum(rejected_slices_node(:,:,2)))/sum(SliceNum))
end

% 
if sum(sum(rejected_slices_node(:,:,3)))==0
    Acceptance_ratio_Node_Sto(T)=1
else
    Acceptance_ratio_Node_Sto(T)=1-(sum(sum(rejected_slices_node(:,:,3)))/sum(SliceNum))
end

xii_subproblem;
