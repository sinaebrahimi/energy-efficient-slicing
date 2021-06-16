%% MOSEK Joint Problem
cvx_begin
cvx_solver MOSEK
%variable xii_var(T,K,M,N) binary
%variable gamma_var(N) binary
variable pi_var(N,N,max(max(possible_paths)),T,K,M,M) binary
%variable theta_var(T,K,M,N,M,N) binary %auxiliary variable for xii_var(t,k,m,n)*xii_var(t,k,mm,nn)

% node_comp_capacity_used=cvx(zeros(1,N));
% for n=1:N
%     node_comp_capacity_used(n)=sum(sum(sum(xii_var(:,:,:,n).*phi_vm(:,:,:,1))));
% end

%BW_consumption_cost=0;
temp_sum_rate=cvx(zeros(N,N));
%% Summing up all rates (part of (11))
for u=1:N
    for uu=1:N
        for t=1:T
            for  k=1:K
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

% %% BW Consumption Cost
% for u=1:N
%     for uu=1:N
%         BW_consumption_cost=BW_consumption_cost+temp_sum_rate(u,uu).*psi(u,uu); %(12)
%     end
% end

%% Objective Function (14) and (15) (A Mixture of the cost of BW consumption & 
%minimize weight_node.*(sum(((P_max-P_idle)./(r_n(1,:))).*node_comp_capacity_used+((gamma_var)'.*P_idle)))+weight_link.*BW_consumption_cost
minimize weight_node.*nodes_problem_cost(T)+weight_link.*sum(sum(temp_sum_rate.*psi))

subject to
% %% C1 (2) (Computing, Memory and Storage Capacity of Cloud Nodes)
% for rt=1:TYPE
%     for n=1:N
%         sum(sum(sum(xii_var(:,:,:,n).*phi_vm(:,:,:,rt)))) <= r_n(rt,n); %const (2) C1
%     end
% end

%% C2 and C3(3) (C2: Each VM can only run on one cloud node ,,, C3: The cloud node must be turned-on to place the VM)
% for t=1:T
%     for  k=1:K
%         for m=1:M
%             if m<=NumReqVMs(t,k) && k<=SliceNum(1,t)
%                 sum(xii_var(t,k,m,:))==1; %const C2
%                 
%                 for n=1:N
%                     xii_var(t,k,m,n)<=gamma_var(n); %const (3) C3
%                 end
%             end
%         end
%     end
% end


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
                                        for b=1:possible_paths(n,nn)
                                            %NO need to relax because
                                                %in this subproblem, xii is
                                                %not a variable!
                                                pi_var(n,nn,b,t,k,m,mm)== xii_var(t,k,m,n).*xii_var(t,k,mm,nn);%C5
                                                
%                                             pi_var(n,nn,b,t,k,m,mm)==theta_var(t,k,m,n,mm,nn); %C5-a (16)
%                                             theta_var(t,k,m,n,mm,nn)<=xii_var(t,k,m,n)+1-xii_var(t,k,mm,nn); %C5-b (16)
%                                             xii_var(t,k,m,n)<=theta_var(t,k,m,n,mm,nn)+1-xii_var(t,k,mm,nn); %C5-c (16)
%                                             theta_var(t,k,m,n,mm,nn)<=xii_var(t,k,mm,nn); %C5-d (16)
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
        temp_sum_rate(u,uu)/2<=BW(u,uu); %C6 (11)
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
                                                    propdelay_path(n,nn,b).*pi_var(n,nn,b,t,k,m,mm)<=Tau_max_vl(m,mm,t,k); % C7 (13)

                                                   % sum(sum(I_l2p(n,nn,b,:,:).*Tau_prop)).*pi_var(n,nn,b,t,k,m,mm)<=Tau_max_vl(m,mm,t,k); % C7 (13)
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

cvx_end

disp('temp sum rate in pi_subproblem');
disp(temp_sum_rate);