%% Solving with MOSEK (xii subproblem) Node selection

cvx_begin
cvx_solver MOSEK
variable xii_var(T,K,M,N) binary
variable gamma_var(N) binary

node_comp_capacity_used=cvx(zeros(1,N));
for n=1:N
    node_comp_capacity_used(n)=sum(sum(sum(xii_var(:,:,:,n).*phi_vm(:,:,:,1))));
end

%minimize sum(temp) %This objective will turn on all cloud nodes!! thus
%this is not energy-efficient enough
minimize sum(((P_max-P_idle)./(r_n(1,:))).*node_comp_capacity_used+((gamma_var)'.*P_idle)) %This objective will try to minimze the number of turned-on cloud nodes

subject to
for rt=1:TYPE
    for n=1:N
        sum(sum(sum(xii_var(:,:,:,n).*phi_vm(:,:,:,rt)))) <= r_n(rt,n); %const (2)
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

cvx_end


%Normalizing xii_var because CVX sometimes approximates the binary value
%with values close to zero or one (e.g. 2.62e-16)
xii_var(xii_var<0.9)=0;
xii_var(xii_var>0.9)=1;

gamma_var(gamma_var<0.9)=0;
gamma_var(gamma_var>0.9)=1;