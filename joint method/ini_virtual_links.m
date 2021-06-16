%%  Virtual links initialization
Vlink_adj=zeros(M,M,T,K);%adj matrix for S_t,k's VM connectivity
Varpi_vl=zeros(M,M,T,K);%matrix for S_t,k's VM m to mm's requested data rate
Tau_max_vl=1000*ones(M,M,T,K);%matrix for S_t,k's VM mto mm's requested data rate

for t=1:T
    for k=1:K
        while true
            flag_connectedvgraph=zeros(NumReqVMs(t,k)-1,1);
            temp_vadj=randi([0 1],NumReqVMs(t,k)); %adj matrix of nfv_enabled nodes of Graph G
            temp_vadj=triu(temp_vadj,1)+triu(temp_vadj)';%making A symmetric (it is necessary for an adj matrix for an undirected graph)
            temp_vadj=temp_vadj-diag(diag(temp_vadj)); %diag entries=0
            %parfor y=2:NumReqVMs(t,k)
            for y=2:NumReqVMs(t,k)
                paths2y=pathbetweennodes(temp_vadj,1,y);
                flag_connectedvgraph(y-1,1)=isempty(paths2y);%if the node is not connected it returns an empty cell, which this answer is 1
            end
            
            if sum(flag_connectedvgraph)==0
                break;
            end
        end
        for m=1:M
            for mm=1:M
                if k<=SliceNum(1,t) && m<=NumReqVMs(t,k) && mm<=NumReqVMs(t,k)
                    Vlink_adj(m,mm,t,k)=temp_vadj(m,mm);
                    if m==mm
                        Tau_max_vl(m,mm,t,k)=0;
                        %10^6 is assumed as infinite throughput (between VMs in one cloud node (Intra-cloud)%%%%%%%%%
                        Varpi_vl(m,mm,t,k)=10^6;
                    elseif m>mm
                        Varpi_vl(m,mm,t,k)=temp_vadj(m,mm)*ceil(10000+rand*100000);
                        if temp_vadj(m,mm)==1
                            Tau_max_vl(m,mm,t,k)=5+rand*9; %In order of 15 to 45 miliseconds
                        end
                    end
                else
                    Vlink_adj(m,mm,t,k)=-1;
                    Varpi_vl(m,mm,t,k)=-1;
                    Tau_max_vl(m,mm,t,k)=-1;
                end
            end
        end
        
        for m=1:M
            for mm=1:M
                if m<mm
                    if Vlink_adj(m,mm,t,k)==1
                        Tau_max_vl(m,mm,t,k)=Tau_max_vl(mm,m,t,k); %In order of 8 to 28 miliseconds
                        Varpi_vl(m,mm,t,k)=Varpi_vl(mm,m,t,k);
                    end
                end
            end
        end
    end
end
