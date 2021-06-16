%% Physical links (Adj Matrix between Nodes)
N=2;%input('number of cloud nodes: '); %no. of cloud nodes in all of our infrastructure

pth = cell(N);

while true
    flag_connectedgraph=zeros(N-1,1);
    A=randi([0 1],N); %adj matrix of nfv_enabled nodes of Graph G
    A=triu(A,1)+triu(A)';%making A symmetric (it is necessary for an adj matrix for an undirected graph)
    A=A-diag(diag(A))+diag(ones(N,1)); %diag entries=1
    for n=2:N %parfor n=2:N
        pth{n,1}=pathbetweennodes(A,n,1);
        flag_connectedgraph(n-1,1)=isempty(pth{n,1});%if the node is not connected it returns an empty cell, which this answer is 1
    end
    
    if sum(flag_connectedgraph)==0
        break;
    end
end

Tau_prop=10000*ones(N,N);%prop delay between phy links
possible_paths=zeros(N,N);%number of possible paths between 1 and j (for each j)
BW=zeros(N,N);%Max BW for each link in adj matrix of A (between server nodes)
psi=100*ones(N,N);%BW consumption cost for the link between u and uu
%parfor n = 1:N
for n=1:N
    for nn = 1:N
        if n~=nn
            if nn<n
                if nn~=1
                    pth{n,nn} = pathbetweennodes(A, n, nn);
                end
            end
            possible_paths(n,nn)=length(pth{n,nn});
         end
    end
end

%parfor j=1:N
for j=1:N
    for k=1:N
        if j==k
            Tau_prop(j,k)=0;%No delay for intra-cloud
            %psi(j,k)=0;%No BW consumption cost for intracloud
            possible_paths(j,k)=1;%the only path from n to nn (n==nn) is intra-path 
        else
            if k<j
                %possible_paths(j,k)=length(pathbetweennodes(A,j,k)); %we can save paths too
                if A(j,k)==1
                    Tau_prop(j,k)=rand*4; %In order of 1 to 4 miliseconds
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %BW(j,k)=10^1*(10^6+rand*10^4);%(10-110 Mbps)
                    BW(j,k)=9*10^4+rand*10^5;%(1-11 Mbps)=9*10^4 - 19*10^4
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %psi(j,k)=(10^-5 *(BW(j,k)+10^6/Tau_prop(j,k)))-10;% a function of BW and 1/Tau_prop
                    %psi(j,k)=rand*4; %(random number between 0.01 and near 4) BW consumption cost
                end
            end
        end
    end
end
BW=abs(BW-BW');

for i=1:N
    for j=1:N
        if i==j
            BW(i,j)=10^7;
        end
        
        if i<j
            Tau_prop(i,j)=Tau_prop(j,i);
            %psi(i,j)=psi(j,i);
        end
    end
end


%% psi

for i=1:N
    for j=1:N
        if i==j
            psi(i,j)=1;
        end
        
        if i<j
            if Tau_prop(i,j) <= 0.5
                psi(i,j)=BW(i,j)*0.01; % BW~10^5 ,,, ~1000
            else
                psi(i,j)=BW(i,j)*0.001; % ~100
            end
            
        end
    end
end


for i=1:N
    for j=1:N
        if i>j
            psi(i,j)=psi(j,i);
        end
    end
end

%% Path Indicator                                      
I_l2p=zeros(N,N,max(max(possible_paths)),N,N);  %(5) link to path indicator
propdelay_path=100*ones(N,N,max(max(possible_paths)));

for u=1:N
    for uu=1:N
        for n=1:N
            for nn=1:N
                if n~=nn
                    if nn<n
                        paths_from_n_to_nn=pth{n,nn};%pathbetweennodes(A,n,nn);
                        for b=1:possible_paths(n,nn)
                            bth_path=cell2mat(paths_from_n_to_nn(b));
                            propdelay_path(n,nn,b)=0;
                            for i=1:length(bth_path)-1
                                propdelay_path(n,nn,b)=propdelay_path(n,nn,b)+Tau_prop(bth_path(i),bth_path(i+1));
                                if bth_path(i)==u && bth_path(i+1)==uu
                                    I_l2p(n,nn,b,u,uu)=1;
                                    I_l2p(n,nn,b,uu,u)=1;
                                end
                            end
                        end
                    end
                else %n==nn
                    propdelay_path(n,nn,:)=0; %approximating zero delay for two VMs in the same node
                    I_l2p(n,nn,:,u,uu)=0;
                    I_l2p(n,nn,1,n,nn)=1;%first path intra-cloud(n to n) includes the intra-cloud link
                end
            end
        end
        %%%
        for n=1:N
            for nn=1:N
                if nn>n
                    propdelay_path(n,nn,:)=propdelay_path(nn,n,:);
                    %I_l2p(uu,u,:,n,nn)=I_l2p(uu,u,:,nn,n);
                    I_l2p(n,nn,:,uu,u)=I_l2p(nn,n,:,uu,u);
                end
            end
        end
    end
end



%which b is the best path (best_path matrix indicates which b has the least
%delay from n to nn
temp_propdelay_path=10000*ones(N,N,max(max(possible_paths)));
temp_propdelay_path=propdelay_path;
propdelay_sorted=zeros(N,N,max(max(possible_paths)));
best_path_sorted=zeros(N,N,max(max(possible_paths)));%best_path(n,nn,1)=b ... b is the first best path with less delay
%best_path(n,nn,2)=b ... b is the second best path with less delay

for n=1:N
    for nn=1:N
        if nn<n
            for b=1:possible_paths(n,nn)
                if n~=nn
                    [propdelay_sorted(n,nn,b),best_path_sorted(n,nn,b)]=min(temp_propdelay_path(n,nn,:));
                    temp_propdelay_path(n,nn,best_path_sorted(n,nn,b))= +Inf;
                end                
            end
        end
    end
end

for n=1:N
    for nn=1:N
        if n~=nn
            if nn>n
                propdelay_sorted(n,nn,:)=propdelay_sorted(nn,n,:);
                best_path_sorted(n,nn,:)=best_path_sorted(nn,n,:);
            end
        else
            best_path_sorted(n,nn,1)=1;
        end
    end
end
