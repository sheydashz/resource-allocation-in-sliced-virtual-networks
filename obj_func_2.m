function [OBJ] = obj_func_2(X,z,H,say,input_size,offload,Alpha,Beta,F_s,P_max,t_max,lambda,services,Nservices,Nchannel,B,noise,Ncell,Nuser,cpu_req,f_l,service);
% save('params.mat')
load('augmented.mat','theta','delta','phi','zeta','eta','say')
[total_rate,update_z] = rate_restated_2(X,z,H,Ncell,Nuser,Nchannel);
OBJ=0;
%% Calculate delay
offloaded=sum(offload,2)-offload(:,1);
%%%%%%%%%%% trans delay
temp=offloaded*input_size;
trans_delay=temp./reshape(total_rate',Ncell*Nuser,1);
%%%%%%%%%%% backhaul delay + augmented for sum(p)-Pmax
obj2=0;
for j=1:Ncell
    indx=(j-1)*Nuser;
    for i=1:Nuser
        temp2(indx+i)=offloaded(indx+i)-offload(indx+i,j+1);
        index1=(j-1)*Nuser*Nchannel+(i-1)*Nchannel;
        index2=(j-1)*Nuser+i;
        obj2=obj2+theta(index2,1)+say*(sum(X(index1+1:index1+Nchannel))-P_max);
    end
end
back_delay=temp2'*input_size*5e-6;
obj2= 1/(2*say)*(max(0,obj2))^2;
%%%%%%%%%%%%%%%%%%%%%%%%%Computation Delay %%%%%%%%%%%%%%%
%%%%%%%%%%% local computation delay
T_l=(offload(:,1)*input_size.*cpu_req')./f_l;
%%%%%%%%%% edge computation delay %%%%%%%%%%%%%%%%%
T_E=zeros(Ncell*Nuser,1);
index0=2*Nuser*Ncell*Nchannel;
for j=1:Nuser
    index=index0+(j-1)*Ncell;
    temp=0;
    for i=1:Ncell
        indx=index+i;
        temp=offload(Nuser*(i-1)+j,i+1)*input_size*cpu_req(1,Nuser*(i-1)+j);
        T_E(Nuser*(i-1)+j,1)=T_E(Nuser*(i-1)+j,1)+temp/X(indx,1);
    end
end
Delay=T_E+trans_delay+T_l+back_delay;
obj1= sum((Delay-t_max).*lambda);
%%% sum(computation)-Beta*S
obj3=0;
comp=zeros(Nservices,1);
index0=Ncell*Nuser*Nchannel;
for k=1:Nservices
    for j=1:Nuser*Ncell
            index=index0+(j-1)*Ncell;
            comp(k,1)=comp(k,1)+services(k,j)*sum(X(index+1:index+Ncell));
    end
    obj3=obj3+delta(k,1)+say*(comp(k,1)-Beta(k)*F_s*Ncell);
end
obj3=1/(2*say)*(max(0,obj3))^2-sum(delta.^2);
%%% p<x*Pmax
obj6=0;
temp=0;
say1=100;
for j=1:Ncell
    for i=1:Nuser
        temp2=0;
        for k=1:Nchannel
            index1=(j-1)*Nuser*Nchannel+(i-1)*Nchannel+k;
            temp2=temp2+X(index1,1)*subchannel(j,i,k);
        end
     temp=temp+eta(j,i)+say1*(temp2-P_max);

    end
end
obj6= 1/(2*say1)*max(0,temp)^2-sum(sum(eta.^2));
%%% 
OBJ=obj1+obj2+obj3+obj4+obj6
save('X.mat')


            
        


