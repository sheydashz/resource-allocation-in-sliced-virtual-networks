clear all
clc

%% Parameters %%%
Nuser=2;
Nservices=2;
Nchannel=4;
Ncell=2;
radius=150;
%%% sevices %%%
temp=1:Nservices;
Lambda=[3,1];
Beta=[0.6,0.4];
Alpha=[0.7,0.3];
for i=1:Ncell
    service(i,:)= randsample(temp,Nuser,true);
end
lambda=Lambda(service);
lambda=reshape(lambda',Nuser*Ncell,1);
service=reshape(service',Nuser*Ncell,1);
T_max=[5*1e-2,1*1e-1,5];
for i=1:Nuser*Ncell
    t_max(i,1)=T_max(find(Lambda==lambda(i)));
end

temp2=reshape(service',Nuser*Ncell,1);
for k=1:Nservices
    for i=1:Nuser*Ncell
        if temp2(i,1)==k
            services(k,i)=1;
        end
    end
end
f_l=0.6*10^9;
F_s=16*1e9;
noise=1e-10;
P_max=0.2;

input_size=0.1*1e6;
A=[1,2,3];
cpu_req=reshape(datasample(A,Nuser*Ncell)*1500,1,Nuser*Ncell);

power=rand(Ncell,Nuser,Nchannel);
computation=ones(Ncell*Nuser,Ncell)*(F_s./Nuser);
%% Path Gain %%%
H=pathgain(Ncell,Nuser,Nchannel,radius);
[offload,subchannel,Delay] = offload_sub_problem_2(Ncell,Nuser,Nchannel,computation,f_l,F_s,cpu_req,input_size,t_max,lambda,P_max,H,power,noise);
[power,channel,computation]= second_sub_problem_2(offload,input_size,H,Nuser,Ncell,Nchannel,lambda,cpu_req,f_l,F_s,noise,B,Nservices,P_max,Alpha,Beta,t_max,services);




