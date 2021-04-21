function [power,channels,computation]= second_sub_problem(offload,input_size,H,Nuser,Ncell,Nchannel,lambda,cpu_req,f_l,F_s,noise,B,Nservices,P_max,Alpha,Beta,t_max,services);

say=32;
z=ones(Ncell,Nuser,Nchannel)*10e1;
power=ones(Ncell*Nuser*Nchannel,1)*1e-1*P_max./Nchannel;
lb=zeros(size(power));
ub=ones(size(power))*P_max;
computation=rand(Ncell*Nuser*Ncell,1)*F_s/100;
lb=vertcat(lb,zeros(size(computation)));
ub=vertcat(ub,ones(Ncell^2*Nuser,1)*F_s);
X=vertcat(power,computation);
z0=0;
power0=zeros(Ncell,Nuser,Nchannel);
computation0=zeros(Ncell*Nuser,Ncell);
while (abs(z-z0)>0.0001)
    %%% augment lagrangian multipliers
    theta=ones(Nuser*Ncell,1)*1e-2; % for sum(p)-Pmax
    delta=ones(Nservices,1)*1e-4; % for sum(f)-Beta*S
%     phi=ones(Nchannel,Ncell)*1e-4; % for x_{u,n}<= 1
%     zeta=ones(Ncell*Nuser,Nchannel)*1e-3; % for x-x^2
    eta=ones(Ncell,Nuser)*1e-2; % for sum p*x <Pmax 
    converge=0;
    while(converge==0)
    save('augmented.mat')
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    nonlcon = [];
    options = optimset('TolFun',1e-4,'TolCon',1e-1,'TolX',1e-8,'MaxFunEvals',10000,'MaxIter',10000);
    [OBJ]=@(X)obj_func_2(X,z,H,say,input_size,offload,Alpha,Beta,F_s,P_max,t_max,lambda,services,Nservices,Nchannel,B,noise,Ncell,Nuser,cpu_req,f_l,services);
    [OBJ] = fmincon(OBJ,X,A,b,Aeq,beq,lb,ub,nonlcon,options);
    load('X.mat')
    power=zeros(Ncell,Nuser,Nchannel);
    channels=zeros(Ncell,Nuser,Nchannel);
    computation=zeros(Ncell*Nuser,Ncell);
    for j=1:Ncell
        for i=1:Nuser
            for k=1:Nchannel
                index=(j-1)*Nuser*Nchannel+(i-1)*Nchannel+k;
                power(j,i,k)=X(index,1);
                channels(j,i,k)=X(index+Nuser*Ncell*Nchannel,1);
                index=(j-1)*Nuser+i;
                index2=2*Nuser*Ncell*Nchannel+(i-1)*Ncell+j;
                computation(index,j)=X(index2,1);
            end
        end
    end
    theta0=theta;
    delta0=delta;
    phi0=phi;
    zeta0=zeta;
    eta0=eta;
    [theta,delta,phi,zeta,eta]=update(theta0,delta0,phi0,zeta0,eta0,P_max,Beta,Alpha,F_s,Ncell,Nchannel,Nuser,Nservices,services,say,channels,power,computation);
%     if sum(sum(sum(abs(power-power0))))+sum(sum(sum(abs(channels0-channels))))+sum(sum(abs(computation-computation0)))<1e-4
    if sum(abs(theta-theta0))+sum(abs(delta-delta0))+sum(sum(abs(phi-phi0)))+sum(sum(sum(abs(zeta-zeta0))))+sum(sum(sum(abs(eta-eta0))))<1e-2
        coverge=1
    end
    power0=power;
    channels0=channels;
    computation0=computation;
    end
   z0=z;
   z=update_z;
end