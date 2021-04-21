function [theta,delta,phi,zeta,eta]=update(theta0,delta0,phi0,zeta0,eta0,P_max,Beta,Alpha,F_s,Ncell,Nchannel,Nuser,Nservices,services,say,channels,power,computation)


%theta: sum p<= pmax
theta=max(0,theta0+say*(reshape(sum(power,3)',Ncell*Nuser,1)-P_max));
% delta: sum f<= beta*S
delta=zeros(Nservices,1);
for k=1:Nservices
    temp=0;
    for j=1:Nuser*Ncell
       temp=temp+services(k,j)*sum(computation(j,:),2);
    end
    delta(k,1)=max(0,delta0(k,1)+say*(temp-Ncell*F_s*Beta(k)));
end
% and eta: p<=x*pmax
eta=zeros(Ncell,Nuser,Nchannel);
for j=1:Ncell
    for k=1:Nchannel
        temp=sum(channels(j,:,k),2);
        for n=1:Nuser
            eta(j,n,k)=max(0,eta0(j,n,k)+say*(power(j,n,k)*subchannel(j,n,k)-P_max));
        end
    end
end



