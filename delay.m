function [Delay]= delay2_new(offload,input_size,cpu_req,Ncell,Nuser,total_rate,f_l,computation);
offloaded=sum(offload,2)-offload(:,1);
%%%%%%%%%%% trans delay
temp=offloaded*input_size;
trans_delay=temp./reshape(total_rate',Ncell*Nuser,1);
%%%%%%%%%%% backhaul delay
for j=1:Ncell
    indx=(j-1)*Nuser;
    for i=1:Nuser
        temp2(indx+i)=offloaded(indx+i)-offload(indx+i,j+1);
    end
end
back_delay=temp2'*input_size*5e-6;
%%%%%%%%%%%%%%%%%%%%%%%%%Computation Delay %%%%%%%%%%%%%%%
%%%%%%%%%%% local computation delay
T_l=(offload(:,1)*input_size.*cpu_req')./f_l;
%%%%%%%%%% edge computation delay %%%%%%%%%%%%%%%%%
T_E=cvx(zeros(Ncell*Nuser,1));
for j=1:Nuser
    temp=0;
    for i=1:Ncell
        indx=(i-1)*Nuser+j;
        temp=offload(indx,i+1)*input_size*cpu_req(1,indx);
        T_E(indx,1)=T_E(indx,1)+temp/computation(indx,i);
    end
end
Delay=T_E+trans_delay+T_l+back_delay;
        