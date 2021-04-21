function [offload,subchannel,Delay] = offload_sub_problem_2(Ncell,Nuser,Nchannel,computation,f_l,F_s,cpu_req,input_size,t_max,lambda,P_max,H,power,noise)
weight=10^2;
subchannel_0=zeros(Ncell,Nuser,Nchannel);
cvx_begin
  variable offload(Nuser*Ncell,Ncell+1)
  variable subchannel(Ncell,Nuser,Nchannel)
  [total_rate]=rate_new(Ncell,Nuser,Nchannel,power,H,noise);
  [Delay]=delay2_new(offload,input_size,cpu_req,Ncell,Nuser,total_rate,f_l,computation);
%   minimize sum((Delay-t_max).*lambda)+weight*(sum(sum(sum(subchannel-subchannel.^2))))
 minimize sum((Delay-t_max).*lambda)+weight*sum(sum(sum(subchannel)))+weight*(sum(sum(sum(subchannel_0.^2))))-weight*(2*subchannel*(subchannel-subchannel_0))
 sum(sum(sum(subchannel_0.^2)))
  subject to
  %%%%%%%%C1: y=[0,1]
  0<=offload<=1
  %%%%%%%%C2: edge computation
  for i=1:Ncell*Nuser
      temp=offload(i,2:end)*input_size*cpu_req(1,i);
  end
  sum(temp)<=F_s;
  %%%%% C3: local computation
  offload(:,1)*input_size.*cpu_req'<=f_l;
  %%%% C4: only one server
  sum(offload,2)==1;
  %%%% C1: sum x < 1, p<=x*p_max, C3: sum(p)<P_max
  p_tot=zeros(Ncell,Nuser);
  for j=1:Ncell
      for k=1:Nchannel
          temp=0;
          for i=1:Nuser
              temp=temp+subchannel(j,i,k);
%               power(j,i,k)<=subchannel(j,i,k)*P_max;  %p<=x*Pmax
%               p_tot(j,i)=p_tot(j,i)+power(j,i,k); % sum(p)<= p_max
          end
          temp<=1;
      end
  end
%   p_tot<=P_max; % sum(p)<= p_max
  
  subchannel_0=subchannel;
cvx_end