function [total_rate,update_z] = rate_restated_2(X,z,H,Ncell,Nuser,Nchannel);
% load('params.mat')
noise=1e-10;
B=20*1e3;
total_rate=zeros(Ncell,Nuser);
update_z=zeros(Ncell,Nuser,Nchannel);

for i=1:Ncell
    for j=1:Nuser
        for k=1:Nchannel
            I=0;
            for ii=1:Ncell
                if (ii ~= i)
                for jj=1:Nuser
                    index=(ii-1)*Nuser*Nchannel+(jj-1)*Nchannel+k;
                    index2=(ii-1)*Nuser+jj;
                    I=I+X(index,1)*H(i,index2,k);
                end
                end
            end
        index=(i-1)*Nuser*Nchannel+(j-1)*Nchannel+k;
        index2=(i-1)*Nuser+j;
        num=sqrt(H(i,index2,k)*X(index,1));
        denum=noise+I;
        total_rate(i,j) = total_rate(i,j)+B*log2(1+2*z(i,j,k)*num-z(i,j,k)^2*denum);
        update_z(i,j,k)=num/denum;

        end
    end
end
