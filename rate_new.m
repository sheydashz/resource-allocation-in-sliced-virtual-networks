function [total_rate]=rate_new(Ncell,Nuser,Nchannel,power,H,noise)
r=zeros(Ncell,Nuser,Nchannel);
total_rate=zeros(Ncell,Nuser);
B=20*1e3;

for cell=1:Ncell
    for i=1:Nuser
        for k=1:Nchannel
            inter=0;
            for cell2=1:Ncell
                if cell2 ~= cell
                    for i2=1:Nuser
                        inter=inter+power(cell2,i2,k)*H(cell,i2,k);
                    end
                end
            end
            sinr=power(cell,i,k)*H(cell,i,k)/(inter+noise);
            r(cell,i,k)=B*log2(1+sinr);
            total_rate(cell,i)=total_rate(cell,i)+r(cell,i,k);
        end
    end
end
end
