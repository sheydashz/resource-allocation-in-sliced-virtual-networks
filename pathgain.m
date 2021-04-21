function [H]= pathgain(Ncell,Nuser,Nchannel,radius);

xBS=[radius,radius*3,radius,radius*3];
yBS=[radius,radius,-radius,-radius];
x=zeros(1,Ncell*Nuser);
y=zeros(1,Ncell*Nuser);

Xsep=[0,22,0,22];
Ysep=[0,0,22,22];

path_comp=zeros(Ncell,Ncell*Nuser);
H=zeros(Ncell,Ncell*Nuser,Nchannel);

%% coordinates
for i=1:Ncell    
    t = 2*pi*rand(Nuser,1);
    r = radius*sqrt(rand(Nuser,1));
    
    indx=Nuser*(i-1);
    x(1,indx+1:indx+Nuser) = xBS(i) + r.*cos(t)+Xsep(i);
    y(1,indx+1:indx+Nuser) = yBS(i) + r.*sin(t)+Ysep(i);
end

%% distance
dist=zeros(Ncell,Ncell*Nuser);
for i=1:Ncell
    dist(i,:)= sqrt((xBS(i) - x).^2 + (yBS(i) - y).^2)/60;
    temp=34+40*log(dist(i,:));
    temp= 10.^(-temp/10);
    path_comp(i,:)=sqrt(temp/2);
end


%% Pathgain in each channel
for i=1:Ncell
    for j=1:Ncell*Nuser
            H(i,j,:)=path_comp(i,j).*(rand(1,Nchannel)./dist(i,j));
    end
end
end