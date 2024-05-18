clc
clear all;
%% SVMD based adaptive clustering method
% general idea, view the clustering issue as a decomposing problem. each time extracts a component.
% the obstacle is to define the value of alf. alf should decrease as the
% decomposition progresses.
% SD defines the similarity(distance) between the original signals after each extraction. if after certain number
% of extraction, SD increases or flatens but not decreases(if sd is
% similarity,otherwise if sd is sample distance, it would decrease or flatten)
% the number of clustering is determinated, the following extraction is the refinement approach.
%written by chen wei
%@2021.1.28


%generate scatters of Gaussian distribution
cx1=-0.5;
cy1=1;
cx2=-1;
cy2=-1;
cx3=1;
cy3=1;
cx4=1;
cy4=-1;

mul=[cx1,cy1];
S1=[0.12 0;0 0.15];
data1=mvnrnd(mul, S1, 400);

mu2=[cx2 cy2];
S2=[0.12 0;0 0.15];
data2=mvnrnd(mu2,S2,320);

mu3=[cx3 cy3];
S3=[0.15 0;0 0.12];
data3=mvnrnd(mu3,S3,300);

mu4=[cx4 cy4];
S4=[0.1 0;0 0.1];
data4=mvnrnd(mu4,S4,250);


figure(1),plot(data1(:,1),data1(:, 2),'b+');hold on; 
plot(data2(:,1),data2(:,2),'r+');hold on
plot(data3(:,1),data3(:,2),'g+'); hold on
plot(data4(:,1),data4(:,2),'k+'); xlim([-3,3]),ylim([-3,3]); hold off
X=[data1(:,1); data2(:,1); data3(:,1);data4(:,1)];
Y=[data1(:,2); data2(:,2); data3(:,2);data4(:,2)];

% figure(2),plot(X,Y,'mo');hold off;
L=length(X);
Z=zeros(L,1);
mx1=10;
mx2=10;
mx3=10;
mx4=10;
%find peaks location,also can use findpeaks function
for i=1:L
    ds1=sqrt((X(i)-cx1)^2+(Y(i)-cy1)^2);
    ds2=sqrt((X(i)-cx2)^2+(Y(i)-cy2)^2);
    ds3=sqrt((X(i)-cx3)^2+(Y(i)-cy3)^2);
    ds4=sqrt((X(i)-cx4)^2+(Y(i)-cy4)^2);
    if mx1>ds1
        mx1=ds1;
        p1=i;
    end
    if mx2>ds2
        mx2=ds2;
        p2=i;
    end
    if mx3>ds3
        mx3=ds3;
        p3=i;
    end
    if mx4>ds4
        mx4=ds4;
        p4=i;
    end
end
D=ones(L,L)*10;  %10 is defined by user,can obtained from the samples. D is the distance matrix
for i=1:L
    for j=1:L
        if i~=j
            D(i,j)=sqrt((X(i)-X(j))^2+(Y(i)-Y(j))^2);
        end
    end
end
n=max(5,floor(L/50)); %n is defined by user, 
% generate the density. the inverse of the average value of the distance
% from the nearest n points is defined as the density of the target sample
for i=1:L
    [drS,pos]=sort(D(i,:),'ascend');
    md=mean(drS(1:n));
    Z(i)=1/md;
end
figure(2),plot3(X,Y,Z,'b.');hold on
figure(3),plot(data1(:,1),data1(:, 2),'b+');hold on;
plot(data2(:,1),data2(:,2),'r+');hold on
plot(data3(:,1),data3(:,2),'g+'); hold on
plot(data4(:,1),data4(:,2),'k+'); hold on,xlim([-3,3]),ylim([-3,3]);
figure(4),plot(data1(:,1),data1(:, 2),'b+');hold on;xlim([-3,3]),ylim([-3,3]);
plot(data2(:,1),data2(:,2),'r+');hold on
plot(data3(:,1),data3(:,2),'g+'); hold on
plot(data4(:,1),data4(:,2),'k+'); hold on
% SS(1)=sum(Z);
ZN=Z./norm(Z);
%% GVMD
u10=zeros(L,1);
u20=zeros(L,1);
u30=zeros(L,1);
u40=zeros(L,1);
u10(p1)=1; %any value does not matter
% u1=u10;
u20(p2)=1;
u30(p3)=1;
u40(p4)=1;
% u1=u10;
m=15;
k=0;
alf=15;
% beta=1;
        c1x=u10'*X/sum(u10);
        c1y=u10'*Y/sum(u10);
        c2x=u20'*X/sum(u20);
        c2y=u20'*Y/sum(u20);
        c3x=u30'*X/sum(u30);
        c3y=u30'*Y/sum(u30);
        c4x=u40'*X/sum(u40);
        c4y=u40'*Y/sum(u40);
%         u1c=Z-u10;
%         crx=u1c'*X/sum(u1c);
%         cry=u1c'*Y/sum(u1c);

%          sr1=sum(((X-crx).^2+(Y-cry).^2))
%          sd1=sum(u1c.^2)
%          st1=sum(((X-c1x).^2+(Y-c1y).^2))
%          alf= 0.8*sqrt(sd1^2/(st1*sr1))   %0.5-1.5,  
         % alf is defined emperically using sr(distance from the center of remaining components)
         % sd(power) st(distance from the center of extracting component)
while k<m
%      for i=1:L
%           u1(i)=(z
% %         u1(i)=Z(i)*(1+beta*((X(i)-crx)^2+(Y(i)-cry)^2))/(1+alf*((X(i)-c1x)^2+(Y(i)-c1y)^2)+beta*((X(i)-crx)^2+(Y(i)-cry)^2));
% 
%     end
      u1=(Z-u20-u30-u40)./(1+alf*((X-c1x).^2+(Y-c1y).^2));
      u2=(Z-u1-u30-u40)./(1+alf*((X-c2x).^2+(Y-c2y).^2));
      u3=(Z-u1-u2-u40)./(1+alf*((X-c3x).^2+(Y-c3y).^2));
      u4=(Z-u1-u2-u3)./(1+alf*((X-c4x).^2+(Y-c4y).^2));
        k=k+1;
         c1x=u1'*X/sum(u1);
         c1y=u1'*Y/sum(u1);
         c2x=u2'*X/sum(u2);
         c2y=u2'*Y/sum(u2);
         c3x=u3'*X/sum(u3);
         c3y=u3'*Y/sum(u3);    
         c4x=u4'*X/sum(u4);
         c4y=u4'*Y/sum(u4);  
%          u1c=Z-u1;
%          crx=u1c'*X/sum(u1c);
%          cry=u1c'*Y/sum(u1c);

% 
%     if sum(sum(abs(u10-u1)))<1e-2  %stop criterion
%        
%         break;
%     else
        
        u10=u1;
        u20=u2;
        u30=u3;
        u40=u4;
end
 
    clust=zeros(L,1);
    uf(1:4,1:L)=[u1';u2';u3';u4'];
    figure(3);
for i=1:L
    
    [vs,p]=max(uf(:,i)); %may use distance from each center of each component(use statistical distance may be better?) much better
    clust(i)=p;
    if p==1
        plot(X(i),Y(i),'bo');hold on;
    elseif p==2
        plot(X(i),Y(i),'ro');hold on;
    elseif p==3
        plot(X(i),Y(i),'go');hold on;
    else
        plot(X(i),Y(i),'ko');hold on;
    end
end

%% SGVMD
%% cycle I
u10=zeros(L,1);
u10(p1)=1; %any value does not matter
u1=u10;

% u1=u10;
m=15;
k=0;
% alf=15;
beta=1;
        c1x=u10'*X/sum(u10);
        c1y=u10'*Y/sum(u10);
        u1c=Z-u10;
        crx=u1c'*X/sum(u1c);
        cry=u1c'*Y/sum(u1c);

         sr1=sum(((X-crx).^2+(Y-cry).^2))
         sd1=sum(u1c.^2)
         st1=sum(((X-c1x).^2+(Y-c1y).^2))
         alf= 0.8*sqrt(sd1^2/(st1*sr1))   %0.5-1.5,  
         % alf is defined emperically using sr(distance from the center of remaining components)
         % sd(power) st(distance from the center of extracting component)
while k<m
     for i=1:L
      
        u1(i)=Z(i)*(1+beta*((X(i)-crx)^2+(Y(i)-cry)^2))/(1+alf*((X(i)-c1x)^2+(Y(i)-c1y)^2)+beta*((X(i)-crx)^2+(Y(i)-cry)^2));

    end
    k=k+1;
         c1x=u1'*X/sum(u1);
         c1y=u1'*Y/sum(u1);
         u1c=Z-u1;
         crx=u1c'*X/sum(u1c);
         cry=u1c'*Y/sum(u1c);


    if sum(sum(abs(u10-u1)))<1e-2  %stop criterion
       
        break;
    else
        
        u10=u1;
    end
    
end

% figure(3);plot3(X,Y,u1,'b.');hold on
XC(1)=c1x;
YC(1)=c1y;


%  hold off
%% cycle II
%% update Z
q2=sum((Z-u1))/sum(Z);  %caculate the decreasing factor of alf, related to the power ratio
Z=Z-u1;

% SS(2)=sum(Z);
ZN1=Z./norm(Z);
% SD(1)=dot(ZN,ZN1); %caculate the similarity
SD(1)=norm(ZN1-ZN); %caculate the distance
% 
% figure(7),plot(1,SD(1),'r.-','linewidth',2);hold on

%% Cycle II decomposition
u20=zeros(L,1);
u20(p2)=1;
u2=u20;


% u1=u10;
m=50;
k=0;
alf=q2*alf %update alf
beta=1;

while k<m

        c2x=u20'*X/sum(u20);
        c2y=u20'*Y/sum(u20);
        u2c=Z-u20;
        crx=u2c'*X/sum(u2c);
        cry=u2c'*Y/sum(u2c);


    for i=1:L

        u2(i)=Z(i)*(1+beta*((X(i)-crx)^2+(Y(i)-cry)^2))/(1+alf*((X(i)-c2x)^2+(Y(i)-c2y)^2)+beta*((X(i)-crx)^2+(Y(i)-cry)^2));
               
    end
    k=k+1;
    if sum(sum(abs(u20-u2)))<1e-2
        %     k;
        break;
    else
        
        u20=u2;
    end
    
end
% 
% figure(4);plot3(X,Y,u2,'b.');hold on
XC(2)=c2x;
YC(2)=c2y;


%%
%% cycle III
%% update Z
q3=sum((Z-u2))/sum(Z);%caculate the decreasing factor, related to the power ratio
Z=Z-u2;

ZN2=Z./norm(Z);
% SD(2)=dot(ZN,ZN2);%caculate the similarity
SD(2)=norm(ZN2-ZN); %caculate the distance
% figure(7),plot(2,SD(2),'r.-','linewidth',2);hold on

%% Cycle III decomposition
u30=zeros(L,1);
u30(p3)=1;
u3=u30;


% u1=u10;
m=50;
k=0;
alf=q3*alf%update alf
beta=1;

while k<m
        c3x=u30'*X/sum(u30);
        c3y=u30'*Y/sum(u30);
        u3c=Z-u30;
        crx=u3c'*X/sum(u3c);
        cry=u3c'*Y/sum(u3c);
    for i=1:L
  
         u3(i)=Z(i)*(1+beta*((X(i)-crx)^2+(Y(i)-cry)^2))/(1+alf*((X(i)-c3x)^2+(Y(i)-c3y)^2)+beta*((X(i)-crx)^2+(Y(i)-cry)^2));
              
    end
    k=k+1;
    if sum(sum(abs(u30-u3)))<1e-4
        %     k;
        break;
    else
        
        u30=u3;
    end
    
end

% figure(5);plot3(X,Y,u3,'b.');hold on
XC(3)=c3x;
YC(3)=c3y;
%%
%% cycle IV
%% update Z
q4=sum((Z-u3))/sum(Z);%caculate the decreasing factor, related to the power ratio
Z=Z-u3;


ZN3=Z./norm(Z);
% SD(3)=dot(ZN,ZN3); %caculate the similarity
SD(3)=norm(ZN3-ZN); %caculate the distance
% figure(7),plot(3,SD(3),'r.-','linewidth',2);hold on

%% Cycle IV decomposition
u40=zeros(L,1);
u40(p4)=1;
u4=u40;


% u1=u10;
m=50;
k=0;
alf=q4*alf; %update alf
beta=1;

while k<m
        c4x=u40'*X/sum(u40);
        c4y=u40'*Y/sum(u40);
        u4c=Z-u40;
        crx=u4c'*X/sum(u4c);
        cry=u4c'*Y/sum(u4c);
    
    for i=1:L
  
       
        u4(i)=Z(i)*(1+beta*((X(i)-crx)^2+(Y(i)-cry)^2))/(1+alf*((X(i)-c4x)^2+(Y(i)-c4y)^2)+beta*((X(i)-crx)^2+(Y(i)-cry)^2));
        
        
    end
    k=k+1;
    if sum(sum(abs(u40-u4)))<1e-4
        %     k;
        break;
    else
        
        u40=u4;
    end
    
end
 
% figure(6);plot3(X,Y,u4,'b.');hold on
XC(4)=c4x;
YC(4)=c4y;
%% update Z
q11=sum(Z-u4)/sum(Z);
Z=Z-u4;


ZN4=Z./norm(Z);
% SD(4)=dot(ZN,ZN4);%caculate the similarity
SD(4)=norm(ZN4-ZN); %caculate the distance
% figure(7),plot(4, SD(4),'r.-','linewidth',2);hold on  %if the similarity flatten or increase(sample distance flatten or decrease), first loop extraction finished.
% figure(2);plot3(X,Y,Z,'r.');hold off
%% if SD increase or flatten(decrease or flatten is sd is sample distance), enter into loop two
%% cycle I
u10=zeros(L,1);
u10(p1)=1;
u11=u10;


% u1=u10;
m=50;
k=0;
alf=alf*q11
beta=1;

while k<m
        c1x=u10'*X/sum(u10);
        c1y=u10'*Y/sum(u10);
        urc=Z-u10;
        crx=urc'*X/sum(urc);
        cry=urc'*Y/sum(urc);
    for i=1:L
 
       
        u11(i)=Z(i)*(1+beta*((X(i)-crx)^2+(Y(i)-cry)^2))/(1+alf*((X(i)-c1x)^2+(Y(i)-c1y)^2)+beta*((X(i)-crx)^2+(Y(i)-cry)^2));
        
        
    end
    k=k+1;
    if sum(sum(abs(u10-u11)))<1e-4
        %     k;
        break;
    else
        
        u10=u11;
    end
    
end

XC(5)=c1x;
YC(5)=c1y;
%% update Z
q21=sum(Z-u11)/sum(Z);
Z=Z-u11;
% SS(6)=sum(Z);
ZN5=Z./norm(Z);
% SD(5)=dot(ZN,ZN5); %caculate the similarity
SD(5)=norm(ZN5-ZN); %caculate the distance
% figure(7),plot(5,SD(5),'r.-','linewidth',2);hold on
% figure(10);plot3(X,Y,Z,'.');
%% cycle II
u20=zeros(L,1);
u20(p2)=10;
u21=u20;


% u1=u10;
m=50;
k=0;
alf=alf*q21
beta=1;

while k<m
        c2x=u20'*X/sum(u20);
        c2y=u20'*Y/sum(u20);
        urc=Z-u20;
        crx=urc'*X/sum(urc);
        cry=urc'*Y/sum(urc);
    for i=1:L
   
        %  dr=sqrt((xs(i)-c1x)^2+(ys(i)-c1y)^2);
        u21(i)=Z(i)*(1+beta*((X(i)-crx)^2+(Y(i)-cry)^2))/(1+alf*((X(i)-c2x)^2+(Y(i)-c2y)^2)+beta*((X(i)-crx)^2+(Y(i)-cry)^2));
        
        
    end
    k=k+1;
    if sum(sum(abs(u20-u21)))<1e-4
        %     k;
        break;
    else
        
        u20=u21;
    end
    
end
XC(6)=c2x;
YC(6)=c2y;

%figure(4);plot3(X,Y,u21,'r.');hold off
%% update Z
q31=sum(Z-u21)/sum(Z);
Z=Z-u21;
% SS(7)=sum(Z);
ZN6=Z./norm(Z);
% SD(6)=dot(ZN,ZN6);%caculate the similarity
SD(6)=norm(ZN6-ZN); %caculate the distance
% figure(7),plot(6,SD(6),'r.-','linewidth',2);hold on
% figure(12);plot3(X,Y,Z,'.');
%% cycle III
u30=zeros(L,1);
u30(p3)=10;
u31=u30;


% u1=u10;
m=50;
k=0;
alf=q31*alf
beta=1;

while k<m
        c3x=u30'*X/sum(u30);
        c3y=u30'*Y/sum(u30);
        urc=Z-u30;
        crx=urc'*X/sum(urc);
        cry=urc'*Y/sum(urc);
    for i=1:L

        %  dr=sqrt((xs(i)-c1x)^2+(ys(i)-c1y)^2);
        u31(i)=Z(i)*(1+beta*((X(i)-crx)^2+(Y(i)-cry)^2))/(1+alf*((X(i)-c3x)^2+(Y(i)-c3y)^2)+beta*((X(i)-crx)^2+(Y(i)-cry)^2));
        
        
    end
    k=k+1;
    if sum(sum(abs(u30-u31)))<1e-4
        %     k;
        break;
    else
        
        u30=u31;
    end
    
end

XC(7)=c3x;
YC(7)=c3y;
%figure(5);plot3(X,Y,u31,'r.');hold off
%%
%% update Z
q41=sum(Z-u31)/sum(Z);
Z=Z-u31;
% SS(8)=sum(Z);
ZN7=Z./norm(Z);
% SD(7)=dot(ZN,ZN7); %caculate the similarity
SD(7)=norm(ZN7-ZN); %caculate the distance
% figure(7),plot(7,SD(7),'r.-','linewidth',2);hold on
% figure(14);plot3(X,Y,Z,'.');
%% cycle III
u40=zeros(L,1);
u40(p4)=10;
u41=u40;


% u1=u10;
m=50;
k=0;
alf=q41*alf
beta=1;

while k<m
        c4x=u40'*X/sum(u40);
        c4y=u40'*Y/sum(u40);
        urc=Z-u40;
        crx=urc'*X/sum(urc);
        cry=urc'*Y/sum(urc);
    for i=1:L
 
        %  dr=sqrt((xs(i)-c1x)^2+(ys(i)-c1y)^2);
        u41(i)=Z(i)*(1+beta*((X(i)-crx)^2+(Y(i)-cry)^2))/(1+alf*((X(i)-c4x)^2+(Y(i)-c4y)^2)+beta*((X(i)-crx)^2+(Y(i)-cry)^2));
        
        
    end
    k=k+1;
    if sum(sum(abs(u40-u41)))<1e-4
        %     k;
        break;
    else
        
        u40=u41;
    end
    
end

XC(8)=c4x;
YC(8)=c4y;

%% update Z
Z=Z-u41;
% SS(9)=sum(Z);
ZN8=Z./norm(Z);
% SD(8)=dot(ZN,ZN8); %caculate the similarity
SD(8)=norm(ZN8-ZN); %caculate the distance
% figure(16);plot3(X,Y,Z,'.');
% figure(7),plot(1:8,SD(1:8),'r.-','linewidth',2); %after one loop of full extraction, the similarity would increase or flatten instead of decrease(the sample distance decrease or flatten), so the number of clusters is determinated
% xlabel('Times','FontSize',15);ylabel('Normalized Similarity','FontSize',15),grid on,hold off
%% final decision, to add similar components togeter, here we consider components with close peak positions belong to one cluster


ub=[u11';u21';u31';u41'];
ut=[u1';u2';u3';u4'];
cxb=zeros(4,1);
cyb=zeros(4,1);
cxt=zeros(4,1);
cyt=zeros(4,1);
for i=1:4
%     cxb(i)=ub(i,:)*X/sum(ub(i,:));
%      cyb(i)=ub(i,:)*Y/sum(ub(i,:));
      [vv,pos]=max(ub(i,:)); %find max
      cxb(i)=X(pos);
      cyb(i)=Y(pos);
      [vv,pos]=max(ut(i,:));
      cxt(i)=X(pos);
      cyt(i)=Y(pos);
%      cxt(i)=ut(i,:)*X/sum(ut(i,:));
%      cyt(i)=ut(i,:)*Y/sum(ut(i,:));
end
uf=ut;
clear vp
% caculate the distance
for i=1:4
    for j=1:4
%         vp(j)=dot(ub(i,:),ut(j,:))/(norm(ub(i,:))*norm(ut(j,:)));
       vp(j)=(cxb(i)-cxt(j))^2+(cyb(i)-cyt(j))^2;  %caculate the distance between the peaks
    end
    [v,p]=min(vp); % find the closest one
    
    
%     figure(2+p),plot3(X,Y,ub(i,:)','r.');hold on  %plot on the same figure
    uf(p,:)= uf(p,:)+ub(i,:);
end

for i=1:4
%     figure(2+i),hold off;
end



figure(4);
clust=zeros(L,1);
for i=1:L
    
    [vs,p]=max(uf(:,i)); %may use distance from each center of each component(use statistical distance may be better?) much better
    clust(i)=p;
    if p==1
        plot(X(i),Y(i),'bo');hold on;
    elseif p==2
        plot(X(i),Y(i),'ro');hold on;
    elseif p==3
        plot(X(i),Y(i),'go');hold on;
    else
        plot(X(i),Y(i),'ko');hold on;
    end
end
xlim([-3,3]);ylim([-3,3]);

hold off
% figure(8),plot(1:9,SS,'b*-');
 figure(5),plot(1:8,SD(1:8),'r*-','linewidth',2); %after one loop of full extraction, the similarity would increase or flatten instead of decrease(the sample distance decrease or flatten), so the number of clusters is determinated
xlabel('Times','FontSize',15);ylabel('Normalized Similarity','FontSize',15),grid on,hold off
