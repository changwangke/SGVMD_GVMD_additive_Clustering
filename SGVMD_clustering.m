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
% cx1=-0.5;
% cy1=1;
% cx2=-0.8;
% cy2=-0.8;
% cx3=1.2;
% cy3=0.7;
% cx4=1;
% cy4=-1;
cx1=0;
cy1=1.3;
cx2=-0.8;
cy2=-0.8;
cx3=1.2;
cy3=0.7;
cx4=1;
cy4=-1;

mul=[cx1,cy1];
S1=[0.12 0.08;0.08 0.12];
mu2=[cx2 cy2];
S2=[0.12 -0.06;-0.06 0.3];
mu3=[cx3 cy3];
S3=[0.25 0.1;0.1 0.10];
mu4=[cx4 cy4];
S4=[0.15 0;0 0.15];
% mul=[cx1,cy1];
% S1=[0.12 -0.05;-0.05 0.10];
data1=mvnrnd(mul, S1, 400);
% 
% mu2=[cx2 cy2];
% S2=[0.12 -0.06;-0.06 0.3];

data2=mvnrnd(mu2,S2,420);

% mu3=[cx3 cy3];
% S3=[0.25 0.1;0.1 0.15];

data3=mvnrnd(mu3,S3,380);

% mu4=[cx4 cy4];
% S4=[0.25 0;0 0.1];

data4=mvnrnd(mu4,S4,350);
LLs=[400,420,380,350];

figure(1),plot(data1(:,1),data1(:, 2),'b+');hold on; 
plot(data2(:,1),data2(:,2),'r+');hold on
plot(data3(:,1),data3(:,2),'g+'); hold on
plot(data4(:,1),data4(:,2),'k+'); xlim([-3,3]),ylim([-3,3]); hold on

X=[data1(:,1); data2(:,1); data3(:,1);data4(:,1)];
Y=[data1(:,2); data2(:,2); data3(:,2);data4(:,2)];

dataTot{1}=data1;
dataTot{2}=data2;
dataTot{3}=data3;
dataTot{4}=data4;

% figure(2),plot(X,Y,'mo');hold off;
L=length(X);
datafile=fopen("data.txt","W");
for i=1:L
fprintf(datafile,"%f %f\n",X(i),Y(i));
end
fclose(datafile);
Satus=zeros(L,1);
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

[X(p1),Y(p1)];
[X(p2),Y(p2)];
[X(p3),Y(p3)];
[X(p4),Y(p4)];
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
% figure(3),plot(data1(:,1),data1(:, 2),'b+');hold on;
% plot(data2(:,1),data2(:,2),'r+');hold on
% plot(data3(:,1),data3(:,2),'g+'); hold on
% plot(data4(:,1),data4(:,2),'k+'); hold on,xlim([-3,3]),ylim([-3,3]);
figure(3),plot(data1(:,1),data1(:, 2),'b+');hold on;xlim([-3,3]),ylim([-3,3]);
plot(data2(:,1),data2(:,2),'r+');hold on
plot(data3(:,1),data3(:,2),'g+'); hold on
plot(data4(:,1),data4(:,2),'k+'); hold on
% SS(1)=sum(Z);
ZN=Z./norm(Z);
Z0=Z;

%%
%% cycle I s1
correctnessS=zeros(1,4);


u10=zeros(L,1);
u10(p2)=1; %any value does not matter
u1=u10;

% u1=u10;
m=10;
k=0;
% alf=15;
beta=1;
        c1x=u10'*X/sum(u10);
        c1y=u10'*Y/sum(u10);
        u1c=Z-u10;
        crx=u1c'*X/sum(u1c);
        cry=u1c'*Y/sum(u1c);

        sig=[0.1,0;0,0.08];
        sigr=[0.1,0;0,0.08];
        % alf= 0.8*sqrt(sd1^2/(st1*sr1))   
        %0.5-1.5,
        alf=2.1;
         % alf is defined emperically using sr(distance from the center of remaining components)
         % sd(power) st(distance from the center of extracting component)
      error=rand(2,2)*1e-4;
      error=(error+error')/2;
while k<m

     for i=1:L      
        u1(i)=Z(i)*(1+beta*([X(i)-crx,Y(i)-cry]*pinv(sigr+error)*[X(i)-crx,Y(i)-cry]'))/(1+alf*([X(i)-c1x,Y(i)-c1y]*pinv(sig+error)*[X(i)-c1x,Y(i)-c1y]')+beta*([X(i)-crx,Y(i)-cry]*pinv(sigr+error)*[X(i)-crx,Y(i)-cry]'));
     end
 k
    k=k+1;
         c1x=u1'*X/sum(u1);
         c1y=u1'*Y/sum(u1);
%          u1c=Z-u1;
%          crx=u1c'*X/sum(u1c);
%          cry=u1c'*Y/sum(u1c);
clear totpar sig sigr cxy crxy
         totpar=SigCal2(u1,Z,X,Y,c1x,c1y);
         sig=totpar{1};
         sigr=totpar{2};
         cxy=totpar{3};
         crxy=totpar{4};
         c1x=cxy(1);
         c1y=cxy(2);
         crx=crxy(1);
         cry=crxy(2);      
    if sum(sum(abs(u10-u1)))<1e-2  %stop criterion
       
        break;
    else
        
        u10=u1;
    end
    
end
k
 figure(5);plot3(X,Y,u1,'b.');hold on
XC(1)=c1x;
YC(1)=c1y;

%[sd(1),pos1]=simNorm(u1,Z0);
%[sd(1),poss(1)]=simDist([c1x,c1y],[mean(X),mean(Y)]);
%belongs{1}=compBel(u1,Z-u1);
% belongs{1}=compBelDis2(u1,Z-u1,c1x,c1y,sig,crx,cry,sigr,L,X,Y,Status);
% Satus(belongs{1})=1;
belongs{1}=compBelAmpF(u1,c1x,c1y,sig,X,Y);
%  hold off
%% cycle I s2
%% update Z
q2=sum((Z-u1))/sum(Z);  %caculate the decreasing factor of alf, related to the power ratio
Z=Z-u1;


% SS(2)=sum(Z);
ZN1=Z./norm(Z);
% SD(1)=dot(ZN,ZN1); %caculate the similarity
SD(1)=norm(ZN1-ZN); %caculate the distance
%%
% figure(7),plot(1,SD(1),'r.-','linewidth',2);hold on

% Cycle II decomposition
u20=zeros(L,1);
u20(p4)=1;
u2=u20;


% u1=u10;
m=10;
k=0;
%alf=q2*alf %update alf
beta=1;
        c2x=u20'*X/sum(u20);
        c2y=u20'*Y/sum(u20);
        u2r=Z-u10;
        crx=u2r'*X/sum(u2r);
        cry=u2r'*Y/sum(u2r);

        sig=[0.1,0;0,0.08];
        sigr=[0.1,0;0,0.08];
        % alf= 0.8*sqrt(sd1^2/(st1*sr1))   
        %0.5-1.5,
        alf=2.3;
         % alf is defined emperically using sr(distance from the center of remaining components)
         % sd(power) st(distance from the center of extracting component)
      error=rand(2,2)*1e-4;
      error=(error+error')/2;

while k<m

    for i=1:L

         u2(i)=Z(i)*(1+beta*([X(i)-crx,Y(i)-cry]*pinv(sigr+error)*[X(i)-crx,Y(i)-cry]'))/(1+alf*([X(i)-c2x,Y(i)-c2y]*pinv(sig+error)*[X(i)-c2x,Y(i)-c2y]')+beta*([X(i)-crx,Y(i)-cry]*pinv(sigr+error)*[X(i)-crx,Y(i)-cry]'));
               
    end
    k=k+1;
          c2x=u2'*X/sum(u2)
          c2y=u2'*Y/sum(u2)
%          u1c=Z-u1;
%          crx=u1c'*X/sum(u1c);
%          cry=u1c'*Y/sum(u1c);
clear totpar sig sigr cxy crxy
         totpar=SigCal2(u2,Z,X,Y,c2x,c2y);
         sig=totpar{1};
         sigr=totpar{2};
         cxy=totpar{3};
         crxy=totpar{4};
         c2x=cxy(1);
         c2y=cxy(2);
         crx=crxy(1);
         cry=crxy(2);      
    if sum(sum(abs(u20-u2)))<1e-2
        %     k;
        break;
    else
        
        u20=u2;
    end
    
end
% 
 figure(5);plot3(X,Y,u2,'r.');hold on
XC(2)=c2x;
YC(2)=c2y;
%[sd(2),pos2]=simNorm(u2,[Z0,u1]);
%[sd(2),pos2]=simNorm(u2,u1);
%[sd(2),poss(2)]=simDist([c2x,c2y],[c1x,c1y]);
%belongs{2}=compBel(u2,Z-u2);
% belongs{2}=compBelDis2(u2,Z-u2,c2x,c2y,sig,crx,cry,sigr,L,X,Y,Satus);
% Satus(belongs{2})=1;
belongs{2}=compBelAmpF(u2,c2x,c2y,sig,X,Y);
[maxC(1),poss(1)]=simCount(belongs{2},belongs,1);

%% cycle I s3

%% update Z
q2=sum((Z-u2))/sum(Z);  %caculate the decreasing factor of alf, related to the power ratio
Z=Z-u2;

% SS(2)=sum(Z);
ZN1=Z./norm(Z);
% SD(1)=dot(ZN,ZN1); %caculate the similarity
SD(2)=norm(ZN1-ZN); %caculate the distance
%% 
% figure(7),plot(1,SD(1),'r.-','linewidth',2);hold on

% Cycle II decomposition
u30=zeros(L,1);
u30(p1)=1;
u3=u30;


% u1=u10;
m=10;
k=0;
%alf=q2*alf %update alf
beta=1;
        c3x=u30'*X/sum(u30);
        c3y=u30'*Y/sum(u30);
        u3r=Z-u30;
        crx=u3r'*X/sum(u3r);
        cry=u3r'*Y/sum(u3r);

        sig=[0.1,0;0,0.08];
        sigr=[0.1,0;0,0.08];
        % alf= 0.8*sqrt(sd1^2/(st1*sr1))   
        %0.5-1.5,
        alf=0.9;
         % alf is defined emperically using sr(distance from the center of remaining components)
         % sd(power) st(distance from the center of extracting component)
      error=rand(2,2)*1e-4;
      error=(error+error')/2;

while k<m

    for i=1:L

         u3(i)=Z(i)*(1+beta*([X(i)-crx,Y(i)-cry]*pinv(sigr+error)*[X(i)-crx,Y(i)-cry]'))/(1+alf*([X(i)-c3x,Y(i)-c3y]*pinv(sig+error)*[X(i)-c3x,Y(i)-c3y]')+beta*([X(i)-crx,Y(i)-cry]*pinv(sigr+error)*[X(i)-crx,Y(i)-cry]'));
               
    end
    k=k+1;
          c3x=u3'*X/sum(u3)
          c3y=u2'*Y/sum(u3)
%          u1c=Z-u1;
%          crx=u1c'*X/sum(u1c);
%          cry=u1c'*Y/sum(u1c);
clear totpar sig sigr cxy crxy
         totpar=SigCal2(u3,Z,X,Y,c3x,c3y);
         sig=totpar{1};
         sigr=totpar{2};
         cxy=totpar{3};
         crxy=totpar{4};
         c3x=cxy(1);
         c3y=cxy(2);
         crx=crxy(1);
         cry=crxy(2);      
    if sum(sum(abs(u30-u3)))<1e-2
        %     k;
        break;
    else
        
        u30=u3;
    end
    
end
% 
 figure(5);plot3(X,Y,u3,'g.');hold on
XC(3)=c3x;
YC(3)=c3y;
%[sd(3),pos3]=simNorm(u3,[Z0,u1,u2]);
%[sd(3),pos3]=simNorm(u3,[u1,u2]);
%[sd(3),poss(3)]=simDist([c3x,c3y],[c1x,c1y;c2x,c2y]);
%belongs{3}=compBel(u3,Z-u3);
% belongs{3}=compBelDis2(u3,Z-u3,c3x,c3y,sig,crx,cry,sigr,L,X,Y,Satus);
% Satus(belongs{3})=1;
belongs{3}=compBelAmpF(u3,c3x,c3y,sig,X,Y);
[maxC(2),poss(2)]=simCount(belongs{3},belongs,2);

%% cycle 1 s4
%% update Z
q3=sum((Z-u3))/sum(Z);  %caculate the decreasing factor of alf, related to the power ratio
Z=Z-u3;

% SS(2)=sum(Z);
ZN1=Z./norm(Z);
% SD(1)=dot(ZN,ZN1); %caculate the similarity
SD(3)=norm(ZN1-ZN); %caculate the distance
%% 
% figure(7),plot(1,SD(1),'r.-','linewidth',2);hold on

% Cycle II decomposition
u40=zeros(L,1);
u40(p3)=1;
u4=u40;


% u1=u10;
m=10;
k=0;
%alf=q2*alf %update alf
beta=1;
        c4x=u40'*X/sum(u40);
        c4y=u40'*Y/sum(u40);
        u4r=Z-u40;
        crx=u4r'*X/sum(u4r);
        cry=u4r'*Y/sum(u4r);

        sig=[0.1,0;0,0.08];
        sigr=[0.1,0;0,0.08];
        % alf= 0.8*sqrt(sd1^2/(st1*sr1))   
        %0.5-1.5,
        alf=0.4;
         % alf is defined emperically using sr(distance from the center of remaining components)
         % sd(power) st(distance from the center of extracting component)
      error=rand(2,2)*1e-4;
      error=(error+error')/2;

while k<m

    for i=1:L

u4(i)=Z(i)*(1+beta*([X(i)-crx,Y(i)-cry]*pinv(sigr+error)*[X(i)-crx,Y(i)-cry]'))/(1+alf*([X(i)-c4x,Y(i)-c4y]*pinv(sig+error)*[X(i)-c4x,Y(i)-c4y]')+beta*([X(i)-crx,Y(i)-cry]*pinv(sigr+error)*[X(i)-crx,Y(i)-cry]'));
               
    end
    k=k+1;
          c4x=u4'*X/sum(u4);
          c4y=u4'*Y/sum(u4);
%          u1c=Z-u1;
%          crx=u1c'*X/sum(u1c);
%          cry=u1c'*Y/sum(u1c);
clear totpar sig sigr cxy crxy
         totpar=SigCal2(u4,Z,X,Y,c4x,c4y);
         sig=totpar{1};
         sigr=totpar{2};
         cxy=totpar{3};
         crxy=totpar{4};
         c4x=cxy(1);
         c4y=cxy(2);
         crx=crxy(1);
         cry=crxy(2);      
    if sum(sum(abs(u40-u4)))<1e-3
        %     k;
        break;
    else
        
        u40=u4;
    end
    
end
% 
 figure(5);plot3(X,Y,u4,'m.');hold off
 XC(4)=c4x;
YC(4)=c4y;
%[sd(4),pos4]=simNorm(u4,[Z0,u1,u2,u3]);
%[sd(4),pos4]=simNorm(u4,[u1,u2,u3]);
% [sd(4),poss(4)]=simDist([c4x,c4y],[c1x,c1y;c2x,c2y;c3x,c3y]);
%belongs{4}=compBel(u4,Z-u4);
%belongs{4}=compBelDis2(u4,Z-u4,c4x,c4y,sig,crx,cry,sigr,L,X,Y,Satus);
%Satus(belongs{4})=1;
belongs{4}=compBelAmpF(u4,c4x,c4y,sig,X,Y);
[maxC(3),poss(3)]=simCount(belongs{4},belongs,3);
%% cycle II s1 update Z
q11=sum(Z-u4)/sum(Z);
Z=Z-u4;


ZN4=Z./norm(Z);
% SD(4)=dot(ZN,ZN4);%caculate the similarity
SD(4)=norm(ZN4-ZN); %caculate the distance
% figure(7),plot(4, SD(4),'r.-','linewidth',2);hold on  %if the similarity flatten or increase(sample distance flatten or decrease), first loop extraction finished.
% figure(2);plot3(X,Y,Z,'r.');hold off
%% if SD increase or flatten(decrease or flatten is sd is sample distance), enter into loop two
%% cycle II s1
[val,peak]=max(Z);

u50=zeros(L,1);
u50(peak)=1;
u5=u50;


% u1=u10;
m=10;
k=0;

beta=1;
        c5x=u50'*X/sum(u50);
        c5y=u50'*Y/sum(u50);
        u5r=Z-u50;
        crx=u5r'*X/sum(u5r);
        cry=u5r'*Y/sum(u5r);

        sig=[0.2,0;0,0.08];
        sigr=[0.1,0;0,0.1];      
        %0.5-1.5,
        alf=0.2;       
      error=rand(2,2)*1e-3;
      error=(error+error')/2;

while k<m
        c5x=u50'*X/sum(u50);
        c5y=u50'*Y/sum(u50);
        urc=Z-u50;
        crx=urc'*X/sum(urc);
        cry=urc'*Y/sum(urc);
        for i=1:L        
u5(i)=Z(i)*(1+beta*([X(i)-crx,Y(i)-cry]*pinv(sigr+error)*[X(i)-crx,Y(i)-cry]'))/(1+alf*([X(i)-c5x,Y(i)-c5y]*pinv(sig+error)*[X(i)-c5x,Y(i)-c5y]')+beta*([X(i)-crx,Y(i)-cry]*pinv(sigr+error)*[X(i)-crx,Y(i)-cry]')); 
        end
    k=k+1;
    k;
    clear totpar sig sigr cxy crxy
         totpar=SigCal2(u5,Z,X,Y,c5x,c5y);
         sig=totpar{1};
         sigr=totpar{2} ;
         cxy=totpar{3};
         crxy=totpar{4};
         c5x=cxy(1);
         c5y=cxy(2);
         crx=crxy(1);
         cry=crxy(2); 

    if sum(sum(abs(u50-u5)))<1e-4
        %     k;
        break;
    else
        
        u50=u5;
    end
    
end
 figure(6);plot3(X,Y,u5,'b.');hold on; plot3(X,Y,Z,'k.');hold on


XC(5)=c5x;
YC(5)=c5y;
%[sd(5),pos5]=simNorm(u5,[Z0,u1,u2,u3,u4]);
%[sd(5),pos5]=simNorm(u5,[u1,u2,u3,u4]);
% [sd(5),poss(5)]=simDist([c5x,c5y],[c1x,c1y;c2x,c2y;c3x,c3y;c4x,c4y]);
% poss(5)
% sd(5)
%belongs{5}=compBel(u5,Z-u5);
% belongs{5}=compBelDis2(u5,Z-u5,c5x,c5y,sig,crx,cry,sigr,L,X,Y);
belongs{5}=compBelAmpF(u5,c5x,c5y,sig,X,Y);
[maxC(4),poss(4)]=simCount(belongs{5},belongs,4);
%%  update Z
q11=sum(Z-u5)/sum(Z);
Z=Z-u5;
ZN5=Z./norm(Z);
SD(5)=norm(ZN5-ZN);

%% cycle II s2
 [val,peak]=max(Z);
u60=zeros(L,1);
u60(peak)=1;
u6=u60;


% u1=u10;
m=10;
k=0;

%beta=1e-3;
beta=1;
        c6x=u60'*X/sum(u60);
        c6y=u60'*Y/sum(u60);
        u6r=Z-u60;
        crx=u6r'*X/sum(u6r);
        cry=u6r'*Y/sum(u6r);

        sig=[0.1,0;0,0.08];
        sigr=[0.1,0;0,0.08];      
        %0.5-1.5,
        alf=0.2;       
      error=rand(2,2)*1e-4;
      error=(error+error')/2;

while k<m
        c6x=u60'*X/sum(u60);
        c6y=u60'*Y/sum(u60);
        urc=Z-u60;
        crx=urc'*X/sum(urc);
        cry=urc'*Y/sum(urc);
        for i=1:L        
u6(i)=Z(i)*(1+beta*([X(i)-crx,Y(i)-cry]*pinv(sigr+error)*[X(i)-crx,Y(i)-cry]'))/(1+alf*([X(i)-c6x,Y(i)-c6y]*pinv(sig+error)*[X(i)-c6x,Y(i)-c6y]')+beta*([X(i)-crx,Y(i)-cry]*pinv(sigr+error)*[X(i)-crx,Y(i)-cry]')); 
        end
    k=k+1;
        clear totpar sig sigr cxy crxy
         totpar=SigCal2(u6,Z,X,Y,c6x,c6y);
         sig=totpar{1};
         sigr=totpar{2} ;
         cxy=totpar{3};
         crxy=totpar{4};
         c6x=cxy(1);
         c6y=cxy(2);
         crx=crxy(1);
         cry=crxy(2); 
    if sum(sum(abs(u60-u6)))<1e-4
        %     k;
        break;
    else
        
        u60=u6;
    end
    
end
 figure(6);plot3(X,Y,u6,'g.');hold on

XC(6)=c6x;
YC(6)=c6y;
%[sd(6),pos6]=simNorm(u6,[Z0,u1,u2,u3,u4,u5]);
%[sd(6),pos6]=simNorm(u6,[u1,u2,u3,u4,u5]);
% [sd(6),poss(6)]=simDist([c6x,c6y],[c1x,c1y;c2x,c2y;c3x,c3y;c4x,c4y;c5x,c5y]);
%belongs{6}=compBel(u6,Z-u6);
% belongs{6}=compBelDis2(u6,Z-u6,c6x,c6y,sig,crx,cry,sigr,L,X,Y);
belongs{6}=compBelAmpF(u6,c6x,c6y,sig,X,Y);
[maxC(5),poss(5)]=simCount(belongs{6},belongs,5);
%%
%%  update Z
q11=sum(Z-u6)/sum(Z);
Z=Z-u6;
ZN6=Z./norm(Z);
SD(6)=norm(ZN6-ZN);

%% cycle II s3

[val,peak]=max(Z);
u70=zeros(L,1);
u70(peak)=1;
u7=u70;


% u1=u10;
m=10;
k=0;

%beta=1e-3;
beta=1;
        c7x=u70'*X/sum(u70)
        c7y=u70'*Y/sum(u70)
        u7r=Z-u70;
        crx=u7r'*X/sum(u7r)
        cry=u7r'*Y/sum(u7r)

        sig=[0.1,0;0,0.08];
        sigr=[0.1,0;0,0.08];      
        %0.5-1.5,
        alf=0.2;       
      error=rand(2,2)*1e-4;
      error=(error+error')/2;

while k<m
        c7x=u70'*X/sum(u70);
        c7y=u70'*Y/sum(u70);
        urc=Z-u70;
        crx=urc'*X/sum(urc);
        cry=urc'*Y/sum(urc);
        for i=1:L        
u7(i)=Z(i)*(1+beta*([X(i)-crx,Y(i)-cry]*pinv(sigr+error)*[X(i)-crx,Y(i)-cry]'))/(1+alf*([X(i)-c7x,Y(i)-c7y]*pinv(sig+error)*[X(i)-c7x,Y(i)-c7y]')+beta*([X(i)-crx,Y(i)-cry]*pinv(sigr+error)*[X(i)-crx,Y(i)-cry]')); 
        end
    k=k+1;
        clear totpar sig sigr cxy crxy
         totpar=SigCal2(u7,Z,X,Y,c7x,c7y);
         sig=totpar{1};
         sigr=totpar{2} ;
         cxy=totpar{3};
         crxy=totpar{4};
         c7x=cxy(1);
         c7y=cxy(2);
         crx=crxy(1);
         cry=crxy(2); 
    if sum(sum(abs(u70-u7)))<1e-4
        %     k;
        break;
    else
        
        u70=u7;
    end
    
end
 figure(6);plot3(X,Y,u7,'r.');hold on

XC(7)=c7x;
YC(7)=c7y;
%[sd(7),pos7]=simNorm(u7,[Z0,u1,u2,u3,u4,u5,u6]);
%[sd(7),pos7]=simNorm(u7,[u1,u2,u3,u4,u5,u6]);
% [sd(7),poss(7)]=simDist([c7x,c7y],[c1x,c1y;c2x,c2y;c3x,c3y;c4x,c4y;c5x,c5y;c6x,c6y]);
%belongs{7}=compBel(u7,Z-u7);
belongs{7}=compBelAmpF(u7,c7x,c7y,sig,X,Y);
% belongs{7}=compBelDis2(u7,Z-u7,c7x,c7y,sig,crx,cry,sigr,L,X,Y);
[maxC(6),poss(6)]=simCount(belongs{7},belongs,6);
%%  update Z
q11=sum(Z-u7)/sum(Z);
Z=Z-u7;
ZN7=Z./norm(Z);
SD(7)=norm(ZN7-ZN);

%% cycle II s4
[val,peak]=max(Z);
u80=zeros(L,1);
u80(peak)=1;
u8=u80;


% u1=u10;
m=10;
k=0;

%beta=1e-3;
beta=1;

        c8x=u80'*X/sum(u80);
        c8y=u80'*Y/sum(u80);
        u8r=Z-u80;
        crx=u8r'*X/sum(u8r);
        cry=u8r'*Y/sum(u8r);

        sig=[0.1,0;0,0.08];
        sigr=[0.1,0;0,0.08];      
        %0.5-1.5,
        alf=1;       
      error=rand(2,2)*1e-4;
      error=(error+error')/2;

while k<m
        c8x=u80'*X/sum(u80);
        c8y=u80'*Y/sum(u80);
        urc=Z-u80;
        crx=urc'*X/sum(urc);
        cry=urc'*Y/sum(urc);
        for i=1:L        
u8(i)=Z(i)*(1+beta*([X(i)-crx,Y(i)-cry]*pinv(sigr+error)*[X(i)-crx,Y(i)-cry]'))/(1+alf*([X(i)-c8x,Y(i)-c8y]*pinv(sig+error)*[X(i)-c8x,Y(i)-c8y]')+beta*([X(i)-crx,Y(i)-cry]*pinv(sigr+error)*[X(i)-crx,Y(i)-cry]')); 
        end
    k=k+1;
      totpar=SigCal2(u8,Z,X,Y,c8x,c8y);
         sig=totpar{1};
         sigr=totpar{2} ;
         cxy=totpar{3};
         crxy=totpar{4};
         c8x=cxy(1);
         c8y=cxy(2);
         crx=crxy(1);
         cry=crxy(2); 
    if sum(sum(abs(u80-u8)))<1e-4
        %     k;
        break;
    else
        
        u80=u8;
    end
    
end
% figure(6);plot3(X,Y,u8,'m.',X,Y,Z-u8,'k.');hold off
figure(6);plot3(X,Y,u8,'m.');hold off
%figure(7);plot3(X,Y,Z-u8,'m.');hold off
XC(8)=c8x;
YC(8)=c8y;
%[sd(8),pos8]=simNorm(u8,[Z0,u1,u2,u3,u4,u5,u6,u7]);
%[sd(8),pos8]=simNorm(u8,[u1,u2,u3,u4,u5,u6,u7]);
%[sd(8),poss(8)]=simDist([c8x,c8y],[c1x,c1y;c2x,c2y;c3x,c3y;c4x,c4y;c5x,c5y;c6x,c6y;c7x,c7y]);
%belongs{8}=compBel(u8,Z-u8);
% belongs{8}=compBelDis2(u8,Z-u8,c8x,c8y,sig,crx,cry,sigr,L,X,Y);
belongs{8}=compBelAmpF(u8,c8x,c8y,sig,X,Y);
[maxC(7),poss(7)]=simCount(belongs{8},belongs,7);
%figure(7),plot(sd(2:end),'r.-');hold off
figure(7),plot(maxC(1:end),'r.-');hold off
utot=[u1,u2,u3,u4,u5,u6,u7,u8];

%% post merging
%from the figure 7 , we can determine the cluster # is 4

for j=4:7
    utot(:,poss(j))=utot(:,poss(j))+utot(j+1);
end
figure(8);
clors=["r.","b.","g.","m."];
for j=1:4
    plot3(X,Y,utot(:,j),clors(j));hold on
end
hold off

%% correctness 闇�瑕侀噸鏂板啓
%uf=[u1,u2,u3,u4];

belongss=classTot(X,Y,utot(:,1:4));
figure(1);
correctnessS=zeros(1,4);
datas{1}=data2;
[rr,cc]=size(data2);
lenths(1)=rr;
datas{2}=data4;
[rr,cc]=size(data4);
lenths(2)=rr;
datas{3}=data1;
[rr,cc]=size(data1);
lenths(3)=rr;
datas{4}=data3; %coincide with the previous extraction order
[rr,cc]=size(data3);
lenths(4)=rr;


for i=1:L
    p= belongss(i);
        clear datatemp
        datatemp=datas{p};
    if p==1

        plot(X(i),Y(i),'ro');hold on;
        for j=1:lenths(p)
            if X(i)==datatemp(j,1) && Y(i)==datatemp(j,2)
             correctnessS(p)= correctnessS(p)+1;
             break;
            end
        end
    elseif p==2
        plot(X(i),Y(i),'ko');hold on;
         for j=1:lenths(p)
            if X(i)==datatemp(j,1) && Y(i)==datatemp(j,2)
             correctnessS(p)= correctnessS(p)+1;
             break;
            end
        end
    elseif p==3
        plot(X(i),Y(i),'bo');hold on;
        for j=1:lenths(p)
            if X(i)==datatemp(j,1) && Y(i)==datatemp(j,2)
             correctnessS(p)= correctnessS(p)+1;
             break;
            end
        end
    else
        plot(X(i),Y(i),'go');hold on;
        for j=1:lenths(p)
            if X(i)==datatemp(j,1) && Y(i)==datatemp(j,2)
             correctnessS(p)= correctnessS(p)+1;
             break;
            end
        end
    end
end
hold off
correctnessS=correctnessS./lenths

figure(7),plot(maxC(1:end),'r*-','linewidth',2);xlabel('Times','FontSize',15);ylabel('Normalized Similarity','FontSize',15),grid on,hold off,hold off
% figure(5),plot(1:8,SD(1:8),'r*-','linewidth',2); %after one loop of full extraction, the similarity would increase or flatten instead of decrease(the sample distance decrease or flatten), so the number of clusters is determinated
% xlabel('Times','FontSize',15);ylabel('Normalized Similarity','FontSize',15),grid on,hold off




% for i=1:L
%     
%     [vs,p]=max(uf(i,:)); %may use distance from each center of each component(use statistical distance may be better?) much better
%     clust(i)=p;
%     if p==1
%         plot(X(i),Y(i),'bo');hold on;
%         for j=1:LLs(p)
%             if X(i)==data1(j,1) && Y(i)==data1(j,2)
%              correctnessS(p)= correctnessS(p)+1;
%              break;
%             end
%         end
%     elseif p==2
%         plot(X(i),Y(i),'ro');hold on;
%          for j=1:LLs(p)
%             if X(i)==data2(j,1) && Y(i)==data2(j,2)
%              correctnessS(p)= correctnessS(p)+1;
%              break;
%             end
%         end
%     elseif p==3
%         plot(X(i),Y(i),'go');hold on;
%         for j=1:LLs(p)
%             if X(i)==data3(j,1) && Y(i)==data3(j,2)
%              correctnessS(p)= correctnessS(p)+1;
%              break;
%             end
%         end
%     else
%         plot(X(i),Y(i),'ko');hold on;
%         for j=1:LLs(p)
%             if X(i)==data4(j,1) && Y(i)==data4(j,2)
%              correctnessS(p)= correctnessS(p)+1;
%              break;
%             end
%         end
%     end
% end
% hold off
