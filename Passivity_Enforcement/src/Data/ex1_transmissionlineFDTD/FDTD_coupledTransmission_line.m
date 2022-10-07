clear
close all
clc
format long e
%%
%%% DUT is a system of two coupled transmission lines: four-port system. 
%%%
%%%          P1  --> |         | --> P3
%%%                  |   DUT   |
%%%          P2  --> |         | --> P4
%%%
%%% The telegrapher's equations are solved directly in the time-domain via
%%% the finite-difference time-domain (FDTD) method. 
%%% Two factors simplify the problem and speed-up the solution of the FDTD:
%%%  1) the per-unit-length (p.u.l.) parameters are frequency-independent
%%%  2) the geometry of the structure allows one to avoid the full
%%%     3-dimensional FDTD formulation 
%%% Note that the approximation 1) limits the ability to accurately 
%%% characterize the DUT at high frequencies


%%
%%% Definition of time-domain source 


%%%Step Source (ADS Data)
load ('StepSignalrise1nts100p.txt')
time = StepSignalrise1nts100p(:,1)';
vs = StepSignalrise1nts100p(:,2)';
nt=length(time); %Number of samples
Deltat=time(2);  %Time step

%%
%%% Physical parameters 
mu0=4*pi*1e-7;
eps0=8.854e-12;
c0=1/sqrt(mu0*eps0);

%%% Characteristics of the MTLs and p.u.l. definition
ncond=2; % Number of conductors
lung=1; % Length of the MTLs expressed in meters. This can be considered as a design parameter.
rW=7.5*2.54e-5;
sigma=5.8e7;

%%% Pul parameters

Rpul=[2/(sigma*pi*rW^2) 1/(sigma*pi*rW^2); 1/(sigma*pi*rW^2) 2/(sigma*pi*rW^2)];   % Values taken from a reference example
Lpul=[0.7485 0.5077; 0.5077 1.0154]*1e-6;
Cpul=[37.432 -18.716; -18.716 24.982]*1e-12;
Gpul=[1e-15 0; 0 1e-15];

%%
%Frequency Response Using FFT all Data

%%% Port resistances
Rs=[50 0;0 50];   % For easier FDTD formulation, the full termination matrix is divided into source and load terminations blocks.
Rl=[50 0;0 50];   % This is equivalent to terminating each port with a 50Ohm resistor (needed to compute the time-domain s-parameters)

%%% Computing time-domain solution via FDTD 

V=zeros(2*ncond,nt);



%vs_MTL(1,:)=vs;  %%% In this configuration, only port 1 is excited. The solution gives the first column of the time-domain s-parameters matrix.
                 %%% This line has to be changed to  
                 %%%     vs_MTL(2,:)=vs; to obtain the second column of the time-domain s-parameters       
                 %%%     vl_MTL(1,:)=vs; to obtain the third column of the time-domain s-parameters 
                 %%%     vl_MTL(2,:)=vs; to obtain the fourth column of the time-domain s-parameters 
                 %%% Size of overall time-domain step s-param matrix: 4x4xnt
                 
% Solving FDTD: nothing to change from here on
Deltaz=1.001*c0*Deltat;   % Condition to enforce passivity of the solution
nz=ceil(lung/Deltaz);
Deltaz=lung/nz;
%FDTD
% Lossy case

S_param_TD=zeros(4,4,nt);
for t=1:4
    V=zeros(2*ncond,nt);
    vs_MTL = zeros(2,nt); 
    vl_MTL= zeros(2,nt);
    V(t,:)=vs;
    vs_MTL = V(1:2,:); 
    vl_MTL=V(3:4,:);
    
    v=zeros(ncond,nz+1,nt);
    i=zeros(ncond,nz,nt);
    
    for n=1:length(time)-1
        v(:,1,n+1)=(Deltaz/Deltat*Rs*Cpul+Deltaz/2*Rs*Gpul+eye(ncond))^(-1)*((Deltaz/Deltat*Rs*Cpul-Deltaz/2*Rs*Gpul-eye(ncond))*squeeze(v(:,1,n))-2*Rs*squeeze(i(:,1,n))+(squeeze(vs_MTL(:,n+1))+squeeze(vs_MTL(:,n))));
        for k=2:nz
           v(:,k,n+1)=(Deltaz/Deltat*Cpul+Deltaz/2*Gpul)^(-1)*(Deltaz/Deltat*Cpul-Deltaz/2*Gpul)*squeeze(v(:,k,n))-(Deltaz/Deltat*Cpul+Deltaz/2*Gpul)^(-1)*(squeeze(i(:,k,n))-squeeze(i(:,k-1,n)));
        end
        v(:,nz+1,n+1)=(Deltaz/Deltat*Rl*Cpul+Deltaz/2*Rl*Gpul+eye(ncond))^(-1)*((Deltaz/Deltat*Rl*Cpul-Deltaz/2*Rl*Gpul-eye(ncond))*squeeze(v(:,nz+1,n))+2*Rl*squeeze(i(:,nz,n))+(squeeze(vl_MTL(:,n+1))+squeeze(vl_MTL(:,n))));
        for k=1:nz
            i(:,k,n+1)=(Deltaz/Deltat*Lpul+Deltaz/2*Rpul)^(-1)*(Deltaz/Deltat*Lpul-Deltaz/2*Rpul)*squeeze(i(:,k,n)) - (Deltaz/Deltat*Lpul+Deltaz/2*Rpul)^(-1)*(squeeze(v(:,k+1,n+1))-squeeze(v(:,k,n+1)));
        end
        v1_FDTD_lossy(:,n+1)=squeeze(v(:,1,n+1));
        v2_FDTD_lossy(:,n+1)=squeeze(v(:,nz+1,n+1));
        %i1_FDTD_lossy(:,n+1)=squeeze(i(:,1,n+1));
        %i2_FDTD_lossy(:,n+1)=squeeze(i(:,nz+1,n+1));
        %Z1=v1_FDTD_lossy(1,:)./i1_FDTD_lossy(1,:);
    end
    for j=1:2
        S_param_TD(j,t,:)=v1_FDTD_lossy(j,:);
    end
    for j=3:4
        S_param_TD(j,t,:)=v2_FDTD_lossy(j-2,:);
    end
        
   
end



%Rect_Sinc all data
frequency=linspace(0,2e9,2000);
%frequency=freq;
for i=1:4
    for j=1:4
        x_br = time;   %using on the curve data
        y_br = S_param_TD(j,i,:);
        num_br = length(time);
        for p=1:num_br-1
            deltaT(p) = x_br(p+1)-x_br(p);
            Amp(p)= (y_br(p+1)-y_br(p))/deltaT(p);   
            D(p) = x_br(p)+(deltaT(p)/2);    
        end
        omega= 2*pi*frequency;
        S_param_sinc_model = Evaluation_Sinc_Model_Arbitrary_Frequencies(num_br,Amp,deltaT,D,omega);
        s_param_FD(j,i,:)=S_param_sinc_model;
        
    end
end

save('./S_param_TD_step.mat','S_param_TD','time');
save('./S_param_FD.mat','s_param_FD','frequency');

%{
figure
%plot(frequency*1e-9,db(abs(squeeze(s_param_data_store_sinc_tot(1,1,:)))),'--',frequency*1e-9,db(abs((squeeze(S_param_FD(1,1,:))))),'linewidth',1.5)
plot(frequency*1e-9,db(abs(squeeze(s_param_data_store_sinc_tot(1,1,:)))),'linewidth',1.5)
hold on
plot(frequency*1e-9,db(abs(squeeze(s_param_data_store_sinc_tot(1,2,:)))),'linewidth',1.5)
legend('RectSinc all data')
xlabel('Frequency(GHz)','FontSize',12)
ylabel('S11','FontSize',12)
set(gca,'FontSize',10,'FontWeight', 'bold')





%check Passivity
Singular_value_store2=[];
for g=1:length(freq)
    
    S=s_param_data_store_FFT(:,:,g);
    [~,v,~]=svd(S);
    Singular_value_store2 =[Singular_value_store2 max(max(v))];
    
end

figure
plot(freq*1e-9,Singular_value_store2,'b','linewidth',1.5)
xlabel('frequency(GHz)','FontSize',12)
ylabel('Max Singular Value real data','FontSize',12)
title('Passivity check','FontSize',12)
set(gca,'FontSize',10,'FontWeight', 'bold')
yline(1,'r--','Threshold');
xlim([0,1.2])
%}
[logic,Idxnonpassive]=ispassive(s_param_FD);






