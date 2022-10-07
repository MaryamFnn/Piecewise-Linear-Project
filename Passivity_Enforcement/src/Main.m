clear all;
close all;
format long e;
%%
%load data
cdw='Data/ex1_transmissionlineFDTD';
load(fullfile(cdw,'S_param_TD_step'));
load(fullfile(cdw,'S_param_FD.mat'));
A=size(S_param_TD);
port_number =A(1); 


figure
plot(time*1e9,abs(squeeze(S_param_TD(1,1,:))),time*1e9,abs(squeeze((S_param_TD(2,1,:)))),time*1e9,abs(squeeze((S_param_TD(3,1,:)))),time*1e9,abs(squeeze((S_param_TD(4,1,:)))),'linewidth',1.5)
legend('Input signal','|S11-TD|','|S21-TD|','|S31-TD|','|S41-TD|')
xlabel('Time(ns)','FontSize',12)
ylabel('SParameters step response','FontSize',12)
set(gca,'FontSize',10,'FontWeight', 'bold')

[logic,idxnonpassive_Data]=ispassive(s_param_FD);

%%
%PWL Model
s_param_store_t=zeros(port_number,port_number,length(frequency));
PWL_store_t=struct('idx_all_store',{},'x_node_store',{},'y_node_store',{},'TDError_store',{},'Error_store',{});
error_store=[];
for i=1:port_number
    for j=1:port_number        
        x = time;
        y = squeeze(S_param_TD(j,i,:))';
        error_store=[];
        x_node=[];
        y_node=[];
        idx_all=[];

        [~,idx] = min(x);
        idx_all(1) = idx;
        [~,idx] = max(x);
        idx_all(2)= idx;


        x_node = x (idx_all);
        y_node = y (idx_all);

        yplw = interp1(x_node, y_node, x);
        
        ABS_error = abs(y - yplw);
        error=sqrt(mean((y-yplw).^2));
        counter = 1;
        

        while error>.001

            [max_error,idxx] = max(ABS_error);
            idx_all = [idx_all idxx];
            idx_all=sort(idx_all);

            counter = counter + 1;
            

            x_node = x (idx_all);
            y_node = y (idx_all);

            yplw = interp1(x_node, y_node, x);
            ABS_error = abs(y - yplw);

            error=sqrt(mean((y-yplw).^2));
            error_store(counter)= error;
        end

        x_br = x_node;   %using on the curve data
        y_br = y_node;
        
        

        %x_br = XI*1e-9    %using off the curve data
        %y_br = YI         %using off the curve data

        num_br = length(x_br);
        for p=1:num_br-1
            deltaT(p) = x_br(p+1)-x_br(p);
            Amp(p)= (y_br(p+1)-y_br(p))/deltaT(p);   
            D(p) = x_br(p)+deltaT(p)/2;    
        end

        omega= 2*pi*frequency;
        S_param_sinc_model = Evaluation_Sinc_Model_Arbitrary_Frequencies(num_br,Amp,deltaT,D,omega); 
        s_param_store_t(j,i,:)=S_param_sinc_model;
        PWL_store_t(end+1)=struct('idx_all_store',{idx_all},'x_node_store',{x_br},'y_node_store',{y_br},'TDError_store',{error},'Error_store',{error_store});
    end
end

ypwl=interp1(PWL_store_t(1).x_node_store,PWL_store_t(1).y_node_store,x);

figure
plot(time*1e9,abs(squeeze(S_param_TD(1,1,:))),PWL_store_t(1).x_node_store*1e9,PWL_store_t(1).y_node_store,'ko','linewidth',1.5)
hold on 
plot(time*1e9,ypwl,'r--','linewidth',1.5)
hold off
legend('Real Data','PWL nodes','PWL Fitting')
xlabel('time(ns)','FontSize',12)
ylabel('S11 Step response','FontSize',12)
set(gca,'FontSize',10,'FontWeight', 'bold')
xlim([0,200])
figure
Abserror=abs(abs(squeeze(S_param_TD(1,1,:)))'-ypwl);
plot(time*1e9,Abserror)
legend('|Abs.Error|')
xlabel('time(ns)','FontSize',12)
ylabel('Abs.Error','FontSize',12)
set(gca,'FontSize',10,'FontWeight', 'bold')
xlim([0,200])


Singular_value_store_pwl_Sinc=[];
for g=1:length(frequency)
    
    S=s_param_store_t(:,:,g);
    [~,v,~]=svd(S);
    
    Singular_value_store_pwl_Sinc =[Singular_value_store_pwl_Sinc max(max(v))];        
end


PWL_Store_model = PWL_store_t;
s_param_store_pwl_sinc = s_param_store_t;
[logic,idxnonpassive_pwl_model]=ispassive(s_param_store_pwl_sinc);

figure
plot(frequency*1e-9,Singular_value_store_pwl_Sinc,'b','linewidth',1.5)
xlabel('frequency(GHz)','FontSize',12)
ylabel('Max Singular Value modeled data','FontSize',12)
title('Passivity check','FontSize',12)
yline(1,'r--','Threshold');
set(gca,'FontSize',10,'FontWeight', 'bold')


%%
%%%Aproach1

PWL_store_t = PWL_Store_model;
s_param_store_t = s_param_store_pwl_sinc;


passivity_counter = [length(idxnonpassive_pwl_model)];
iterrr =0;
iter_store=[iterrr];
error_store=[];
idxnonpassive_update1 = idxnonpassive_pwl_model;
tic
while logic~=1
        iterrr=iterrr+1;  
        iter_store=[iter_store iterrr];
        Frange=[min(idxnonpassive_update1) max(idxnonpassive_update1)];
        max_error_F=0;
        %see which s parameters has the most error in that ferequency range
        for i=1:4
            for k= 1:4
                Cal_Error = max(db(s_param_store_t(k,i,(Frange(1):Frange(2))))- db(abs(s_param_FD(k,i,(Frange(1):Frange(2))))));
                if Cal_Error > max_error_F
                    idxmax1=k;
                    idxmax2=i;
                    max_error_F =Cal_Error;
                end
                               
            end
        end
        
        x = time ;
        y = squeeze(S_param_TD(idxmax1,idxmax2,:))';
        g=(idxmax2-1)*4+idxmax1;
        x_br_t = PWL_store_t(g).x_node_store;
        y_br_t = PWL_store_t(g).y_node_store;
        yplw = interp1(x_br_t, y_br_t, x);
        ABS_error = abs(y - yplw);
        [max_error,idxx] = max(ABS_error);        
        idx_all_t=PWL_store_t(g).idx_all_store;       
        idx_all_t = [idx_all_t idxx];
        idx_all_t=sort(idx_all_t);      
        x_br = x (idx_all_t);
        y_br = y (idx_all_t);
        yplw = interp1(x_br_t, y_br_t, x);
        error=sqrt(mean((y-yplw).^2));

        num_br = length(x_br);
        for p=1:num_br-1
             deltaT(p) = x_br(p+1)-x_br(p);
             Amp(p)= (y_br(p+1)-y_br(p))/deltaT(p);   
             D(p) = x_br(p)+deltaT(p)/2;    
        end

        omega= 2*pi*frequency;
        S_param_sinc_model = Evaluation_Sinc_Model_Arbitrary_Frequencies(num_br,Amp,deltaT,D,omega); 
        s_param_store_t(idxmax1,idxmax2,:)=S_param_sinc_model;
        
        s_param_store_t(idxmax2,idxmax1,:)=S_param_sinc_model;
        
        PWL_store_t(g).idx_all_store=idx_all_t;
        PWL_store_t(g).x_node_store=x_br;
        PWL_store_t(g).y_node_store=y_br;
        PWL_store_t(g).TDError_store=error;
        
        z=(idxmax1-1)*4+idxmax2;
        PWL_store_t(z).x_node_store=x_br;
        PWL_store_t(z).y_node_store=y_br;
        PWL_store_t(z).idx_all_store=idx_all_t;
        PWL_store_t(z).TDError_store=error;
        
        [logic,idxnonpassive_update1]=ispassive(s_param_store_t);
        passivity_counter = [passivity_counter length(idxnonpassive_update1)];
        error_store=[error_store PWL_store_t(2).TDError_store];
                
end         

PWL_Store_update1 = PWL_store_t;
s_param_store_update1 = s_param_store_t;
[logic,idxnonpassive_update1]=ispassive(s_param_store_update1);

%%
%%%Aproach2
PWL_store_t = PWL_Store_model;
s_param_store_t = s_param_store_pwl_sinc;




[logic,idxnonpassive_update2]=ispassive(s_param_store_t);
passivity_counter2 = [length(idxnonpassive_update2)];
iter2 =0;
iter_store_2=[iter2];
FMaxError_store=[];
FMaxError_store2=[];

tic
while logic~=1
        iter2=iter2+1;  
        iter_store_2=[iter_store_2 iter2];
        Frange=[min(idxnonpassive_update2) max(idxnonpassive_update2)];
        max_error_F=0;
        %see which s parameters has the most error in that ferequency range
        for i=1:4
            for k= 1:4
                Cal_Error = max(db(s_param_store_t(k,i,(Frange(1):Frange(2))))- db(abs(s_param_FD(k,i,(Frange(1):Frange(2))))));
                if Cal_Error > max_error_F
                    idxmax1=k;
                    idxmax2=i;
                    max_error_F =Cal_Error;
                end
                               
            end
        end
      
        x = time ;
        y = squeeze(S_param_TD(idxmax1,idxmax2,:))';
        g=(idxmax2-1)*4+idxmax1;
        
        x_br = PWL_store_t(g).x_node_store;
        y_br = PWL_store_t(g).y_node_store;
        idx_all=PWL_store_t(g).idx_all_store; 
        
        x_temp=x;
        y_temp=y;
        
        x_temp(idx_all)=[];
        y_temp(idx_all)=[];
        
        for q=1:length(x_temp)
            
            x_br_t = [x_br x_temp(q)];
            y_br_t = [y_br y_temp(q)];
            
           
            
            [x_br_t,idxxx]=sort(x_br_t);
            y_br_t=y_br_t(idxxx);
            
            num_br = length(x_br_t);
            for p=1:num_br-1
                deltaT(p) = x_br_t(p+1)-x_br_t(p);
                Amp(p)= (y_br_t(p+1)-y_br_t(p))/deltaT(p);   
                D(p) = x_br_t(p)+(deltaT(p)/2);    
            end
           
            frequency_t=frequency(Frange(1):Frange(2));
            omega= 2*pi*frequency_t;
            S_param_sinc_model = Evaluation_Sinc_Model_Arbitrary_Frequencies(num_br,Amp,deltaT,D,omega);                     
            FMaxError=max(db(abs(S_param_sinc_model))- (db(abs(squeeze(s_param_FD(idxmax1,idxmax2,Frange(1):Frange(2))))))');
            FMaxError_store=[FMaxError_store FMaxError];
        end
        [MinError,idx]=min(FMaxError_store);
        idx5=find(x==x_temp(idx));
        idx_all=[idx_all idx5];
        idx_all= sort(idx_all);
        
        PWL_store_t(g).x_node_store=x(idx_all);
        PWL_store_t(g).y_node_store=y(idx_all);
        PWL_store_t(g).idx_all_store=idx_all;
        z=(idxmax1-1)*4+idxmax2;
        PWL_store_t(z).x_node_store=x(idx_all);
        PWL_store_t(z).y_node_store=y(idx_all);
        PWL_store_t(z).idx_all_store=idx_all;
        
        x_br = PWL_store_t(g).x_node_store;
        y_br = PWL_store_t(g).y_node_store;
        FMaxError_store=[];
        num_br = length(x_br);
        for p=1:num_br-1
             deltaT(p) = x_br(p+1)-x_br(p);
             Amp(p)= (y_br(p+1)-y_br(p))/deltaT(p);   
             D(p) = x_br(p)+deltaT(p)/2;    
        end

        omega= 2*pi*frequency;
        S_param_sinc_model = Evaluation_Sinc_Model_Arbitrary_Frequencies(num_br,Amp,deltaT,D,omega); 
        s_param_store_t(idxmax1,idxmax2,:)=S_param_sinc_model; 
        s_param_store_t(idxmax2,idxmax1,:)=S_param_sinc_model; 
        
        FMaxError2=max(db(abs(S_param_sinc_model))- (db(abs(squeeze(s_param_FD(idxmax1,idxmax2)))))');
        FMaxError_store2=[FMaxError_store2 FMaxError2];
        [logic,idxnonpassive_update2]=ispassive(s_param_store_t);
        passivity_counter2 = [passivity_counter2 length(idxnonpassive_update2)];
        %error_store=[error_store PWL_store(2).TDError_store];              
end   

PWL_Store_update2 = PWL_store_t;
s_param_store_update2 = s_param_store_t;
[logic,idxnonpassive_update2]=ispassive(s_param_store_update2);

figure 
plot(frequency,db(abs(squeeze(s_param_FD(2,1,:)))),frequency,db(abs(squeeze(s_param_store_pwl_sinc(2,1,:)))),'k',frequency,db(abs(squeeze(s_param_store_update1(2,1,:)))),'r--','linewidth',1.5)
legend('Real data','Sinc Representation','Approach 1')

figure 
plot(frequency,db(abs(squeeze(s_param_FD(2,1,:)))),frequency,db(abs(squeeze(s_param_store_pwl_sinc(2,1,:)))),'k',frequency,db(abs(squeeze(s_param_store_update2(2,1,:)))),'r--','linewidth',1.5)
legend('Real data','Sinc Representation','Approach 2')

