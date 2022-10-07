

function S_param_sinc = Evaluation_Sinc_Model_Arbitrary_Frequencies(Num_break,Amp_sinc,DeltaT,Di,omega) %,Size_Spar_data)



Aux=0;
S_par=0;
for g=1:(Num_break-1)
    Aux= Amp_sinc(g)*DeltaT(g)*(sin(omega*0.5*DeltaT(g))./(omega*0.5*DeltaT(g))).*exp(-1i*omega*Di(g));
    S_par=S_par+Aux;
end
S_param_sinc=S_par;


%%%% DC Estimation
clear Aux S_par
if omega(1)==0
    
    Aux=0;
    S_par=0;
    for g=1:(Num_break-1)
        Aux= Amp_sinc(g)*DeltaT(g);
        S_par=S_par+Aux;
    end
    
    S_param_sinc(1)=S_par;
    
end

    
% if Size_Spar_data == size(S_param_sinc,1)
% 
% else
%     
%     S_param_sinc=S_param_sinc.';
%     
% end

end
