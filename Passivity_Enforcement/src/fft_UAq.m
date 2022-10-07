function [Trasformata]=fft_UAq(t,x)
% tale funzione trasforma secondo Fourier il segnale aperiodico
Timetot=t(size(t,2))-t(1); % lungh. finestra temporale
L=size(t,2);    % n. di camp. nel tempo (2^15)
Ts=t(2)-t(1);      % intervallo tra due camp. consecutivi
%dtem=Timetot/lenght;    % nel dominio del tempo

Fs=1/Ts;% freq di campionamento
fremax=Fs/2;    % freq di Nyquist
frefond=1/Timetot;    % Interv. fra due campioni consec. nel
windowtime=.5e-9;                     % dominio della freq.
windowL=windowtime*Fs;
t=(0:L-1)*Ts;

% vett è una funzione APERIODICA

% faccio la FFT
% della funzione APERIODICA

%XAA=spectrogram(x,hamming(windowL));
%[XAA,f]=pwelch(x,hamming(windowL),[],[],Fs);
win=hamming(windowL);
WinNum=L/windowL;

% reshape(x,
XAA=fft(x.*win')*Ts;
%XA=XAA(1:lenght/2+1);
% Per mettere d'accordo il risultato col metodo dei fasori è necessario
% raddoppiare lo spettro dalla seconda armonica
%XA(2:end)=2*XA(2:end);
% f=(0:lenght/2)*frefond; 
Trasformata=[f,XAA];
