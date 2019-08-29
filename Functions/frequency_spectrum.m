
function [db,normal,f]=frequency_spectrum(xp,Fs)
Parseval=0;
if Parseval==1 %preserve energy
    dt=1/Fs;
    L=length(xp);
    df=Fs/L;
    xp=xp-mean(xp);
%     xp = xp(:).*hanning(L); 
    energy_a=sum(xp.*conj(xp)*dt);  %Check Energy for Parseval Theorem a=b;
    Y = fft(xp,L)*dt;
    energy_b=sum(Y.*conj(Y)*df);
    f = Fs/2*linspace(0,1,round(L/2)+1);
    
    db=20*log10(2*abs(Y(1:round(L/2)+1)));
    normal=2*abs(Y(1:round(L/2)+1));
else %Preserve Amplitud
    L=length(xp);
    xp=xp-mean(xp);   %substract average
    xp = xp(:).*hanning(L); 
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    Y = fft(xp,NFFT)/NFFT;
    f = Fs/2*linspace(0,1,NFFT/2+1);
    % Plot single-sided amplitude spectrum.
    db=20*log10(2*abs(Y(1:NFFT/2+1)));
    normal=2*abs(Y(1:NFFT/2+1));
end



end