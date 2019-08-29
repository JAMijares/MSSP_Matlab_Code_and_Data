
N=length(t); %Obtain number of tests
for i=1:N
%     t_lub(:,i)=t{i};
%     rpm_lub(:,i)=rpm_raw{:,i};
    vib = x{:,i};
    vib=vib-mean(vib);
    [dbr,mag,fr]= frequency_spectrum(vib,Fs);
%     waterfall(gpuArray(fr(start:final)/1000),i+n0,gpuArray(mag(start:final)')); hold on
    waterfall(fr(start:final)/1000,i+n0,mag(start:final)'); hold on
    if i==round(N/2)
            str1{i+n0}=label;
    else
            str1{i+n0}=' ';
    end    
end