%MAIN_Envelope_Analysis_FTF_BPFO - Creates the Envelope Analysis with the
%comparison of BPFO and FTF
%
% Syntax:  MAIN_Envelope_Analysis_FTF_BPFO
%
% Outputs:
%    FTF_BPFO\[Bearing_ID].fig
%    FTF_BPFO\[Bearing_ID].png
%
% Other m-files required: Fast_kurtogram_noPlots
% Subfunctions: Fast_kurtogram_noPlots
% MAT-files required: 
%       HFAccel_Lub_Full_[Bearing_ID]
%       HFAccel_Dry_[Bearing_ID]
%
% Author: Jorge Mijares
% email: jorge.mijares@tamu.edu
% Aug 2019; Last revision: 27-Aug-2019

clc
clear
close all
%% Add Data and Functions folder
addpath(genpath('../data'))
addpath(genpath('../Functions'))

LoadData=false;  %%Enable to load previous processed data to save time.
%% Fast Kurtogram

nlevel = 7;     % number of decomposition levels
% Pre-whitening of the signal by AR filtering (optional)
prewhiten = 1;  % (always helpful in detection problems)


set1={'B71','B7','B8','B9','B10','B11'};
Fc_DB1={17770, 17770 , 18750, 17770, 17770, 17770}; 
set2={'B26','B27','B28','B31','B32','B34'};
Fc_DB2={17770, 17770 , 8500, 17770, 17770, 17770};

Fc_DB= {Fc_DB1{:}, Fc_DB2{:}};
Bearing = {set1{:}, set2{:}};

for k=1:1
    Bearing{k}
    Tests={['Lub_Full_' Bearing{k}];
            ['Dry_' Bearing{k}]
            };
    if(~LoadData)
        for j=1:length(Tests)

            filename=['HFAccel_' Tests{j}];
            load(filename);

            if (any(strcmp(Bearing{k},set2)))
                vibR=vibR_Y;
            end

            Ns=length(vibR);
            for i=1:Ns 
                filename=[Tests{j} '_' Tag{i}];
                speed=mean(rpm_raw{i})/60;

                x=vibR{:,i};
                if prewhiten == 1
                   x = x - mean(x);
                   Na = 100;
                   a = lpc(x,Na);
                   x = fftfilt(a,x);
                   x = x(Na+1:end);		% it is very important to remove the transient of the whitening filter, otherwise the SK will detect it!!
                end

                %Fix Central Frequency to filter, use 0 to choose the automatic
                Fc=Fc_DB{k};%17770;
                lv=3.5;

                %Apply Fast-Kurtogram to signal
                Folder=[Bearing{k} '_' Tests{j} '_Test' num2str(i)];
                [cL,levL]=Fast_kurtogram_noPlots(x,nlevel,Fs,Folder,Fc,lv);


                Sc=levL;
                level = fix(Sc) + (rem(Sc,1)>=0.5)*(log2(3)-1);
                nfft = 2*ceil(length(cL)/2);
                env = abs(cL).^2;
                S(:,i) = abs(fft((env(:)-mean(env)).*hanning(length(env)),nfft)/nfft);

            end
            S_mean{j,k}=mean(S,2);
            f=speed;
            bar=5;
            fr = linspace(0,.5*Fs/2^level,nfft/2);
            index = find(fr>=f); 
            [y2,I]=max(S_mean{j,k}(index(1)-bar:index(1)+bar));
            Fix_f=fr(I+index(1)-bar-1)/f;
            error_1X=100*(1-Fix_f);
            if(abs(error_1X) <= 2)
                f=f*Fix_f;
            end
            f_rot{j,k}=f;

        end

        save('SpectrumData','S_mean','speed','Fs','level','nfft','f_rot','fr')
    end
    %% FFT of Envelope
    
    load('SpectrumData')
    figure(1)
    set(gcf,'Position',[  250 200        800         500])
    
    %Find 1x and correct speed
    
    Bearing6204_10_Frequencies

    % Lub
    p0{1}= semilogy(fr/f_rot{1,k},S_mean{1,k}(1:nfft/2),'b','LineWidth',1), hold on
    % Dry
    p0{2}= semilogy(fr/f_rot{2,k},S_mean{2,k}(1:nfft/2),'r','LineWidth',1), hold off
    % Preload
%     p0{3}= semilogy(fr/f_rot{3,k},S_mean{3,k}(1:nfft/2),'b','LineWidth',1), hold off
    
    xlim([0 20])
%     yyaxis left
    ylabel('FFT Envelope (log_{10})')
    xlabel('Orders')
  

    %% FIND Frequencies
      
    p1=line([1 1],ylim,'Color','k','LineStyle','--')
    tag0='1X';

    Faults={'FTF' ,'BPFO', 'BPFI', 'BSF' ,'BDF'};
    Faults_val={FTF, BPFO, BPFI, BSF, 2*BSF};
    
    bar=10;
    index = find(fr>=Faults_val{1}*f_rot{2,k}); 

    offset=[0.004,-0.0038,-0.02450,-0.001,-0.036,0,...
            +0.002,-0.0225,-0.002,-0.004,+0.002,-0.002];
    
    [y2,I]=max(S_mean{2,k}(index(1):index(1)+bar));
    x=fr(I+index(1)-1)/f_rot{2,k}+offset(k);
    error=round(100*(Faults_val{1}-x)/Faults_val{1},1);
    Faults_val{1}=x;

    for h=1:100
        p2=line([h*x h*x],ylim,'Color',[77,190,238]/255,'LineStyle',':','LineWidth',1)
        tag1='Fundamental Train Frequency (FTF)';
    end
    
    for h=1:5
        p3=line([8*h*x 8*h*x],ylim,'Color',[119,200,48]/255,'LineStyle','-.','LineWidth',1)
        tag2='Ball Pass Frequency on the Outer Ring (BPFO)';
    end
    
    title(['Experiment #' num2str(k)])
    legend([p0{1}(1) p0{2}(1) p1(1) p2(1) p3(1)],'Lubricated Bearing','Dry Bearing',tag0,tag1,tag2)

    Folder='FTF_BPFO\';
    if(not(exist(Folder,'dir')))
        mkdir(Folder)
    end
%  

    imageName = [Folder  Bearing{k}];
%     print(imageName,'-depsc','-r1000')
%     print(imageName,'-dtiffn','-r1000')
%     print(imageName,'-dpng','-r1000')%
    saveas(gcf,imageName)
    

end

imageName = 'EnvelopeSpectrum';
print(imageName,'-depsc','-r1000')
print(imageName,'-dtiffn','-r1000')
print(imageName,'-dpng','-r1000')%
saveas(gcf,imageName)

