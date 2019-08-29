%Acce_Fast_Kurtogram - Creates a Kurtogram with the central frequency
%and the bandwidth
%
% Syntax:  Acce_Fast_Kurtogram
%
% Outputs:
%    KurtogramComparison.fig
%    KurtogramComparison.png
%
% MAT-files required: 
%   HFAccel_Dry_B7
%
% Author: Jorge Mijares
% email: jorge.mijares@tamu.edu
% Aug 2019; Last revision: 27-Aug-2019

%Program to Analize with Kurtogram accelerometer data
clc
clear
close all
%% Add Data and Functions folder
addpath(genpath('../data'))
addpath(genpath('../Functions'))

%% Fast Kurtogram
               % sampling frequency

nlevel = 7;     % number of decomposition levels
% Pre-whitening of the signal by AR filtering (optional)
prewhiten = 1;  % (always helpful in detection problems)

% load Accel_Dry_E3; i=13;
% load Accel_Lub_E3; i=15;
Bearing='B7';
Tests={['Dry_' Bearing]};
for k=1:length(Tests)
    

    filename=['HFAccel_' Tests{k}];
    load(filename); 
    vib=vibR;
    Ns=length(vib);
    for i=1:Ns 
        filename=[Tests{k} '_' Tag{i}];
        speed=mean(rpm_raw{i})/60;

        x=vib{:,i};

        if prewhiten == 1
           x = x - mean(x);
           Na = 100;
           a = lpc(x,Na);
           x = fftfilt(a,x);
           x = x(Na+1:end);		% it is very important to remove the transient of the whitening filter, otherwise the SK will detect it!!
        end

        path=[pwd '\' Bearing '_' Tests{k}];
        mkdir(path)
        str=[path '\KG_' Tag{i} ];
        Fc=0;  %Use zero for variable Envelope, >0 for fix Fc and lc
        lv=4;
        [cL{i},levL(i)]=Fast_kurtogram4(x,nlevel,Fs,str,Fc,lv);

        close 
        print([path '\Test_' num2str(i)],'-dpng','-r600')   
        saveas(gcf,[path '\Test_' num2str(i)])
    end
end
close all

%% Edit Figure 10
close all
addpath(genpath('./'))
fig=openfig('\Test_10');

line([17812 17812],ylim,'Color','r','LineStyle','--')
text(17812,2*2,['\leftarrow Fc = 17,812Hz' ],'Color','w')
text(18700,4*2,'\leftarrow Bw = 1,875Hz','Color','w')
text(15800,4*2,'\rightarrow','Color','w')
% text(15000,4*2-1,['Bw = 1,875Hz' ],'Color','w')
title('')
colorbar off

imageName = '.\KurtogramResults';
print(imageName,'-depsc','-r1000')
print(imageName,'-dtiffn','-r1000')
print(imageName,'-dpng','-r1000')%
saveas(gcf,imageName)