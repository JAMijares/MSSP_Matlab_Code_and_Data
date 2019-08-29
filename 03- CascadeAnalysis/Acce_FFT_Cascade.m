%Acce_FFT_Cascade - Creates the base figure that plots all the FFT in a
%Cascade, later AddAnnotations add all the labels and annotations to the
%generated figure
%
% Syntax:  Acce_FFT_Cascade
%
% Outputs:
%    Bearing [BearingID].fig
%    Bearing [BearingID].png
%
% Other m-files required: PlotCascade
% Subfunctions: PlotCascade
% MAT-files required: 
%    Bearings={'B7'};    %Bearing ID
%    str{1}=['HFAccel_Dry_' Bearings{k}];
%    str{2}=['HFAccel_Dry_P1_' Bearings{k} ];
%    str{3}=['HFAccel_Lub_min_' Bearings{k} ];
%    str{4}=['HFAccel_Lub_full_' Bearings{k} ];
%
% Author: Jorge Mijares
% email: jorge.mijares@tamu.edu
% Aug 2019; Last revision: 27-Aug-2019

clc
clear
close all
%% Add Data
addpath(genpath('../data'))
addpath(genpath('../functions'))

%% Select Bearings to test
% Bearings={'B71','B7','B8','B9','B10','B11'};
Bearings={'B7'};

Labels={'Dry' ;  'Dry + \newline Interference';...
            'Lub 5%'; 'Lub 100%'};

%% Load Data
for k=1:length(Bearings)

    str{1}=['HFAccel_Dry_' Bearings{k}];
    str{2}=['HFAccel_Dry_P1_' Bearings{k} ];
    str{3}=['HFAccel_Lub_min_' Bearings{k} ];
    str{4}=['HFAccel_Lub_full_' Bearings{k} ];
    
    final=87300;
    start=final/20;
    n0=1;
    
    for j=1:length(str)
        load(str{j})
        x=vibR;
        str1{n0}=' ';
        label=Labels{j};
        PlotCascade
        n0=n0+N+1;
        clear rpm_raw vibR vibL Fs Acc_Data t filename vib dbr mag fr Temp
    end
    hold off
    ylim([0 length(str1)])
    set(gca,'YTick',0:length(str1)-1)
    set(gca,'YTickLabel',str1)

    xlabel('Frequency (kHz)')
    ylabel('Conditions')
    zlabel('Magnitude')

    set(gcf,'Position', [300 250      833         517])
    az = 20;
    el = 70;
    view(az, el);
    imageName = ['Bearing ' Bearings{k}];
    print(imageName,'-depsc','-r1000')
    print(imageName,'-dtiffn','-r1000')
    print(imageName,'-dpng','-r1000')%
    saveas(gcf,imageName)

end


