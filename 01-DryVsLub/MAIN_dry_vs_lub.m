%MAIN_dry_vs_lub - Creates the figure for the comparison of dry vs lub
%bearing
%
% Syntax:  MAIN_dry_vs_lub
%
% Outputs:
%    DryVsLub_[Experiment].fig
%    DryVsLub_[Experiment].png
%
% MAT-files required: 
%       Bearings={'B26','B27','B28','B31','B32','B34'};
%       'HFAccel_Dry_Bearings{:}'
%       'HFAccel_Lub_full_Bearings{:}'
%
% Author: Jorge Mijares
% email: jorge.mijares@tamu.edu
% Aug 2019; Last revision: 27-Aug-2019

clc
clear
close all
%% Add data path

addpath(genpath('../data'))
Bearings={'B26','B27','B28','B31','B32','B34'};

%% Select labels and markers
str2={'Y Axis','X Axis','Z Axis'};
markers={'s','^','o','x','d','<','>','v'};
Col={'r' 'b'};

%% Load Bearing data and plot Lub vs Dry
for k=1:1

    str{1}=['HFAccel_Dry_' Bearings{k}];
    str{2}=['HFAccel_Lub_full_' Bearings{k} ];
    
    n0=0;
    figure
    set(gcf,'Position',[300 250 800 200])

    for j=1:length(str)
        
        load(str{j})
        N=length(t);
        %Calculate the impulsivity of each signal
        for i=1:N 
            vib=vibR_Y{:,i};
            vib=vib-mean(vib);
            Kurt(i,j,1)=round(kurtosis(vib),2);

            vib=vibR_X{:,i};
            vib=vib-mean(vib);
            Kurt(i,j,2)=round(kurtosis(vib),2);

            vib=vibR_Z{:,i};
            vib=vib-mean(vib);
            Kurt(i,j,3)=round(kurtosis(vib),2);
        end

        if j==1 
            [mx,indx]=max(Kurt(:,j,2));
            vib=vibR_Y{:,indx};
            plot(t{indx},vib-mean(vib),'k');hold on       
        else


            [mx,indx]=max(Kurt(:,j,2));
            vib=vibR_Y{:,indx};
            plot(t{indx},vib-mean(vib),'Color',[0.6 0.6 0.6]);hold on
            ylabel('Acceleration (g)')
            [leg,hObj]=legend('Dry Bearing','Lubricated Bearing','Orientation','horizontal');
            leg.Position=[0.3922    0.8932    0.2482    0.0488];
            hL=findobj(hObj,'type','line');  % get the lines, not text
            set(hL,'linewidth',2)            % set their width property
        end
        xlabel('Time (s)')

    end
%     
%     imageName = ['DryVsLub_' num2str(k)];
%     print(imageName,'-depsc','-r1000')
%     print(imageName,'-dtiffn','-r1000')
%     print(imageName,'-dpdf','-r1000')
%     print(imageName,'-dpng','-r1000')%
%     print(imageName, '-djpeg','-r1000')   % save as COLOR pdf file
%     saveas(gcf,imageName)
end
imageName = 'DryVsLub' ;
print(imageName,'-depsc','-r1000')
print(imageName,'-dtiffn','-r1000')
print(imageName,'-dpng','-r1000')%
print(imageName, '-djpeg','-r1000')   % save as COLOR pdf file
saveas(gcf,imageName)
