%Main_IndicatorsComparison_4Conditions - Creates the figures for the
%comparison of RMS, Crest Factor and Kurtosis
%
% Syntax:  Main_IndicatorsComparison_4Conditions
%
% Outputs:
%    RMSComparison.fig
%    RMSComparison.png
%    CrestComparison.fig
%    CrestComparison.png
%    KurtComparison.fig
%    KurtComparison.png
%
% MAT-files required: 
%   Bearing_ID = 'B71','B7','B8','B9','B10','B11'
%   HFAccel_Dry_[Bearing_ID]
%   HFAccel_Dry_P1_[Bearing_ID]
%   HFAccel_Lub_min_[Bearing_ID];
%   HFAccel_Lub_full_[Bearing_ID];
%
% Author: Jorge Mijares
% email: jorge.mijares@tamu.edu
% Aug 2019; Last revision: 27-Aug-2019

clc
clear
close all
%% Add Data
addpath(genpath('../Data'))

%% Select Bearings to test
Bearings={'B71','B7','B8','B9','B10','B11'};

%% Select Markers
markers={'s','^','o','x','d','<','>'};

%% Load Data and plot mean and uncertainties 
for k=1:length(Bearings)

    str{1}=['HFAccel_Dry_' Bearings{k}];
    str{2}=['HFAccel_Dry_P1_' Bearings{k} ];
    str{3}=['HFAccel_Lub_min_' Bearings{k} ];
    str{4}=['HFAccel_Lub_full_' Bearings{k} ];

    for j=1:length(str)
        load(str{j})
        N=length(t);
        for i=1:N
            if (length(Bearings{k})>2 && (all(Bearings{k} == 'B26') || all(Bearings{k} == 'B27')))
                vib=vibR_Y{:,i};
            else
                vib = vibR{:,i};
            end
            vib=vib-mean(vib);
            RMS_v(i,j)=rms(vib);
            Kurt(i,j)=round(kurtosis(vib),2);
            CF(i,j)=max(abs(vib))/rms(vib);
        end
    end

figure(1)
alpha=0.05;
n=length(RMS_v);
nu=n-1;
u_rnd=abs(tinv(alpha,nu))*(std(RMS_v)/sqrt(n));
errorbar([1:2:7]+0.2*(k-1),mean(RMS_v),u_rnd,markers{k},'MarkerSize',10); hold on

figure(2)
alpha=0.05;
n=length(Kurt);
nu=n-1;
u_rnd=abs(tinv(alpha,nu))*(std(Kurt)/sqrt(n));
errorbar([1:2:9]+0.2*(k-1),[-1 mean(Kurt)],[0 u_rnd],markers{k},'MarkerSize',10); hold on

figure(3)
alpha=0.05;
n=length(CF);
nu=n-1;
u_rnd=abs(tinv(alpha,nu))*(std(CF)/sqrt(n));
errorbar([1:2:7]+0.2*(k-1),mean(CF),u_rnd,markers{k},'MarkerSize',10); hold on

end

%% Modify RMS plot
figure(1)
set(gcf,'Position',[-1036         268         729         553])


xlim([0.1*(k-1) 8+0.1*(k-1)])
limits=ylim;
% ylim([0 0.8])
set(gca,'xtick',[1:2:7]+0.1*(k-1))
Test_name={'Dry' ;  'Dry + \newline Interference';...
            'Lub 5%'; 'Lub 100%'};
set(gca,'xticklabel',Test_name)
ylabel('RMS (g)')
xlabel('Conditions')
for j=0:2:8
line([j j]+0.1*(k-1), ylim,'Color',[0.7 0.7 0.7])
end
legend('Experiment #1','Experiment #2','Experiment #3','Experiment #4',...
        'Experiment #5','Experiment #6')

imageName = 'RMSComparison';
print(imageName,'-depsc','-r1000')
print(imageName,'-dtiffn','-r1000')
print(imageName,'-dpng','-r1000')%
saveas(gcf,imageName)

%% Modify Kurt plot
figure(2)
set(gcf,'Position',[250         200         729         553])
xlim([0 10+0.1*(k-1)])
limits=ylim;
ylim([0 75])
set(gca,'xtick',[1:2:9]+0.1*(k-1))
Test_name={'Levels'; 'Dry' ;  'Dry + \newline Interference';...
            'Lub 5%'; 'Lub 100%'};
set(gca,'xticklabel',Test_name)
ylabel('Kurtosis')
xlabel('Conditions')
fontsize =12;
line(xlim,[3 3],'LineStyle','-','Color',[0.7 0.7 0.7])%[0 0.5 0])
text(0.1,1.6,'Good','Color',[0 0.5 0],'FontSize',fontsize)
line(xlim,[7 7],'LineStyle','-','Color',[0.7 0.7 0.7])%[0.85 0.85 0])
text(0.1,5.,'Satisfactory','Color',[0.8 0.8 0],'FontSize',fontsize)
line(xlim,[10 10],'LineStyle','-','Color',[0.7 0.7 0.7])%[0.9 0.4 0])
text(0.1,8.75,'Unsatisfactory','Color',[0.9 0.4 0],'FontSize',fontsize)
text(0.1,45,'Unacceptable','Color',[0.9 0 0],'FontSize',fontsize)

for j=2:2:10
line([j j]+0.1*(k-1), ylim,'Color',[0.7 0.7 0.7])
end

legend('Experiment #1','Experiment #2','Experiment #3','Experiment #4',...
        'Experiment #5','Experiment #6')

imageName = 'KurtComparison';
print(imageName,'-depsc','-r1000')
print(imageName,'-dtiffn','-r1000')
print(imageName,'-dpdf','-r1000','-fillpage')
print(imageName,'-dpng','-r1000')%
saveas(gcf,imageName)

%% Modify Kurt plot
figure(3)
set(gcf,'Position',[-1036         268         729         553])
xlim([0.1*(k-1) 8+0.1*(k-1)])
limits=ylim;
% ylim([0 0.8])
set(gca,'xtick',[1:2:7]+0.1*(k-1))
Test_name={'Dry' ;  'Dry + \newline Interference';...
            'Lub 5%'; 'Lub 100%'};
set(gca,'xticklabel',Test_name)
ylabel('CF')
xlabel('Conditions')

for j=0:2:8
line([j j]+0.1*(k-1), ylim,'Color',[0.7 0.7 0.7])
end

legend('Experiment #1','Experiment #2','Experiment #3','Experiment #4',...
        'Experiment #5','Experiment #6')
    
imageName = 'CrestComparison';
print(imageName,'-depsc','-r1000')
print(imageName,'-dtiffn','-r1000')
print(imageName,'-dpng','-r1000')%
saveas(gcf,imageName)

% set(gca,'ygrid','on')
