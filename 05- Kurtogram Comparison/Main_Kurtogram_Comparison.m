%Main_Kurtogram_Comparison - Creates the figure for the comparison of
%Kurtogram
%
% Syntax:  Main_Kurtogram_Comparison
%
% Outputs:
%    KurtogramComparison.fig
%    KurtogramComparison.png
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

clc
clear
close all
%% Fast Kurtogram
               % sampling frequency

nlevel = 7;     % number of decomposition levels
% Pre-whitening of the signal by AR filtering (optional)
prewhiten = 1;  % (always helpful in detection problems)


%% Add Data and Functions folder
addpath(genpath('../data'))
addpath(genpath('../Functions'))

%% Select Bearings to test and tests
Bearings={'B71','B7','B8','B9','B10','B11'};
Test_name={ 'Levels'; 'Dry'; '   Dry + \newline Interference';...
            'Lub 5%'; 'Lub 100%'};

%% Select Markers
markers={'s','^','o','x','d','<','>'};

%% Load Data and plot mean and uncertainties 
figure(1)
set(gcf,'Position',[ 867   227   698   573])
for k=1:length(Bearings)  %%Read each bearing experiment

    file{1}=['HFAccel_Dry_' Bearings{k}];
    file{2}=['HFAccel_Dry_P1_' Bearings{k} ];
    file{3}=['HFAccel_Lub_min_' Bearings{k} ];
    file{4}=['HFAccel_Lub_full_' Bearings{k} ];
    K_max=[];
    fc_r=[];
    BW_r=[];
    for j=1:length(file) %Iterate over test conditions
    
        load(file{j}); 
        Ns=length(vibR);
        for i=1:Ns  %Iterate over repetitions
            filename=[file{j} '_' Tag{i}];
            speed=mean(rpm_raw{i})/60;
            x=vibR{:,i};
            
            if prewhiten == 1
               x = x - mean(x);
               Na = 100;
               a = lpc(x,Na);
               x = fftfilt(a,x);
               x = x(Na+1:end);		% it is very important to remove the transient of the whitening filter, otherwise the SK will detect it!!
            end

            str=[pwd '\' Bearings{k}];
%             if(~exist(str,'dir'))
%                 mkdir(str)
%             end
            str=[str '\KG_' Tag{i} ];
            Fc=0;   %Not used due lack of envelope analysis
            lv=0;
            %This function returns the max kurtosis in all kurtogram
            [cL{i},levL(i),K_max(i,j),fc_r(i,j),BW_r(i,j)]=Fast_kurtogram_KurtMaX(x,nlevel,Fs,str,Fc,lv);
        end

    end

figure(1)
alpha=0.05;
n=length(K_max);
nu=n-1;
u_rnd=abs(tinv(alpha,nu))*(std(K_max)/sqrt(n));
errorbar([1:2:9]+0.2*(k-1),[-1 mean(K_max)],[0 u_rnd],markers{k},'MarkerSize',10); hold on

end
   

%% Modify Plot

xlim([0 10+0.1*(k-1)])
limits=ylim;
ylim([0 130])
set(gca,'xtick',[1:2:9]+0.1*(k-1))
set(gca,'xticklabel',Test_name)
ylabel('Kurtosis')
xlabel('Conditions')
%%
fontsize = 12;
%Add separator lines  
line(xlim,[3 3],'LineStyle','-','Color',[0.7 0.7 0.7])%[0 0.5 0])
text(0.1,2,'Good','Color',[0 0.5 0],'FontSize',fontsize)
line(xlim,[7 7],'LineStyle','-','Color',[0.7 0.7 0.7])%[0.85 0.85 0])
text(0.1,5,'Satisfactory','Color',[0.8 0.8 0],'FontSize',fontsize)
line(xlim,[10 10],'LineStyle','-','Color',[0.7 0.7 0.7])%[0.9 0.4 0])
text(0.1,9,'Unsatisfactory','Color',[0.9 0.4 0],'FontSize',fontsize)
text(0.1,14,'Unacceptable','Color',[0.9 0 0],'FontSize',fontsize)
%%
for j=2:2:12
line([j j]+0.1*(k-1), ylim,'Color',[0.7 0.7 0.7])
end
%  breakyaxis([130 190]);

legend('Experiment #1','Experiment #2','Experiment #3','Experiment #4',...
        'Experiment #5','Experiment #6')

 
imageName = 'KurtogramComparison';
print(imageName,'-depsc','-r1000')
print(imageName,'-dtiffn','-r1000')
print(imageName,'-dpng','-r1000')%
saveas(gcf,imageName)

%%
beep
pause(1)
beep