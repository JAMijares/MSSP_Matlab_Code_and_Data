%Main_RunToFailure_P1 - Creates the figure for the first part of the run to
%failure indicators comparison
%
% Syntax:  Main_RunToFailure_P1
%
% Outputs:
%    IndicatorsComparison_P1.fig
%    IndicatorsComparison_P1.png
%
% Other m-files required: Fast_kurtogram_KurtMaX
% Subfunctions: Fast_kurtogram_KurtMaX
% MAT-files required: 
%       'Accel_Lub_B4'
%       'Accel_Lub_B4_Part2'
%       'Accel_Dry_B4'
%       'Accel_Dry_B4_I5'
%       'Accel_Dry_B4_I5_Part2'
%       'Accel_Dry_B4_I10'
%       'Accel_Dry_B4_I10_Part2'
%       'Accel_Dry_B4_I15'
%       'Accel_Dry_B4_I15_Part2'
%
% Author: Jorge Mijares
% email: jorge.mijares@tamu.edu
% Aug 2019; Last revision: 27-Aug-2019

clc
clear
close all
%% Add Data
addpath(genpath('../Data'))
addpath(genpath('../Functions'))
load Filters
%% Select Markers
% markers={'-s','-^','-o','-x','-d','-<','->'};
Colors={[135,206,250]/256,...     %Blue
        [ 1 1 0],...    %Yellow
        [1 0.5 0],...   %orange
        [ 1 0 0],...  %orangered
        [0.75 0 0]};      %Red
%% Filename
files{1}='Accel_Lub_B4';
files{2}='Accel_Lub_B4_Part2';
files{3}='Accel_Dry_B4';
files{4}='Accel_Dry_B4_I5';
files{5}='Accel_Dry_B4_I5_Part2';
files{6}='Accel_Dry_B4_I10';
files{7}='Accel_Dry_B4_I10_Part2';
files{8}='Accel_Dry_B4_I15';
files{9}='Accel_Dry_B4_I15_Part2';


%%Kurtogram parameters
nlevel = 7;     % number of decomposition levels
% Pre-whitening of the signal by AR filtering (optional)
prewhiten = 1;  % (always helpful in detection problems)

%% Load Data and plot mean and uncertainties 
x=1;
RMS_v=[];
FK=[];
Kurt=[];
CF=[];
fig=figure(1);
hold on

set(gcf,'Position',[300         273         967         350])
left_color = [1 1 1];
right_color = [0.75 0.75 0.75];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
j=1;

for i=1:5
    RMS_v(x)= NaN;
    Kurt(x)= NaN;
    FK(x)=NaN;
    CF(x)= NaN;
    x = x +1;
end
for k=1:length(files)
    load(files{k})
    N=length(t);
    labels{2*k} = Tag{end}(1:7);
    labels{2*k-1} = Tag{1}(1:7);
    ticks_start(k) = x;
    for i=1:N
        vib = vibR{:,i};
        vib=vib-mean(vib);
        vib=filtfilt(hpFilt_1k,vib);
        vib=filtfilt(bsFilt_15k,vib);
        vib=filtfilt(lpFilt_20k,vib);
        RMS_v(x)=rms(vib);
        Kurt(x)=round(kurtosis(vib),2);
        CF(x)=max(abs(vib))/rms(vib);
        %Kurtogram
        if prewhiten == 1
           vib = vib - mean(vib);
           Na = 100;
           a = lpc(vib,Na);
           vib = fftfilt(a,vib);
           vib = vib(Na+1:end);		% it is very important to remove the transient of the whitening filter, otherwise the SK will detect it!!
        end
        Fc=0;   %Not used due lack of envelope analysis
        lv=0;
        [cL,levL,K_max,fc_r,BW_r]=Fast_kurtogram_KurtMaX(vib,nlevel,Fs,files{k},Fc,lv);
        FK(x)=round(K_max,2);
        yyaxis left
        p{k}=bar(x,RMS_v(x),'FaceColor',Colors{1});
        %labels{x}=Tag{i}(1:7);
        x=x+1;
    end
    ticks_end(k) = x;
    if(k~=1 )
        
        for i=1:5
            RMS_v(x)= NaN;
            Kurt(x)= NaN;
            FK(x)=NaN;
            CF(x)= NaN;
            x = x +1;
        end
    end
    
end

%% Plot Kurtosis and CF
ylabel('RMS (g)')
yyaxis right
p1 = plot(Kurt,'ks',...
    'LineWidth',2,...
    'MarkerFaceColor',[.5 1 .5],...
    'MarkerSize',8);
p2 = plot(CF,'ko',...
    'LineWidth',2,...
    'MarkerFaceColor',[0 0 1],...
    'MarkerEdgeColor',[0 0 1],...
    'MarkerSize',8);
p3 = plot(FK,'kd',...
    'LineWidth',2,...
    'MarkerFaceColor',[1 .5 1],...
    'MarkerSize',8);
% legend([p{1}(1) p{2}(1) p{3}(1) p{4}(1) p{5}(1) p1(1) p2(1)],...
%         'RMS Lub',...
%         'RMS Dry',...
%         'RMS Dry + 5mils',...
%         'RMS Dry + 10mils',...
%         'RMS Dry + 15mils',...
%         'Kurtosis','CF','Orientation','Horizontal',...
%         'Location','northoutside')
    

%% Configure Plot
xlabel('Operating Time (minutes)')
ylabel('Kurtosis, CF')


ylim([0 140])
ticks_start(2)=[]; %Eliminate start of second lub file
ticks_end(1)=[];   %Eliminate end of first lub file
ticks=[ticks_start; ticks_end-1];
ticks=ticks(:)';
ticks(end)=ticks(end)-2;
xticks(ticks)
%xticklabels(time)
labels(2)=[];
labels(2)=[];
minutes=[1 3 1 2 1 2 13 14 1 2 4 5 1 2 14 15];
xticklabels(minutes)

xlim([1.5 70.5])
ylim([0 40])

grid on
%% Plot Gray Areas
% 
yl =  ylim;
%fill(ticks(2:3), [yl(2) yl(2)], [0.8 0.8 0.8]);
pauseText={sprintf('Mounting of new lubricated bearing \n + 1 min running to warm up'), ...
           sprintf('Pause to remove lubricant  \n + 12hrs to wait for drying \n + 1 min running to warm up'), ...
           sprintf('Pause to increase\n interference to 5 mils\n + 1 min running to warm up'), ...
           sprintf('Continue operation')};
 ticks=[1 ticks];      
for k=1:4
    bar_xlim = ticks(2*k-1:2*k);
    bar_xlim(2) = bar_xlim(2)-0.5;
    bar_xlim(1) = bar_xlim(1)+0.5;
    [X,Y] = hatch_coordinates( bar_xlim , yl , 0.33 , 0.33 ) ; %// this return coordinates to plot a hatch pattern
    plot(X,Y,'Color',[.6 .6 .6],'Marker','none','linewidth',1,'LineStyle',':')
    fill([bar_xlim(1) bar_xlim(1) bar_xlim(2) bar_xlim(2)], [yl(1) yl(2) yl(2) yl(1)], [0.8 0.8 0.8],...
        'EdgeColor','none','FaceAlpha',.3,'LineStyle','None');
    % Add Text to gray regions
    h=text(0.5*(bar_xlim(2)+bar_xlim(1)),yl(2)/2,pauseText{k},...
        'HorizontalAlignment','center','VerticalAlignment','middle');
    set(h,'Rotation',90);
end



%% Add Legends
 legend([p{1}(1) p1(1) p2(1) p3(1)],'RMS','Kurtosis','CF','Enhanced Kurtosis','Orientation','Horizontal',...
          'Location','northoutside')
%%
imageName = 'IndicatorsComparison_P1';
print(imageName,'-depsc','-r1000')
print(imageName,'-dtiffn','-r1000')
print(imageName,'-dpng','-r1000')%
print(imageName, '-djpeg','-r1000')   % save as COLOR pdf file
saveas(gcf,imageName)