%AddAnnotations - Add the annotations to the CascadeAnalysis
%
% Syntax:  AddAnnotations
%
% Outputs:
%    CascadeAnalysis.fig
%    CascadeAnalysis.png
%
% MAT-files required: 
%       Bearing [BearingID].fig
%
% Author: Jorge Mijares
% email: jorge.mijares@tamu.edu
% Aug 2019; Last revision: 27-Aug-2019
clc
clear
f=openfig("Bearing B7")  %Select Bearing ID
%% Annotations
% Delete annotations
delete(findall(gcf,'type','annotation'))
fontsize = 11;
% Increase in stiffness

x = [0.2 0.38];
y = [0.65 0.48];
str = sprintf('Shift due to\nincreased stiffness');
annotation('textarrow',x,y,'String',str,'HorizontalAlignment','center','FontSize',fontsize)
% Surface roughness excitation
x = [0.2 0.26];
y = [0.15 0.33];
str = sprintf('Surface roughness\nexcitation');
annotation('textarrow',x,y,'String',str,'HorizontalAlignment','center','FontSize',fontsize)
% 2nd arrow
x(2) = 0.41;
y(2) = 0.35;        
annotation('textarrow',x,y)
% 3rd arrow
x(2) = 0.49;
y(2) = 0.37;
annotation('textarrow',x,y)

% High pass filter
x = [0.4 0.35];
y = [0.95 0.75];
str = sprintf('High pass filter effect');
annotation('textarrow',x,y,'String',str,'HorizontalAlignment','center','FontSize',fontsize)
% Band pass filter
x = [0.8 0.73];
y = [0.9 0.64];
str = sprintf('Band stop filter for\nVFD frequency');
annotation('textarrow',x,y,'String',str,'HorizontalAlignment','center','FontSize',fontsize)
% FK central frequency
x = [0.882 0.838];
y = [0.78 0.62];
str = sprintf('Central frequency\ndetected by FK');
annotation('textarrow',x,y,'String',str,'HorizontalAlignment','center','FontSize',fontsize)
% Rotor Resonances
x = [0.52 0.46];
y = [0.88 0.78];
str = sprintf('System\nresonances');
annotation('textarrow',x,y,'String',str,'HorizontalAlignment','center','FontSize',fontsize)
x(2) = 0.65;
y(2) = 0.7;
annotation('textarrow',x,y)
%%
imageName ='CascadeAnalysis';
print(imageName,'-depsc','-r1000')
print(imageName,'-dtiffn','-r1000')
print(imageName,'-dpng','-r1000')%
saveas(gcf,imageName)
