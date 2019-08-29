clc
clear
fileList = dir('*.png');

n = length(fileList);
%%
for i = 1:n
    fileName = fileList(i).name; % your FILE NAME as string
    [filepath,name,ext] = fileparts(fileName)
    A = imread(fileName);
    set(gcf,'visible','off') %suppress figure
    image(A);                
    axis image               % resolution based on image
    axis off                 % avoid printing axis 
    set(gca,'LooseInset',get(gca,'TightInset')); % removing extra white space in figure
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    %saveas(gcf,name,'epsc');   
    %saveas(gcf,name,'tiffn');   
    print(fig, name, '-dpdf','-r600')   % save as COLOR pdf file
end
