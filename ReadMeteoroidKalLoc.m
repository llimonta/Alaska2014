function []=ReadMeteoroidKalLoc(imageFolder,Frame,day,kkk,iii)
% kkk = number meteoroids from event at night
% iii = number frame we are interested in for position
load('/home/limo/PHD/PokerFlat2014/MatlabCode/Data/20140331Meteoroids.mat');
% Writing the position for future use
xpixel = round(Meteoroid{kkk}.Loc(iii,1));
ypixel = round(Meteoroid{kkk}.Loc(iii,2));
fileID = fopen([imageFolder,'/Frame',num2str(Frame),day,'XpixelMeteoroidPosition.txt'],'w');
fprintf(fileID,'%3d',[xpixel]);
fileID = fopen([imageFolder,'/Frame',num2str(Frame),day,'YpixelMeteoroidPosition.txt'],'w');
fprintf(fileID,'%3d',[ypixel]);
fclose(fileID);
