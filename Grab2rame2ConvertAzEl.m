function Grab2rame2ConvertAzEl(BegFrame,NightToBeAnalyzed,day)
%
% Beg Frame = beginning of Frame (char)
% NightToBeAnalyzed : day to be analyzed (char)
% Day to be analyzed (char)

if class(BegFrame)=='char'
    FrameInit=str2num(BegFrame);
end

[data_baseline,~,tUT] = rawDMCreader(NightToBeAnalyzed,512,512,1,1,FrameInit:FrameInit+350,0,[100,1100],'auto','auto');
data_baseline=mean(double(data_baseline),3);
fitswrite(data_baseline,['/home/limo/PHD/PokerFlat2014/AzDecImages/Frame',BegFrame,day,'.fits'])