function GrabFrame2fits(BegFrame,NFRAMES,NightToBeAnalyzed,day,imageFolder,imageName,dataName)
%
% Beg Frame = beginning of Frame (char)
% NightToBeAnalyzed : day to be analyzed (char)
% Day to be analyzed (char)

BegFrame
NightToBeAnalyzed
day

if strcmp(class(BegFrame),'char')
    FrameInit=str2num(BegFrame);
else
    FrameInit=BegFrame
end

if strcmp(class(NFRAMES),'char')
    Nframes=str2num(NFRAMES);
else
    Nframes=NFRAMES;
end

%cd('~/PHD/PokerFlat2014/MatlabCode/')
FrameInit-round(Nframes/2)
FrameInit+round(Nframes/2)
[data_baseline,~,tUT] = rawDMCreader(NightToBeAnalyzed,512,512,1,1,FrameInit-round(Nframes/2):FrameInit+round(Nframes/2),0,[100,1100],'auto','auto');
data_baseline=mean(double(data_baseline),3);
% Writing the Fits file
fitswrite(data_baseline,[imageFolder,imageName])   
figura1=figure;
% Saving the matlab file
save([imageFolder,dataName],'data_baseline')
% Writing the time tUTC for future use
fileID = fopen([imageFolder,'/Frame',num2str(BegFrame),day,'Time.txt'],'w');
fprintf(fileID,'%24.24f\n',tUT);
fclose(fileID);
