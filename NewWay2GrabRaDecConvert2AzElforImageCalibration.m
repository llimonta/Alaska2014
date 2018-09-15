function NewWay2GrabRaDecConvert2AzElforImageCalibration(filename,folder,imageName,xpixmax,ypixmax,lat,long,tUTCfilename)
X= dlmread(filename,'',2);
folder
imageName(2:end-5)
format long g
if class(xpixmax)=='char'
    xpixmax=str2num(xpixmax)
end
if class(ypixmax)=='char'
    ypixmax=str2num(ypixmax)
end

if class(lat)=='char'
    lat=str2num(lat);
end


if class(long)=='char'
    long=str2num(long);
end

tUTC = load(tUTCfilename);
tUTC =  tUTC(median([1:length(tUTC)]))
date=datestr(tUTC,'yyyy-mm-dd HH:MM:SS.FFF');

Ra=X(:,2);
Decl=X(:,3);

% Ra`
lat
long
date

Az=zeros(size(Ra));
El=zeros(size(Ra));

for iii=1:size(Ra,1)
    [Az(iii) El(iii)] = RaDec2AzEl(Ra(iii),Decl(iii),lat,long,date);
end
size(Az)
size(El)
xpixmax
ypixmax
% Saving Right ascension Image
Ra=reshape(Ra,xpixmax,ypixmax);
immagine=figure;
imagesc(Ra);
colorbar;
set(gca,'YDir','normal');
scrittaRa=sprintf([folder,'/Ra',imageName(2:end-5)]);
saveas(immagine,scrittaRa,'eps');
saveas(immagine,scrittaRa,'fig');

% Saving Declination Image
Decl=reshape(Decl,xpixmax,ypixmax);
immagine=figure;
imagesc(Decl);
colorbar;
set(gca,'YDir','normal');
scrittaDecl=sprintf([folder,'/Decl',imageName(2:end-5)]);
saveas(immagine,scrittaDecl,'eps');
saveas(immagine,scrittaDecl,'fig');


% Saving Azimuth Image
Az=reshape(Az,xpixmax,ypixmax);
immagine=figure;
imagesc(Az);
colorbar;
set(gca,'YDir','normal');
scrittaAz=sprintf([folder,'/Az',imageName(2:end-5)]);
saveas(immagine,scrittaAz,'eps');
saveas(immagine,scrittaAz,'fig');

% Saving Elevatio Image
El=reshape(El,xpixmax,ypixmax);
immagine=figure;
imagesc(El);
colorbar;
set(gca,'YDir','normal');
scrittaEl=sprintf([folder,'/El',imageName(2:end-5)]);
saveas(immagine,scrittaEl,'eps');
saveas(immagine,scrittaEl,'fig');

[folder,'/RaDecAzEl',imageName(2:end-5),'.mat']

% Save RA,Decl,Az,Elevation matrix variables
save([folder,'/RaDecAzEl',imageName(2:end-5),'.mat'],'Az','El','Decl','Ra')

end