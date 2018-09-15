function MeteorMagnitudeEstimation(folder,day,file_saved_stars,RaDecAzElfile,inizio,fine)

format long g
% Getting the stars for magnitude calibration of our meteor
load(file_saved_stars); %this file contain a thing that is star{1}.Etc
load(RaDecAzElfile);

% Prepare for estimating parameters for the estimate of the meteor
% Magnitude
Vmag = zeros(length(star),1);Bmag=zeros(length(star),1);Lum=zeros(length(star),1); B_V=zeros(length(star),1);
AirMass=zeros(length(star),1);
for iii=1:length(star) 
    Vmag(iii)=star{iii}.Vmag;
    Bmag(iii)=star{iii}.Bmag;
    Lum(iii)=star{iii}.Luminosity(1);
    B_V(iii)=star{iii}.B_V;
    AirMass(iii)=star{iii}.AirMass(1);
end
LumPrimed = -100^(1/5)*log10(Lum);
X = [ones(length(star),1),LumPrimed,AirMass]; %(aka x, what we regress on
Y = Vmag; % aka y, our target
Y1 = Bmag;
% We fit 
%Finding the fitting coefficients for the meteor direct method
[b,~,~,rint] = regress(Y,X,0.05);
contain0 = (rint(:,1)<0 & rint(:,2)>0);
idx = find(contain0==false);
Y(idx)=NaN;   
[b,~,~,~,~] = regress(Y,X,0.05);
%Finding fitting coefficients for BMAG
[bbmag,~,~,rintbmag] = regress(Y1,X,0.05);
contain0alt = (rintbmag(:,1)<0 & rintbmag(:,2)>0);
idxbmag = find(contain0alt==false);
Y1(idxbmag)=NaN;
[bbmag,~,~,~,~] = regress(Y1,X,0.05);



% Loading our Meteoroids and their parameters

%This load a structure called Meteoroid
% Meteoroid{iii}.MeteorIntensity
% Meteoroid{iii}.SNR
% Meteoroid{iii}.tUT
% Meteoroid{iii}.KalLoc --> x,y position of our meteoroid
['/home/limo/PHD/PokerFlat2014/MatlabCode/Data/31-03-2014luminosity_by_hand_for_orbit_determination_meteor_number',[num2str(inizio),'_',num2str(fine)],'.mat']
load(['/home/limo/PHD/PokerFlat2014/MatlabCode/Data/31-03-2014luminosity_by_hand_for_orbit_determination_meteor_number',[num2str(inizio),'_',num2str(fine)],'.mat'])
load(['/home/limo/PHD/PokerFlat2014/MatlabCode/Data/31-03-2014location_by_hand_for_orbit_determination_meteor_number',[num2str(inizio),'_',num2str(fine)],'.mat'])
meteor_instrumental_luminosity = luminosity;%Meteoroid{meteoroid_number}.MeteorIntensity;
meteor_time = tUT;
meteor_location = position_meteor; % given as meteor_location(1,:) = xpixel,ypixel
VmagMeteor=zeros(length(meteor_instrumental_luminosity),1);
BmagMeteor=VmagMeteor;
% elevation = load(exlevation_file) %loading the pixel value of elavation, which is constant over time.

for iii=1:length(meteor_instrumental_luminosity)
    if (meteor_location(iii,1))>0.5 && (meteor_location(iii,2))>0.5 && meteor_instrumental_luminosity(iii)>50 && (meteor_location(iii,1))<512.49 && (meteor_location(iii,2))<512.49
        %written this way because first is row, then is column, hence preserving x-y pixel coordinates
        meteor_location(iii,2);
        meteor_location(iii,1);
        meteor_El = El(round(meteor_location(iii,2)),round(meteor_location(iii,1)))*pi/180;
        meteor_zenith = pi/2-meteor_El;
        AirMassMeteor = 1/cos(meteor_zenith);
        meteor_instrumental_luminosity(iii);
        VmagMeteor(iii) = b(1)+b(2)*(-100^(1/5)*log10(meteor_instrumental_luminosity(iii)))+b(3)*AirMassMeteor;
        VmagMeteor(iii)
%         VmagTLS(iii)=bTLS(1)+bTLS(2)*(-100^(1/5)*log10(meteor_instrumental_luminosity(iii)))+bTLS(3)*AirMassMeteor;
        BmagMeteor(iii) = bbmag(1)+bbmag(2)*(-100^(1/5)*log10(meteor_instrumental_luminosity(iii)))+bbmag(3)*AirMassMeteor;
        BmagMeteor(iii)
%         VmagTLSalt(iii) = fzero(@(VmagAlt) -100^(1/5)*log10(meteor_instrumental_luminosity(iii))-bTLSalt(1)-bTLSalt(2)*VmagAlt-bTLSalt(3)*AirMassMeteor,VmagMeteor(iii))
    else
        VmagMeteor(iii) = NaN;
        BmagMeteor(iii) = NaN;
%         VmagTLS(iii) = NaN;
    end
end
[folder,'/BandVmagMeteorHand',[num2str(inizio),'_',num2str(fine)],day,'.mat']
save([folder,'/BandVmagMeteorHand',[num2str(inizio),'_',num2str(fine)],day,'.mat'],'VmagMeteor','BmagMeteor','meteor_time');
end


function [Xnew,coeffs]=TotalLeastSquare(X,Y)

%Total Last Square
aug_mat=[X,Y];
[~,~, V] = svd(aug_mat,0);
[~, n_col_X]=size(X);
VXY=V(1:n_col_X,n_col_X+1:end);
VYY=V(n_col_X+1:end,n_col_X+1:end);
coeffs = -VXY*inv(VYY);
Xtyt = - [X,Y] * V(:,n_col_X+1) * V(:,n_col_X+1)'; % [ X-tilde  y-tilde ] eq’n (15)
%Xtyt_norm_frob = norm(Xtyt,'fro');% Frobeneus norm of [Xt yt] (5)
Xt = Xtyt(:,1:n_col_X);
Xnew = Xt+X;
end
