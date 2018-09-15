function [star] = GetMagAndPixelsForBrightStars(inputFile,Nxpix,Nypix,limiting_mag,frame_data_name,star_day_saved,tUTCfilename,lat,long)

%   Input
%   inputFile= absolute path of where the astronometry bright star path txt
%   baseline_data_file= path to where the data for the given frame is being
%   Nxpix,Nypix Number of xpixels and ypixels in the image
%   analyzed
%   file is
%   tycho.fits = absolute path of where the tycho2 fits table is 
%       Output
% star.NameTycho1=First Digit of the tycho denomination
% star.NameTycho2=Second Digit of the tycho denomination
% star.NameTycho3=Third Digit of the tycho denomination
% star.Xpixel,Ypixel= location of the star at the given frame
% star.Frame = where we are detecting the star
% 
%   NB to visualize a specific array for all the stars do the following:  
%   cellfun(@(S) S.BTmag, star)   
%
%
%
% Load the Tycho-2 fits columns with the tychos numbers so that we can get
% the BT, VT magnitudes (columns from 25 to 30) 

TychosNameTable=fitsread('/home/limo/PHD/astrometry.net-0.73/catalogs/tycho2.fits','binarytable','TableColumns',[1:5,25:30]);
SpectralData=importdata('/home/limo/PHD/PokerFlat2014/catalog.dat');
tUTC = load(tUTCfilename);
tUTC =  tUTC(median([1:length(tUTC)]));
% AllRatios=[];

if nargin<6,star_day_saved='babbo'; end

% if exist(star_day_saved,'file') ==2
% %     load(star_day_saved);
%     flagStarsCatalaogExist=true;
% else
     flagStarsCatalaogExist=false;
% end

% if you cannot load the entire table uncomment the following and comment
% out the line above
% TychosNameTable=fitsread('/home/limo/PHD/Astrometry/astrometry.net-0.50/catalogs/tycho2.fits','binarytable','TableColumns',[1:3]);

% The structure of the above fits file is the following 

%  C1   C2   C3  C4 C5   C6        C7      C8        C9        C10          
% TYC1 TYC2 TYC3 RA DEC MEAN_RA MEAN_DEC SIGMA_RA SIGMA_DEC SIGMA_MEAN_RA     {1×1 struct}    {1×1 struct}    {1×1 struct}    {1×1 struct}    {1×1 struct}    {1×1 struct}    

%      C11         C12    C13     C14         C15          C16      C17  
% SIGMA_MEAN_DEC PM_RA  PM_DEC SIGMA_PM_RA SIGMA_PM_DEC EPOCH_RA EPOCH_DEC

%   C18              C19        C20     C21            C22            C23
% EPOCH_MEAN_RA EPOCH_MEAN_DEC NOBS GOODNESS_MEAN_ GOODNESS_MEAN_ GOODNESS_PM_RA 

%   C24             C25         C26      C27         C28       C29     C30
% GOODNESS_PM_DE   MAG_BT   SIGMA_MAG_BT MAG_VT SIGMA_MAG_VT  MAG_HP  SIGMA_MAG_HP 

%  C31         C32       C33      C34    C35 
% PROX    CORRELATION HIPPARCOS_I CCD    FLAGS 


infile=fileread(inputFile);
[tempFrame]=regexp(infile,'[*](\d+)[*]','match');
[tempFrame]=regexp(tempFrame,'(\d+)','match'); %excluding the '*' characters
Frame=(char(tempFrame{1}));
if flagStarsCatalaogExist==false
counter_star=0;
end
        TychoNames=regexp(infile,'(\d+)-(\d+)-(\d+)','match'); %looking in the file where the tycho-2 name is since we know that the format is number-number-number we look for that!
        [dummy,BegXp]=regexp(infile,'"], "pixelx": '); %looking in the file where the pixelx word begins and end. Looking at the inputFile where xpixel begins is where tycho-2 reference number ends
        BegXp=BegXp+1; %Adding one because the above line gives the position to the space
        [EndXp,BegYp]=regexp(infile,', "pixely": '); %looking in the file where the pixely word begins and end. Looking at the inputFile where ypixel ends is where xpixel number ends
        EndXp=EndXp-1; %Subtractiong minus one because the above line gives location of the ' , ' character
        BegYp=BegYp+1;
       
       
        for iii=1:length(TychoNames)
         
            Name=TychoNames(iii);
            [tempNumbers]=regexp(Name,'(\d+)','match');
            Tycho1=str2double(tempNumbers{1}(1));
            Tycho2=str2double(tempNumbers{1}(2));
            Tycho3=str2double(tempNumbers{1}(3));
            
            dummy1=find(TychosNameTable{1}(:)==Tycho1);
            dummy2=find(TychosNameTable{2}(dummy1)==Tycho2);  % This allow me to find in which row is the star located in the table
            
            StarTemp=[TychosNameTable{1}(dummy1(dummy2)),TychosNameTable{2}(dummy1(dummy2)),TychosNameTable{3}(dummy1(dummy2)), ...
                      TychosNameTable{6}(dummy1(dummy2)),TychosNameTable{7}(dummy1(dummy2)),TychosNameTable{8}(dummy1(dummy2)), ...
                      TychosNameTable{9}(dummy1(dummy2)),TychosNameTable{10}(dummy1(dummy2)),TychosNameTable{11}(dummy1(dummy2)),...
                      TychosNameTable{4}(dummy1(dummy2)),TychosNameTable{5}(dummy1(dummy2))];
            %
            %If you cannot load the entire table comment out the above line
            %and uncomment the following two
%             StarDummy=fitsread('/home/limo/PHD/Astrometry/astrometry.net-0.50/catalogs/tycho2.fits','binarytable','TableColumns',[1:3,25:30],'TableRows',dummy1(dummy2));
%             StarTemp=[StarDummy{:}];
            % star temp has the following format:
              %C1    C2  C3   C4       C5         C6       C7            C8       C9         C10              C11
            % TYC1 TYC2 TYC3 MAG_BT SIGMA_MAG_BT MAG_VT SIGMA_MAG_VT  MAG_HP SIGMA_MAG_HP     RA              DEC
             
            xpix=str2double(infile(BegXp(iii):EndXp(iii)));
            ypix=str2double(infile(BegYp(iii):BegYp(iii)+12));
            if (StarTemp(4)<limiting_mag) && (StarTemp(6)<limiting_mag) && StarTemp(8)<limiting_mag && (xpix>10 && xpix<(Nxpix-10)) && (ypix>10 && ypix<(Nypix-10))
%             StarTemp(4)
%             StarTemp(6)
                if flagStarsCatalaogExist==0
                    counter_star=counter_star+1;
                     star{counter_star}.Xpixel=[];
                     star{counter_star}.Ypixel=[];
                     star{counter_star}.Luminosity=[];
                     star{counter_star}.Frame=[];
                     star{counter_star}.AirMass=[];
                     star{counter_star}.SNR=[];
                else
                    %
                    % find counter_star
                    %
                    ExistingStarNames=cellfun(@(S) S.Name, star,'UniformOutput',false);
                    LocationExistingStar=regexp(ExistingStarNames,TychoNames(iii),'match');
                    counter_star=find((cellfun('isempty',LocationExistingStar))==0);
                    if isempty(counter_star)
                        counter_star=length(star)+1;
                        star{counter_star}.Xpixel=[];
                        star{counter_star}.Ypixel=[];
                        star{counter_star}.Luminosity=[];
                        star{counter_star}.SNR=[];
                        star{counter_star}.Frame=[];
                        star{counter_star}.AirMass=[];
                    end
                end
                  counter_star;
                star{counter_star}.Name=char(TychoNames(iii));
                star{counter_star}.Tycho1=StarTemp(1);
                star{counter_star}.Tycho2=StarTemp(2);
                star{counter_star}.Tycho3=StarTemp(3);
                star{counter_star}.Xpixel=[star{counter_star}.Xpixel,round(str2double(infile(BegXp(iii):EndXp(iii))))];
                star{counter_star}.Ypixel=[star{counter_star}.Ypixel,round(str2double(infile(BegYp(iii):BegYp(iii)+12)))];
                if StarTemp(4)==0 && StarTemp(6)==StarTemp(4)
                    StarTemp(4)=NaN;
                    StarTemp(6)=NaN;
                end
                star{counter_star}.BTmag=StarTemp(4);
                star{counter_star}.BTsigma=StarTemp(5);
                star{counter_star}.VTmag=StarTemp(6);
                star{counter_star}.VTsigma=StarTemp(7);
                star{counter_star}.HPmag=StarTemp(8);
                star{counter_star}.HPsigma=StarTemp(9);
                star{counter_star}.Frame=[star{counter_star}.Frame,str2double(Frame)];
                % correcting VT and BT mag into standard BV mag, http://www.aerith.net/astro/color_conversion.html
                star{counter_star}.Vmag = star{counter_star}.VTmag + 0.00097 - 0.1334 * (star{counter_star}.BTmag - star{counter_star}.VTmag) + ...
                    0.05486 * (star{counter_star}.BTmag - star{counter_star}.VTmag)^2 - 0.01998 * (star{counter_star}.BTmag - star{counter_star}.VTmag)^3 ;
                if star{counter_star}.BTmag - star{counter_star}.VTmag>0.5
                    star{counter_star}.Bmag = star{counter_star}.Vmag + (star{counter_star}.BTmag - star{counter_star}.VTmag) - 0.007813 * (star{counter_star}.BTmag - star{counter_star}.VTmag) ...
                        - 0.1489 * (star{counter_star}.BTmag - star{counter_star}.VTmag)^2 + 0.03384 * (star{counter_star}.BTmag - star{counter_star}.VTmag)^3;
                else
                    star{counter_star}.Bmag = star{counter_star}.Vmag + (star{counter_star}.BTmag - star{counter_star}.VTmag) - 0.006 - 0.1069 * (star{counter_star}.BTmag - star{counter_star}.VTmag) ...
                        + 0.1459 * (star{counter_star}.BTmag - star{counter_star}.VTmag)^2;
                end
                star{counter_star}.B_V=star{counter_star}.Bmag-star{counter_star}.Vmag;
                % Adding air mass
                % From pag.29
 %of A practical Guide to LightcuvePhotometry
                % and Analysis
                Ra = StarTemp(10);
                Decl = StarTemp(11);
                tempo=datestr(tUTC,'yyyy-mm-dd HH:MM:SS.FFF');
                [~, El]=RaDec2AzEl(Ra,Decl,lat,long,tempo);
                El = El*pi/180; %putting it in radians
                zenith = pi/2-El;
                star{counter_star}.AirMass = [star{counter_star}.AirMass,1/cos(zenith)];
                
                % Adding the stellar spectrum
                counter_match = find(strncmp( SpectralData,['TYC ',num2str(StarTemp(1)),' 0',num2str(StarTemp(2))],14));
                if isempty(counter_match) == 0
                    star{counter_star}.Spectrum = SpectralData{counter_match}(94);
                    clear counter_match
                else
                    star{counter_star}.Spectrum = 'Z';
                end                
                 star{counter_star}.Spectrum
                [luminosity,SNR]=StarPixelCount(xpix,ypix,star{counter_star}.VTmag,Nxpix,Nypix,frame_data_name, star{counter_star}.Spectrum);
%                 AllRatios=[AllRatios;ratios];
                star{counter_star}.Luminosity=luminosity;
                star{counter_star}.SNR=SNR;
                star{counter_star}.SNR;
                
            end
        end
        save(star_day_saved,'star')%RBrEV4Ixon300320141-100000.mat','RBrEv4')

end