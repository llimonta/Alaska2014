close all
clear all
% Compute bolometric output power for zero magnitude object in specified
% band

% Reference Star Vega, %  True Zero % Units: erg/sec/cm^2/A
load('Uspectrum'); Uspectrum=flipud(Uspectrum); Uspectrum(Uspectrum(:,2)<0,2)=0;
load('Bspectrum'); Bspectrum=flipud(Bspectrum); Bspectrum(Bspectrum(:,2)<0,2)=0; Bspectrum(Bspectrum(:,1)>540,2)=0;
load('Vspectrum'); Vspectrum=flipud(Vspectrum); Vspectrum(Vspectrum(:,2)<0,2)=0;
load('Rspectrum'); Rspectrum=flipud(Rspectrum); Rspectrum(Rspectrum(:,2)<0,2)=0;
VegaData=fitsread('alpha_lyr_stis_008.fits','binarytable');
VegaSpectrum = VegaData{1}; %Anstrongms
VegaFlux = VegaData{2}; % Erg/(s*cm^2*A)

% Find common values of flux;
VegaSpectrumBands=VegaSpectrum(VegaSpectrum>min(Uspectrum(:,1))*10 & VegaSpectrum<max(Uspectrum(:,1))*10);
VegaFlux = VegaFlux(VegaSpectrum>min(Uspectrum(:,1))*10 & VegaSpectrum<max(Uspectrum(:,1))*10);
%Interpolate U,B,V,R,spectrum into VegaValues
Ufilter=interp1(Uspectrum(:,1)*10,Uspectrum(:,2),VegaSpectrumBands,'spline');
Bfilter=interp1(Bspectrum(:,1)*10,Bspectrum(:,2),VegaSpectrumBands,'spline');
Vfilter=interp1(Vspectrum(:,1)*10,Vspectrum(:,2),VegaSpectrumBands,'spline');
Rfilter=interp1(Rspectrum(:,1)*10,Rspectrum(:,2),VegaSpectrumBands,'spline');
%Camera QE
BVBFQE=csvread('/home/limo/Documents/Thesis/CodeForPlots/CameraQEBvBvf.csv',1);

lambda=Uspectrum(:,1)*1e-8*10; % in cm

% Constants
k = 1.38064852*10^(-16); % erg/K
h = 6.6260755*10^(-27); % erg*s
c = 2.99792458*10^10; % cm/s
sigma = 2*pi^5*k^4/(15*h^3*c^2);% erg/(cm^2*s*K^4)
Height = 10^7; % cm

BV=[];
VR=[];
BR=[];
TOTAL=[];
temp=[13200:5:14850];
% ffff=figure('Visible', 'off'); 
for T=4500
    T
    planck=2*h*c^2./lambda.^5*1./(exp(h*c./(lambda.*k*T))-1);
    figure(1)
    yyaxis left
    plot(lambda*1e8,planck./max(planck),'k','LineWidth',3)    
    hold on
    plot(VegaSpectrumBands,VegaFlux./max(VegaFlux),'.y','MarkerSize',15)
    yyaxis right
    plot(VegaSpectrumBands,Bfilter,'-.b',VegaSpectrumBands,Vfilter,'-.g',VegaSpectrumBands,Rfilter,'-.r','LineWidth',3)
    hold on
    plot(BVBFQE(:,1)*1e1,BVBFQE(:,2),'color',[0.8500, 0.3250, 0.0980],'LineWidth',3)
    legend('Planck','Vega','B-Filter','V-filter','R-Filter','Camera QE')
    xAX = get(gca,'XAxis');set(xAX,'FontSize', 12);yAX = get(gca,'YAxis');set(yAX,'FontSize', 12);
    xlabel(['Wavelength [',char(197),']'],'FontSize',20)
    ylabel('Efficiency','FontSize',20);hold on
    title('Calibration Properties','FontSize',20,'FontWeight','normal')
    yyaxis left
    ylabel('Normalized Spectrum','FontSize',20);hold on
    % Following outputs are all in Watts and correct, compare with Ceplecha
    % 1998 and Weryk 2012
    IboloU=4*sigma*T^4*Height^2*trapz(VegaSpectrumBands,Ufilter.*VegaFlux)/(pi*trapz(lambda,Uspectrum(:,2).*planck))*10^(-7);
    IboloB=4*sigma*T^4*Height^2*trapz(VegaSpectrumBands,Bfilter.*VegaFlux)/(pi*trapz(lambda,Bspectrum(:,2).*planck))*10^(-7)
    IboloV=4*sigma*T^4*Height^2*trapz(VegaSpectrumBands,Vfilter.*VegaFlux)/(pi*trapz(lambda,Vspectrum(:,2).*planck))*10^(-7)
    IboloR=4*sigma*T^4*Height^2*trapz(VegaSpectrumBands,Rfilter.*VegaFlux)/(pi*trapz(lambda,Rspectrum(:,2).*planck))*10^(-7)
    BV=[BV,abs(IboloB-IboloV)];
    VR=[VR,abs(IboloV-IboloR)];
    BR=[BR,abs(IboloB-IboloR)];
    figure(2)
    TOTAL=[TOTAL,sqrt(abs(IboloB-IboloR)^2+abs(IboloV-IboloR)^2+abs(IboloB-IboloV))];
    plot(T,IboloB,'.b',T,IboloV,'.g',T,IboloR,'.r','MarkerSize',15)
%     set(ffff,'Visible', 'off'); 
    hold on
% pause 
end
axis([min(temp) max(temp) 500 5000])
xlabel('Temperature [K]');ylabel('Power [W]')
legend('B band','V band','R band')
figure(); semilogy(temp,BV,temp,VR,temp,BR,temp,TOTAL,'LineWidth',3);legend('BV difference','VR difference','BR difference','Total Discrepancy')
xlabel('Temperature [K]');ylabel('Difference [W]')