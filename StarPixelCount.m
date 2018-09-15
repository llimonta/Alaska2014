function [luminosity, SNR,ratios]=StarPixelCount(xpixel,ypixel,mag,Nx,Ny,framen_data_name,spectrum)
% [luminosity, SNR,ratios]=StarPixelCount(xpixel,ypixel,mag,Nx,Ny,framen_data_name,spectrum)
% close all
load(framen_data_name)

Xcenter=xpixel;
Ycenter=ypixel;
NpixelX=Nx;
NpixelY=Ny;

% for MaskRadius =1:10
MaskRadius = 4.5;
[X,Y]=meshgrid(-(Xcenter-1):(NpixelX-Xcenter),-(Ycenter-1):(NpixelY-Ycenter));
c_mask=((X.^2+Y.^2)<=MaskRadius^2);
temp=data_baseline(c_mask);

if Xcenter>26 && Ycenter>26 && Xcenter<485 && Xcenter< 485
    annulus1=((X.^2+Y.^2)<=25*MaskRadius^2);
    annulus2=((X.^2+Y.^2)<=9*MaskRadius^2);
    annulus= data_baseline(logical(abs(annulus1-annulus2)));
    night_sky_baseline = mode(annulus(annulus<prctile(annulus,50)));
    temp = temp(temp(:)>night_sky_baseline);
    night_sky=length(temp(temp(:)>0))*night_sky_baseline;
    luminosity=sum(temp(temp(:)>0))-night_sky;
else
% %     temp(temp<mean(data_baseline(data_baseline(:)<prctile(data_baseline(:),90))))=0;
%     length(temp);
%     length(temp(temp~=0));
%     night_sky=(length(temp(temp~=0))*mean(data_baseline(data_baseline(:)<prctile(data_baseline(:),90))))
%     luminosity=sum(temp(:))-night_sky;
% temp(temp<mean(data_baseline(data_baseline(:)<prctile(data_baseline(:),90))))=0;
% length(temp);
% length(temp(temp~=0));
% luminosity=sum(temp(:))-(length(temp(temp~=0))*mean(data_baseline(data_baseline(:)<prctile(data_baseline(:),90))));
     luminosity = NaN;
     SNR = NaN;
end
% pause
SNR = luminosity/sqrt(luminosity)*sqrt(2);
end 

%
% %
% % If you are cundcting the aperture size determination, use this code for
% % plotting it all
% %
% color = containers.Map;
% color('O')=[0.9100    0.4100    0.1700];color('B')='c';color('A')='b';color('F')='y';color('G')='g';color('K')='r';color('M')='m';color('Z')='k';
% magnitude = containers.Map;
% magnitude('5') = 'y'; magnitude('5.5')=[0.9100    0.4100    0.1700]; magnitude('6')='r'; magnitude('6.5')='m';magnitude('7') = 'g'; magnitude('7.5')='c'; magnitude('8') = 'b'; magnitude('8.5')='k';
% magnitude('9')='k';magnitude('NaN')='k'
% 
% magnitude('5') = []; magnitude('5.5')=[]; magnitude('6')=[]; magnitude('6.5')=[];magnitude('7') = []; magnitude('7.5')=[]; magnitude('8') = []; magnitude('8.5')=[];
% magnitude('9')=[];magnitude('NaN')='k'
% lunghezza = 10;
% yellow = [1, 0, 0];
% blue = [0 0 1];
% colors_p = [linspace(yellow(1),blue(1),lunghezza)', linspace(yellow(2),blue(2),lunghezza)', linspace(yellow(3),blue(3),lunghezza)'];
% colors_p= varycolor(10);
% colors_p = flipud(jet(10));
% dummy_counter=0;
% for aaa=5:0.5:9
%     dummy_counter=dummy_counter+1;
%     magnitude(num2str(aaa))=colors_p(dummy_counter,:);
% end
% magnitude('NaN')='k'
% 
% ratios=luminosity(1:end-1)./luminosity(2:end);
% if randn(1)>.4
% 
%     figure(1)
%     plot(1:10,luminosity,'.','MarkerSize',15,'color',color(spectrum))
%     hold on
%     plot(1:10,luminosity,'LineWidth',3,'color',color(spectrum))
%     xAX = get(gca,'XAxis');set(xAX,'FontSize', 15);yAX = get(gca,'YAxis');set(yAX,'FontSize', 15);
%     xlabel('Aperture size in pixels','FontSize',20)
%     ylabel('Raw count','FontSize',20)
%     title('Count as a function of aperture','FontWeight','Normal','FontSize',20)
% %     legend('O','B','A','F','G','K','M')
%     h=zeros(8,1);h(1) = plot(NaN,NaN,'color',[0.9100    0.4100    0.1700]); h(2) = plot(NaN,NaN,'c'); h(3) = plot(NaN,NaN,'b'); h(4) = plot(NaN,NaN,'y'); h(5) = plot(NaN,NaN,'g'); h(6) = plot(NaN,NaN,'r'); h(7) = plot(NaN,NaN,'m');h(8) = plot(NaN,NaN,'color','k');
%     lgd = legend(h, 'O','B','A','F','G','K','M','N/A','Location','Best');
%     lgd.FontSize = 12;
% 
%     for zzz=6:10
%         
%         figure(zzz-4)
%         plot(1:zzz,luminosity(1:zzz)./max(luminosity(1:zzz)),'.','MarkerSize',15,'color',magnitude(num2str(round(mag*2)/2)))
%         hold on
%         plot(1:zzz,luminosity(1:zzz)./max(luminosity(1:zzz)),'LineWidth',2,'color',magnitude(num2str(round(mag*2)/2)))
% %         colorbar(jet(10))
%         xlabel('Aperture size in pixels')
%         ylabel('Raw Count(i)/Raw Count(end)','Interpreter','latex')
%         title('Percentage of total luminosity as function of aperture','FontWeight','Normal')
%         hold on
%     end
%     figure(zzz+1)
%     plot(1:10,luminosity(1:10),'.','MarkerSize',15,'color',magnitude(num2str(round(mag*2)/2)))
%     hold on
%     plot(1:10,luminosity(1:10),'LineWidth',3,'color',magnitude(num2str(round(mag*2)/2)))
%     xlabel('Aperture size in pixels')
%     ylabel('Raw count')
%     title('Count as a function of aperture','FontWeight','Normal')
%     %         colorbar(jet(10))
%     hold on
% %     figure(3)
% %     plot(1:6,luminosity(1:6)./max(luminosity(1:6)),'.','MarkerSize',10,'color',magnitude(num2str(round(mag*2)/2)))
% %     hold on
%     % Xcenter
%     % Ycenter
%     % pause
% end
% 
% 
% 
% %   Following part is for fixing axis and colorbars, use it in the main
% %   matlab window not here
% for zzz=6:10
% figure(zzz-4)
% hold on
% colormap(jet(10))
% colorbar('Ticks',[0,0.25,0.5,0.75,1],'TickLabels',{'9','8','7','6','5'})
% axis([1 zzz 0 1])
% end
% figure(zzz+1)
% hold on
% colormap(jet(10))
% colorbar('Ticks',[0,0.25,0.5,0.75,1],'TickLabels',{'9','8','7','6','5'})
% axis([1 zzz 0 25000])
% % xla
% 
% 
