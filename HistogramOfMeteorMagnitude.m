clear all
close all

set(0,'DefaultAxesTitleFontWeight','normal');


Fparameter=[];
%load meteor magnitude for 31st
dummy_counter=1;
% Fparameter events
%[1,2,6,12,13,17,21,23,24,34,35,36,38,41,43,46,47,50,54,55,60,61,67,68,70,72,76]
[1,3,24,32,46,50,52,60,65,66,72]

% [1,24,46,50,60,72]c;
for iii=[1:76]
    dummy_counter    
    iii
    figure(1)
    load(['/home/limo/PHD/PokerFlat2014/AzDecImages/VmagMeteor',num2str(iii),'Day31032014.mat'])
%     dummy_variable=[1:length(VmagMeteor)]';
%     f = fit(dummy_variable(~isnan(VmagMeteor)),VmagMeteor(~isnan(VmagMeteor)),'smoothingspline','SmoothingParam',0.1);
%     luminosity_smooth_spline = feval(f,dummy_variable);
%     [~,idx_min]=min(luminosity_smooth_spline);
%     Fparameter = [Fparameter;(idx_min-2)/(length(luminosity_smooth_spline)-4)];
    VmagMin(dummy_counter)=min(VmagMeteor);
    VmagAltMin(dummy_counter)=min(VmagAlt);
    dummy_counter=dummy_counter+1;
%     plot(VmagMeteor)
%     set(gca,'ydir','reverse');
%     hold on
%     plot(luminosity_smooth_spline)
%     set(gca,'ydir','reverse');
%     pause
    close all
end
%
for iii=[4,5,7,9,11,18,23,26,28]
    load(['/home/limo/PHD/PokerFlat2014/AzDecImages/VmagMeteor',num2str(iii),'Day30032014.mat'])
    VmagMin(dummy_counter)=min(VmagMeteor);
    VmagAltMin(dummy_counter)=min(VmagAlt);
%     dummy_counter=dummy_counter+1;  
%     dummy_variable=[1:length(VmagMeteor)]';
%     f = fit(dummy_variable(~isnan(VmagMeteor)),VmagMeteor(~isnan(VmagMeteor)),'smoothingspline','SmoothingParam',0.1);
%     luminosity_smooth_spline = feval(f,dummy_variable);
%     [~,idx_min]=min(luminosity_smooth_spline);
%     Fparameter = [Fparameter;(idx_min-2)/(length(luminosity_smooth_spline)-4)];
    VmagMin(dummy_counter)=min(VmagMeteor);
    VmagAltMin(dummy_counter)=min(VmagAlt);
    dummy_counter=dummy_counter+1;
%     plot(VmagMeteor)
%     set(gca,'ydir','reverse');
%     hold on
%     plot(luminosity_smooth_spline)
%     set(gca,'ydir','reverse');
%     pause
    close all
end 

figure,hist(VmagMin)
xlabel('Meteor Minimum Magnitude')
ylabel('#')
title('Meteor magnitude distribution')

figure,hist(VmagAltMin)
xlabel('Meteor Minimum Magnitude')
ylabel('#')
title('Meteor magnitude distribution - indirect method')

% 30/03:[4,5,7,9,11,17,23,26,28]
% 31/03:[1,3,24,32,40,46,50,52,60,65,66,72]

dummy_counter=1;

Fparameter2=[];
for iii=[1,3,24,32,40,46,50,52,60,65,66,72]
    load(['/home/limo/PHD/PokerFlat2014/AzDecImages/VmagMeteor',num2str(iii),'Day31032014.mat'])
    DualVmagMin(dummy_counter)=min(VmagMeteor);
    DualVmagAltMin(dummy_counter)=min(VmagAlt);
    dummy_counter=dummy_counter+1;
    
    dummy_variable=[1:length(VmagMeteor)]';
    f = fit(dummy_variable(~isnan(VmagMeteor)),VmagMeteor(~isnan(VmagMeteor)),'smoothingspline','SmoothingParam',0.1);
    luminosity_smooth_spline = feval(f,dummy_variable);
    [~,idx_min]=min(luminosity_smooth_spline);
    Fparameter2 = [Fparameter2;(idx_min-2)/(length(luminosity_smooth_spline)-4)];
end

%[4,5,7,9,11,18,23,26,28]
for iii=[4,5,7,9,11,18,23,26,28]
    load(['/home/limo/PHD/PokerFlat2014/AzDecImages/VmagMeteor',num2str(iii),'Day30032014.mat'])
    DualVmagMin(dummy_counter)=min(VmagMeteor);
    DualVmagAltMin(dummy_counter)=min(VmagAlt);
    dummy_counter=dummy_counter+1;
    dummy_variable=[1:length(VmagMeteor)]';
    f = fit(dummy_variable(~isnan(VmagMeteor)),VmagMeteor(~isnan(VmagMeteor)),'smoothingspline','SmoothingParam',0.1);
    luminosity_smooth_spline = feval(f,dummy_variable);
    [~,idx_min]=min(luminosity_smooth_spline);
    Fparameter2 = [Fparameter2;(idx_min-2)/(length(luminosity_smooth_spline)-4)];
end 

figure,hist(DualVmagMin)
xlabel('Meteor Minimum Magnitude')
ylabel('#')
title('Dual observations meteor magnitude distribution')

figure,hist(DualVmagAltMin)
xlabel('Meteor Minimum Magnitude')
ylabel('#')
title('Dual observations meteor magnitude distribution - indirect method')