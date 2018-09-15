% Computing the Magnitude Law to relate different bodies
clear all
% close all
set(0,'DefaultAxesTitleFontWeight','normal');

dummy_counter=0;
start=100000;
passo = 20000;
fine =1000000;
ballwithoutliers=zeros(length(start:passo:fine),5);
ballnooutliers=zeros(length(start:passo:fine),5);
for kkk=start:passo:fine
    dummy_counter=dummy_counter+1;
    load(['/home/limo/PHD/PokerFlat2014/AzDecImages/',num2str(kkk),'Day31032014.mat'])
%     load('113000Day31032014.mat')
%     load('455000Day31032014.mat')
    %loading data
    Vmag = zeros(length(star),1);Bmag=zeros(length(star),1);Lum=zeros(length(star),1); B_V=zeros(length(star),1);
    AirMass=zeros(length(star),1);
    Groups=[];
    for iii=1:length(star)
        Vmag(iii)=star{iii}.Vmag;
        Bmag(iii)=star{iii}.Bmag;
        Lum(iii)=star{iii}.Luminosity(1);
        Groups=[Groups;star{iii}.Spectrum];
        B_V(iii)=star{iii}.B_V;
        AirMass(iii)=star{iii}.AirMass(1);
        Sigma(iii)=1.0875/star{iii}.SNR(1);
    end
    %(GMM). We report it in fig.\ref{fig:TLS} the use of TLS. As we can see in the top row, the in-sample estimate for TLS is excellent: the error is almost non-existent. Unfortunately, as shown in the bottom row, the out-of-sample estimate performs much worse than a simple multilinear regression. This is caused by the fact that TLS in-sample computation allows you to estimate a correction on your regression variables: $\tilde{Y}\approx(X+X_t)\theta$, where $X$ is the matrix containing the regressors variables as columns, $X_t$ its correction, and $\theta$ is the vector containing the TLS coefficients of our model. The value of $X_t$ is not known for out-of-sample values, which --- in our case --- leads to the large residuals shown in fig.\ref{fig:TLS} sub-figure d. While over-fitting could be the cause of this large out-of-sample residual, a minimum can be found which still provides a statistically significantly worse fit than that proposed by multilinear case.

    % Sorting according to star type
    %
    [GroupsSorted, bbb]=sort(Groups);
    Lum=Lum(bbb);
    GroupsSorted=GroupsSorted(~isnan(Lum));%,'descend');
    Vmag=Vmag(bbb);Vmag=Vmag(~isnan(Lum));
    Bmag=Bmag(bbb);Bmag=Bmag(~isnan(Lum));  
    B_V=B_V(bbb);B_V=B_V(~isnan(Lum));  
    AirMass=AirMass(bbb);AirMass=AirMass(~isnan(Lum));
    a=1:3:length(Vmag);
    b=1:length(Vmag);
    c=ismember(b,a);
    d=b(c==0);
    AirMassInSample=AirMass(d);
    AirMassOutOfSample=AirMass(a);
    Lum=Lum(~isnan(Lum));
    LumPrimed=-100^(1/5)*log10(Lum);    
    IndependentVar=LumPrimed(d); %(aka x, what we regress on
    OutOfSampleIndep=LumPrimed(a);
    DependentVar=Vmag(d); % aVa y, our target
    OutOfSamplDep=Vmag(a);
%     ResidualsX=Vmag;
    scrittay1=['V-mag'];
    scrittax=['Instrumental mag'];
    scrittatitle1=['m as function of V-mag'];
    scrittax2=scrittay1;
    scrittay2=['Residual'];
    scrittatitle2=['Residual of V-mag vs m'];
    color_scritta=['WithColor'];
    outliers{1}='WithOutliers';
    outliers{2}='NoOutliers';
    for zzz=1:2
        % Equation:
        % InstrMag = c_0 + c1*V_mag+c2*AirMass+c3*CI+c4*CI*Airmass
        % CI = B-B Magntiude
        %
        if zzz == 2
            contain0 = (rint(:,1)<0 & rint(:,2)>0);
            idx = find(contain0==false);
%             figure(1)
%             gscatter(IndependentVar,DependentVar,GroupsSorted,'bcygrmk')
%             hold on
%             scatter(IndependentVar(idx),DependentVar(idx),'rx')
%             print([scrittax,'vs',scrittay1,color_scritta,outliers{zzz-1},'.eps'],'-depsc')
            DependentVar(idx)=NaN;           
        end
        X=[ones(length(IndependentVar),1),IndependentVar,AirMassInSample];%,B_V,B_V.*AirMass];
        [b,bint,r,rint,stats] = regress(DependentVar,X,0.05);
        stats
%         gca=figure()
        %plotting the data
%         gscatter(IndependentVar,DependentVar,GroupsSorted,'bcygrmk')
%         hold on
%         scatter(IndependentVar,sum(b.*X'),'bo')
        % Using TLS
%         [Xnew,Xprime,Y,coeffs]=TotalLeasSquare(X(1:10:end,:),DependentVar(1:10:end));
%         scatter(X(:,2),sum(coeffs.*X'),'r+')
%         xlabel(scrittax)
%         ylabel(scrittay1)
%         title(scrittatitle1)
%         print([scrittax,'vs',scrittay1,color_scritta,outliers{zzz},'.eps'],'-depsc')
        %plotting the residuals
%         figure()
%         gscatter(ResidualsX,r,GroupsSorted,'bcygrmk')
%         xlabel(scrittax2)
%         ylabel(scrittay2)
%         title(scrittatitle2)
%         print(['Residuals',scrittax,'vs',scrittay1,color_scritta,outliers{zzz},'.eps'],'-depsc')
        if zzz==1
            ballwithoutliers(dummy_counter,1:size(X,2))=b';
            r2_stats_with_outliers(dummy_counter)=stats(1);
            error_variance_with_outliers(dummy_counter)=stats(4);
            oos_error_var_with_out(dummy_counter)=nanvar(sum(b'.*[ones(length(OutOfSampleIndep),1),OutOfSampleIndep,AirMassOutOfSample],2)-OutOfSamplDep)
            figure(1)
            plot(OutOfSampleIndep,sum(b'.*[ones(length(OutOfSampleIndep),1),OutOfSampleIndep,AirMassOutOfSample],2),'o');
            hold on
            plot(OutOfSampleIndep,OutOfSamplDep,'.');
        else
            ballnooutliers(dummy_counter,1:size(X,2))=b';
            r2_stats_no_outliers(dummy_counter)=stats(1);
            error_variance_no_outliers(dummy_counter)=stats(4);
            oos_error_var_no_out(dummy_counter)=nanvar(sum(b'.*[ones(length(OutOfSampleIndep),1),OutOfSampleIndep,AirMassOutOfSample],2)-OutOfSamplDep)
            figure(2)
            plot(OutOfSampleIndep,sum(b'.*[ones(length(OutOfSampleIndep),1),OutOfSampleIndep,AirMassOutOfSample],2),'o');
            hold on
            plot(OutOfSampleIndep,OutOfSamplDep,'.');
        end
%     saveas(
    end

    pause
        close all
end


function [Xnew,Xprime,Y,coeffs]=TotalLeasSquare(X,Y)
X=X(~isnan(Y),:);
Xprime=X;
Y=Y(~isnan(Y));
sum(isnan(Y));
sum(isnan(X));
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