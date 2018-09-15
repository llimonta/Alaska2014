import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate
from scipy.integrate import trapz
from scipy.optimize import curve_fit
import itertools


#dictionary for plotting
color_line = {'Iron':'k','Normal':'orange', 'Na-poor':'.g','Na-enhanced':'.y','Na-rich':'.b', 'Na-free':'.r','Fe-poor':'cyan'}
#dictionary for creating mean spectral line
meteor_types_counter = {'Iron':0,'Normal':0, 'Na-poor':0,'Na-enhanced':0,'Na-rich':0, 'Na-free':0,'Fe-poor':0,'Other':0} # NB this is an updatable dictionary used later on in the code

def main():
    text_file_sporadic = open('/home/limo/PHD/PokerFlat2014/MeanSpectrum/MeteorShowerVsSporadic.txt', 'r')
    text_file_characteristic = open('/home/limo/PHD/PokerFlat2014/MeanSpectrum/MeteorSpecifications.txt', 'r')
    print 'Plotting all spectra'
    min_wave,max_wave = plotting_all_spectra(text_file_sporadic,text_file_characteristic)
    wavelengths_2_interp = np.linspace(min_wave,max_wave,1000)
    # obtain expected spectra with associated variance
    print 'Mean spectrum'
    mean_spectra , mean_spectra_std, mean_spectra_cov =plotting_mean_spectra(text_file_sporadic,text_file_characteristic,wavelengths_2_interp)
    # obtain response spectra by Camera
    qeBVdata = pd.read_csv('/home/limo/Documents/Thesis/CodeForPlots/CameraQEBvBvf.csv')
    qeBVx_points = qeBVdata['x']*10
    qeBVy_points = qeBVdata['Curve2']
    params1,cov1=curve_fit(QECurveFit,qeBVx_points,qeBVy_points,)
    camera_QE =QECurveFit(wavelengths_2_interp,*params1)/100
    #integral of mean spectra
    # Deterministic method
    print 'Computing integral and associated error'
    integral_mean_spectra = trapz(mean_spectra,wavelengths_2_interp)
    delta_x_integral = wavelengths_2_interp[1]-wavelengths_2_interp[0]
    deltaL = delta_x_integral
    integral_mean_spectra_std = np.sqrt((deltaL/2)**2*(mean_spectra_std[0]**2+mean_spectra_std[-1]**2)+deltaL**2*(sum(np.square(mean_spectra_std))))
    minimum_integral_mean_spectra = trapz(mean_spectra-mean_spectra_std,wavelengths_2_interp)
    maximum_integral_mean_spectra = trapz(mean_spectra+mean_spectra_std,wavelengths_2_interp)
    # adding co-variance, very time consuming
    cov_int = np.sqrt(integral_mean_spectra_std**2 + deltaL**2*sum([mean_spectra_cov[i,j] for i,j in itertools.product(range(1,len(mean_spectra_cov)-1), range(1,len(mean_spectra_cov)-1)) if i!=j]) + deltaL**2/2*sum([mean_spectra_cov[i,j] for i,j in itertools.product(range(1,2), range(1,len(mean_spectra_cov))) if i!=j]) + deltaL**2/2*sum([mean_spectra_cov[i,j] for i,j in itertools.product(range(0,len(mean_spectra_cov)), range(len(mean_spectra_cov)-1,len(mean_spectra_cov))) if i!=j]))

    # compute variance of integral_mean_spectra
    # MonteCarlo method
    print 'Montecarlo'
    nsamples = 1000
    integral_mean_spectra_mc,integral_mean_spectra_std_mc, integral_mean_spectra_mc_mvn, integral_mean_spectra_std_mc_mvn = MonteCarloIntegration(mean_spectra,mean_spectra_std,mean_spectra_cov,wavelengths_2_interp,nsamples)
    # Integral of the camera response
    integral_camera = trapz(camera_QE*mean_spectra,wavelengths_2_interp)
    integral_camera_spectra_mc,integral_camera_spectra_std_mc, integral_camera_spectra_mc_mvn, integral_camera_spectra_std_mc_mvn = MonteCarloIntegrationCamera(mean_spectra,mean_spectra_std,mean_spectra_cov,camera_QE,wavelengths_2_interp,nsamples)


    print '\nExpected True Integral: ',integral_mean_spectra,'Associated Std:',integral_mean_spectra_std/1226.3706307,'If cov present:',cov_int/1226.3706307# normalizing by exact integral
    print 'Mean Integral From MonteCarlo:', integral_mean_spectra_mc,'Associated Std:',integral_mean_spectra_std_mc
    print 'Mean Integra From MC MVN:',integral_mean_spectra_mc_mvn,'Associated Std:', integral_mean_spectra_std_mc_mvn
    print 'Minimum Integral:', minimum_integral_mean_spectra,'Maximum Integral:',maximum_integral_mean_spectra
    print 'Expected Camera Integral:',integral_camera
    print 'Expected Camera MC:',integral_camera_spectra_mc,'Associated Std:',integral_camera_spectra_std_mc
    print 'Expected Camera MC MVN:',integral_camera_spectra_mc_mvn,'Associated Std:',integral_camera_spectra_std_mc_mvn
    print 'Expected Percentage Seen:',integral_camera/integral_mean_spectra*100 #,' Associated Std:', integral_camera/integral_mean_spectra*integral_mean_spectra_std
    plt.show()
    return 0


def plotting_all_spectra(text_file_sporadic,text_file_characteristic):
    # put file at beginning
    text_file_sporadic.seek(0)
    text_file_characteristic.seek(0)
    x_text_file_sporadic = text_file_sporadic.readlines()
    x_text_file_characteristic = text_file_characteristic.readlines()
    # finding the minimum and maximum wavelengths on which to interpolate for mean_spectra
    meteor_min_spectrum = np.zeros(1)
    meteor_max_spectrum = np.zeros(1)
    for iii in range(0,len(x_text_file_sporadic)):
        is_sporadic = x_text_file_sporadic[iii][:].split('|')[-1].rstrip() #rstrip removes the \n at the end, spaces et al
        if is_sporadic == str('SPO'):
            # print iii
            name_meteor = x_text_file_sporadic[iii][:].split('|')[0].rstrip()
            meteor_spectrum = pd.read_csv("".join(['/home/limo/PHD/PokerFlat2014/MeanSpectrum/',name_meteor,'.dat']),delim_whitespace = True, skiprows =[-1], header = None) # 0 are the lines at which spectrum is # 2 are where the calibrated spectral response is
            meteor_min_spectrum = np.append(meteor_min_spectrum,min(meteor_spectrum[0][:]))
            meteor_max_spectrum = np.append(meteor_max_spectrum,max(meteor_spectrum[0][:]))
            spo_characteristic = x_text_file_characteristic[iii][:].split('|')[-1].rstrip()
            # #plt.figure()
            if spo_characteristic in color_line:
                color = color_line[spo_characteristic]
                meteor_types_counter[spo_characteristic]=meteor_types_counter[spo_characteristic]+1
            else:
                color = 'magenta'
                meteor_types_counter['Other']=meteor_types_counter['Other']+1
            plt.figure(1)
            plt.plot(meteor_spectrum[0][:],meteor_spectrum[2][:]/max(meteor_spectrum[2][:]),color,label=spo_characteristic)
            plt.xlabel('Wavelength $\AA$')
            plt.ylabel('Calibrated Raw Units')
            plt.title('Distribution of expected spectral lines')
            if spo_characteristic == str('Normal'):
                plt.figure(10)
                plt.plot(meteor_spectrum[0][:],meteor_spectrum[2][:],color,label=spo_characteristic)
                plt.xlabel('Wavelength $\AA$')
                plt.ylabel('Calibrated Raw Units')
                plt.title('Distribution of expected spectral lines')

    return min(meteor_min_spectrum[1::]),max(meteor_max_spectrum[1::])

def plotting_mean_spectra(text_file_sporadic,text_file_characteristic,wavelengths_2_interp):
    # put file at beginning
    text_file_sporadic.seek(0)
    text_file_characteristic.seek(0)
    x_text_file_sporadic = text_file_sporadic.readlines()
    x_text_file_characteristic = text_file_characteristic.readlines()
    # total number of meteors
    total_meteors = 0
    dummy_counter = 0
    for key in meteor_types_counter:
        total_meteors = total_meteors + meteor_types_counter[key]
    # print total_meteors
    mean_spectra = np.zeros((total_meteors,len(wavelengths_2_interp)))
    for iii in range(0,len(x_text_file_sporadic)):
        is_sporadic = x_text_file_sporadic[iii][:].split('|')[-1].rstrip() #rstrip removes the \n at the end, spaces et al
        if is_sporadic == str('SPO'):
            name_meteor = x_text_file_sporadic[iii][:].split('|')[0].rstrip()
            meteor_spectrum = pd.read_csv("".join(['/home/limo/PHD/PokerFlat2014/MeanSpectrum/',name_meteor,'.dat']),delim_whitespace = True, skiprows =[-1], header = None) # 0 are the lines at which spectrum is # 2 are where the calibrated spectral response is
            meteor_min_spectrum = min(meteor_spectrum[0][:])
            meteor_max_spectrum = max(meteor_spectrum[0][:])

            # Find closest point to min_spectrum, find closest point to max spectrum
            starting_index_wave_2_interp = [iii for iii in range(0,len(wavelengths_2_interp)) if wavelengths_2_interp[iii]- meteor_min_spectrum>=0][0]
            ending_index_wave_2_interp = [iii for iii in range(0,len(wavelengths_2_interp)) if wavelengths_2_interp[iii]- meteor_max_spectrum<=0][-1]
            # Interpolate Over those Points, call common_points_interpolation
            wavelengths_2_interp_common_elements = wavelengths_2_interp[starting_index_wave_2_interp:ending_index_wave_2_interp+1]
            meteor_spectrum_interpolated = common_points_interpolation(wavelengths_2_interp_common_elements,meteor_spectrum[0][:], meteor_spectrum[2][:]/max(meteor_spectrum[2][:]))
            # NaN padding to reach length of wavelengths_2_interp
            meteor_spectrum_interpolated_padded = np.append(np.append(np.zeros(starting_index_wave_2_interp)+np.nan, meteor_spectrum_interpolated), np.zeros(len(wavelengths_2_interp)-ending_index_wave_2_interp-1)+np.nan)
            # take out negative values because unphysical
            meteor_spectrum_interpolated_padded[meteor_spectrum_interpolated_padded<0]=np.nan
            # #plt.figure(2)
            # #plt.plot(wavelengths_2_interp,meteor_spectrum_interpolated_padded,meteor_spectrum[0][:],meteor_spectrum[2][:]/max(meteor_spectrum[2][:]),linewidth = 1)
            mean_spectra[dummy_counter][:] = meteor_spectrum_interpolated_padded
            dummy_counter += 1
    mean_spectra_std = np.nanstd(mean_spectra, axis = 0)
    dummy_mean_spectra = np.nan_to_num(mean_spectra)
    mean_spectra_cov = np.cov(np.transpose(dummy_mean_spectra))
    mean_spectra = np.nanmean(mean_spectra, axis=0)
    fig, ax = plt.subplots()
    ax.plot(wavelengths_2_interp,mean_spectra,'g',linewidth = 5,label='Mean Spectrum')
    ax.fill_between(wavelengths_2_interp, mean_spectra+mean_spectra_std, mean_spectra-mean_spectra_std, color=[1.0, 0.5, 0.5], alpha=0.4, label = 'Standard Deviation ')
    plt.xlabel('Wavelength, $\AA$')
    plt.ylabel('Arbitrary Units')
    plt.title('Mean Spectrum Emission')
    plt.legend(loc='upper right')
    return mean_spectra,mean_spectra_std, mean_spectra_cov

    # Weight all points accordingly to their expected appearance
    # Plot mean_spectra

def common_points_interpolation(wavelengths_2_interp_common_elements,meteor_wavelengths, meteor_spectral_response):
    f = interpolate.interp1d(meteor_wavelengths, meteor_spectral_response,kind = 'cubic')
    return f(wavelengths_2_interp_common_elements)

def QECurveFit(x,a,b,c,d,e,f,g,h,i,l,m):
    return a+b*x**2+c*x**3+d*x**4+e*x**5+f*x**6+g*x**7+h*x**8+i*x**9+l*x**10+m*x**11

def  MonteCarloIntegration(mean_spectra,mean_spectra_std,mean_spectra_cov,wavelengths_2_interp,nsamples):
     # using uncorrelated values
     storage_of_spectra = np.array([np.random.normal(mean_spectra[iii],mean_spectra_std[iii],nsamples) for iii in range(0,len(mean_spectra))])
     storage_of_integral = np.array([trapz(storage_of_spectra[0:len(storage_of_spectra),iii],wavelengths_2_interp) for iii in range(0,nsamples) ])/1226.3706307# normalizing by exact integral
     # using correlation matrix
     storage_of_spectra_mvn = np.random.multivariate_normal(mean_spectra,mean_spectra_cov,nsamples)
     storage_of_integral_mvn = np.array([trapz(storage_of_spectra_mvn[iii,0:],wavelengths_2_interp) for iii in range(0,nsamples)])/1226.3706307# normalizing by exact integral
     return np.nanmean(storage_of_integral),np.nanstd(storage_of_integral),np.nanmean(storage_of_integral_mvn),np.nanstd(storage_of_integral_mvn)

def  MonteCarloIntegrationCamera(mean_spectra,mean_spectra_std,mean_spectra_cov,camera_QE,wavelengths_2_interp,nsamples):
     # using uncorrelated values
     storage_of_spectra = np.array([np.random.normal(mean_spectra[iii],mean_spectra_std[iii],nsamples) for iii in range(0,len(mean_spectra))])
     storage_of_integral = np.array([trapz(storage_of_spectra[0:len(storage_of_spectra),iii]*camera_QE,wavelengths_2_interp) for iii in range(0,nsamples) ])/1226.3706307# normalizing by exact integral
     # using correlation matrix
     storage_of_spectra_mvn = np.random.multivariate_normal(mean_spectra,mean_spectra_cov,nsamples)
     storage_of_integral_mvn = np.array([trapz(storage_of_spectra_mvn[iii,0:]*camera_QE,wavelengths_2_interp) for iii in range(0,nsamples)])/1226.3706307# normalizing by exact integral
     # For plotting!!!, yes the mean_spectra_std_uncorr needs to stay there for plotting
     mean_spectra_std_uncorr = np.nanstd(storage_of_spectra, axis = 1)
     fig, ax = plt.subplots()
     ax.plot(wavelengths_2_interp,mean_spectra*camera_QE,'g',linewidth = 5,label='Mean Spectrum')
     ax.fill_between(wavelengths_2_interp, mean_spectra*camera_QE+mean_spectra_std_uncorr, mean_spectra*camera_QE-mean_spectra_std_uncorr, color=[1.0, 0.5, 0.5], alpha=0.4, label = 'Standard Deviation ')
     plt.xlabel('Wavelength, $\AA$')
     plt.ylabel('Arbitrary Units')
     plt.title('Mean Spectrum Response Camera')
     plt.legend(loc='upper right')
     plt.show()
     return np.nanmean(storage_of_integral),np.nanstd(storage_of_integral),np.nanmean(storage_of_integral_mvn),np.nanstd(storage_of_integral_mvn)



if __name__ == '__main__':
    mat=main()
