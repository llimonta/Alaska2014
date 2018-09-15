import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from scipy import interpolate
from scipy.integrate import trapz
from scipy.optimize import curve_fit
from scipy.stats.stats import pearsonr
from scipy import stats
import itertools


#dictionary for plotting
color_line = {'Iron':'k','Normal':'orange', 'Na-poor':'g','Na-enhanced':'y','Na-rich':'b', 'Na-free':'r','Fe-poor':'cyan'}
#dictionary for creating mean spectral line
meteor_types_counter = {'Iron':0,'Normal':0, 'Na-poor':0,'Na-enhanced':0,'Na-rich':0, 'Na-free':0,'Fe-poor':0,'(Atmospheric lines)':0} # NB this is an updatable dictionary used later on in the code

meteors_camera_response_by_type = pd.DataFrame(np.NaN*np.random.randint(low=0, high=10, size=(50, 8)), columns=['Iron','Normal', 'Na-poor','Na-enhanced','Na-rich', 'Na-free','Fe-poor','(Atmospheric lines)'])
def main():
    text_file_sporadic = '/home/limo/PHD/PokerFlat2014/MeanSpectrum/MeteorShowerVsSporadic.txt'
    text_file_characteristic ='/home/limo/PHD/PokerFlat2014/MeanSpectrum/MeteorSpecifications.txt'
    # obtain response spectra by Camera
    # obtain response spectra by Camera
    qeBVdata = pd.read_csv('/home/limo/Documents/Thesis/CodeForPlots/CameraQEBvBvf.csv')
    qeBVx_points = qeBVdata['x']*10
    qeBVy_points = qeBVdata['Curve2']
    params1,cov1=curve_fit(QECurveFit,qeBVx_points,qeBVy_points,)
    plot_meteor_spectra_velocity(text_file_sporadic,text_file_characteristic,params1)

    return 0


# def computing_integrals_all_spectra(text_file_sporadic,text_file_characteristic,params1):
def plot_meteor_spectra_velocity(text_file_sporadic,text_file_characteristic,params1):
    # put file at beginning
    color_meteor=[]
    meteor_types_legend=[]
    dummy_counter=0
    fig = plt.figure(figsize=(8, 6))
    all_sporadic=pd.read_csv(text_file_sporadic,delimiter='|',header=0)
    meteor_char=pd.read_csv(text_file_characteristic,delimiter='|',header=0)
    dummy_counter=0
    all_sporadic_velocities=[]
    all_sporadic_intensity_ratios=[]
    meteor_list_type=[]
    meteor_list_integral_ratio=[]
    print
    # print meteor_char
    for iii in range(0,len(all_sporadic)):
         if all_sporadic['Shower'][iii] in 'SPO':
             dummy_counter+=1
             spo_characteristic = meteor_char['Class'][iii]
             # print meteor_spectrum,meteor_spectral_response
             if spo_characteristic in color_line:
                 print spo_characteristic
                 color_meteor.append(color_line[spo_characteristic])
                 meteor_types_counter[spo_characteristic]=meteor_types_counter[spo_characteristic]+1
                 meteor_types_legend.append(spo_characteristic)
             else:
                 print spo_characteristic
                 color_meteor.append('magenta')
                 meteor_types_counter['(Atmospheric lines)']=meteor_types_counter['(Atmospheric lines)']+1
                 meteor_types_legend.append('(Atmospheric lines)')
             # print iii
             # print meteor_char['elinf'][iii]
             # print color_meteor
             print meteor_types_counter[spo_characteristic]
             name_meteor = all_sporadic['SpFile '][iii].rstrip()
             meteor_data = pd.read_csv("".join(['/home/limo/PHD/PokerFlat2014/MeanSpectrum/',name_meteor,'.dat']),delim_whitespace = True, skiprows =[-1], header = None)
             meteor_spectrum=np.float64(meteor_data[0][:])
             meteor_spectral_response=np.float64(meteor_data[2][:])
             # Compute spectral response seen by camera
             Meteor_Camera_Response=meteor_spectral_response*QECurveFit(meteor_spectrum[:],*params1)/100
             # Compute ingteral
             MeteorIntegral = trapz(meteor_spectral_response,meteor_spectrum)
             CameraMeteorIntegral = trapz(Meteor_Camera_Response,meteor_spectrum)
             Meteor_Integral_Ratio = CameraMeteorIntegral/MeteorIntegral
             # print spo_characteristic,type(meteors_camera_response_by_type)
             meteors_camera_response_by_type.iloc[meteor_types_counter[spo_characteristic]-1,[list(meteors_camera_response_by_type.columns.values).index(spo_characteristic)][0]]=CameraMeteorIntegral/MeteorIntegral
             meteor_list_type.append(spo_characteristic)
             all_sporadic_intensity_ratios.append(Meteor_Integral_Ratio)
             all_sporadic_velocities.append(meteor_char['elinf'][iii])
             plt.plot(meteor_char['elinf'][iii],CameraMeteorIntegral/MeteorIntegral,color_meteor[dummy_counter-1],marker='.',markersize=10)
    custom_lines = [Line2D([0], [0], color=color_line['Normal'],marker='.',linestyle='None'),
                    Line2D([0], [0], color=color_line['Iron'],marker='.',linestyle='None'),
                    Line2D([0], [0], color=color_line['Fe-poor'],marker='.',linestyle='None'),
                    Line2D([0], [0], color=color_line['Na-rich'], marker='.',linestyle='None'),
                    Line2D([0], [0], color=color_line['Na-enhanced'],marker='.',linestyle='None'),
                    Line2D([0], [0], color=color_line['Na-poor'], marker='.',linestyle='None'),
                    Line2D([0], [0], color=color_line['Na-free'], marker='.',linestyle='None'),
                    Line2D([0], [0], color='magenta',marker='.',linestyle='None')]
    plt.xlabel('Meteor Velocity [km/s]',fontsize=20)
    plt.ylabel('Integral Ratio',fontsize=20)
    plt.title('Camera Response vs. Meteor Velocity ',fontsize=20)
    plt.tight_layout
    # plt.tight_layout()
    # plt.legend(custom_lines, ['Normal','Fe','Fe-poor','Na-rich','Na-enhanced','Na-poor','Na-free', '(Atmospheric lines)'],fontsize=15)
    print 'Correlation coefficient', pearsonr(all_sporadic_velocities, all_sporadic_intensity_ratios)[0]**2
    x = np.array(all_sporadic_velocities)
    y = np.array(all_sporadic_intensity_ratios)
    # z = np.polyfit(x, y, 1)
    # p = np.poly1d(z)
    # xp=np.linspace(min(x),max(x),100)
    # plt.plot(xp,p(xp),'-y')
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    print slope, intercept, r_value, p_value, std_err
    plt.figure()
    boxplot = meteors_camera_response_by_type.boxplot()

    print meteor_list_type, all_sporadic_intensity_ratios
    # plt.figure()
    # plt.plot(x,p(x)-y,'.')
    # plt.title('Residuals')
    plt.show()
    return 0 #min(meteor_min_spectrum[1::]),max(meteor_max_spectrum[1::])


def QECurveFit(x,a,b,c,d,e,f,g,h,n,l,m):
    # print type(x)
    # print a,b,c,d,e,f,g,h,n,l,m
    # print a*x,b*x**2,c*x**3,d*x**4,e*x**5,f*x**6,g*x**7,h*x**8,n*x**9,l*x**10,m*x**11
    return  a*x+b*x**2+c*x**3+d*x**4+e*x**5+f*x**6+g*x**7+h*x**8+n*x**9+l*x**10+m*x**11

if __name__ == '__main__':
    mat=main()
