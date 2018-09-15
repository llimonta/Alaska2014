import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from scipy import interpolate
from scipy.integrate import trapz
from scipy.optimize import curve_fit
import itertools


#dictionary for plotting
color_line = {'Iron':'k','Normal':'orange', 'Na-poor':'g','Na-enhanced':'y','Na-rich':'b', 'Na-free':'r','Fe-poor':'cyan'}
#dictionary for creating mean spectral line
meteor_types_counter = {'Iron':0,'Normal':0, 'Na-poor':0,'Na-enhanced':0,'Na-rich':0, 'Na-free':0,'Fe-poor':0,'Other':0} # NB this is an updatable dictionary used later on in the code

def main():
    text_file_sporadic = open('/home/limo/PHD/PokerFlat2014/MeanSpectrum/MeteorShowerVsSporadic.txt', 'r')
    text_file_characteristic = open('/home/limo/PHD/PokerFlat2014/MeanSpectrum/MeteorSpecifications.txt', 'r')
    # obtain response spectra by Camera
    qeBVdata = pd.read_csv('/home/limo/Documents/Thesis/CodeForPlots/CameraQEBvBvf.csv')
    qeBVx_points = qeBVdata['x']*10
    qeBVy_points = qeBVdata['Curve2']
    params1,cov1=curve_fit(QECurveFit,qeBVx_points,qeBVy_points,)
    # print params1
    # plt.plot(qeBVx_points,qeBVy_points)
    # plt.plot(np.linspace(min(qeBVx_points),max(qeBVx_points),1000),QECurveFit(np.linspace(min(qeBVx_points),max(qeBVx_points),1000),*params1))
    # plt.plot(qeBVx_points,QECurveFit(qeBVx_points,*params1))
    # plt.show()
    computing_integrals_all_spectra(text_file_sporadic,text_file_characteristic,params1)

    return 0


def computing_integrals_all_spectra(text_file_sporadic,text_file_characteristic,params1):
    # put file at beginning
    text_file_sporadic.seek(0)
    text_file_characteristic.seek(0)
    x_text_file_sporadic = text_file_sporadic.readlines()
    x_text_file_characteristic = text_file_characteristic.readlines()
    # finding the minimum and maximum wavelengths on which to interpolate for mean_spectra
    meteor_min_spectrum = np.zeros(1)
    meteor_max_spectrum = np.zeros(1)
    all_meteor_integral_ratios=np.zeros(1)
    color_meteor=[]
    meteor_types_legend=[]
    dummy_counter=0
    fig = plt.figure(figsize=(8, 6))
    for iii in range(0,len(x_text_file_sporadic)):
        is_sporadic = x_text_file_sporadic[iii][:].split('|')[-1].rstrip() #rstrip removes the \n at the end, spaces et al
        if is_sporadic == str('SPO'):
            # print iii
            dummy_counter+=1
            name_meteor = x_text_file_sporadic[iii][:].split('|')[0].rstrip()
            spo_characteristic = x_text_file_characteristic[iii][:].split('|')[-1].rstrip()
            meteor_data = pd.read_csv("".join(['/home/limo/PHD/PokerFlat2014/MeanSpectrum/',name_meteor,'.dat']),delim_whitespace = True, skiprows =[-1], header = None)
            meteor_spectrum=np.float64(meteor_data[0][:])
            meteor_spectral_response=np.float64(meteor_data[2][:])
            # meteor_spectral_response[meteor_spectral_response<0]=np.NaN
            # meteor_spectral_response[meteor_spectral_response<0]=np.nan
            # print meteor_spectrum,meteor_spectral_response
            if spo_characteristic in color_line:
                color_meteor.append(color_line[spo_characteristic])
                meteor_types_counter[spo_characteristic]=meteor_types_counter[spo_characteristic]+1
                meteor_types_legend.append(spo_characteristic)
            else:
                color_meteor.append('magenta')
                meteor_types_counter['Other']=meteor_types_counter['Other']+1
                meteor_types_legend.append('Other')
            # Compute spectral response seen by camera
            Meteor_Camera_Response=meteor_spectral_response*QECurveFit(meteor_spectrum[:],*params1)/100
            # Compute ingteral
            MeteorIntegral = trapz(meteor_spectral_response,meteor_spectrum)
            CameraMeteorIntegral = trapz(Meteor_Camera_Response,meteor_spectrum)
            all_meteor_integral_ratios= np.append(all_meteor_integral_ratios,CameraMeteorIntegral/MeteorIntegral)
            plt.plot(dummy_counter,CameraMeteorIntegral/MeteorIntegral,color_meteor[dummy_counter-1],marker='.',markersize=10)
    custom_lines = [Line2D([0], [0], color=color_line['Normal'],marker='.',linestyle='None'),
                    Line2D([0], [0], color=color_line['Iron'],marker='.',linestyle='None'),
                    Line2D([0], [0], color=color_line['Fe-poor'],marker='.',linestyle='None'),
                    Line2D([0], [0], color=color_line['Na-rich'], marker='.',linestyle='None'),
                    Line2D([0], [0], color=color_line['Na-enhanced'],marker='.',linestyle='None'),
                    Line2D([0], [0], color=color_line['Na-poor'], marker='.',linestyle='None'),
                    Line2D([0], [0], color=color_line['Na-free'], marker='.',linestyle='None'),
                    Line2D([0], [0], color='magenta',marker='.',linestyle='None')]
    plt.xlabel('Meteor #',fontsize=20)
    plt.ylabel('Integral Ratio',fontsize=24)
    plt.title('Ratio between Camera Response and Ideal Response',fontsize=24)
    plt.legend(custom_lines, ['Normal','Fe','Fe-poor','Na-rich','Na-enhanced','Na-poor','Na-free', 'other'],fontsize=15)
    plt.show()
    print all_meteor_integral_ratios
    print np.mean(all_meteor_integral_ratios[1::]),np.std(all_meteor_integral_ratios[1::])
    return 0 #min(meteor_min_spectrum[1::]),max(meteor_max_spectrum[1::])


def QECurveFit(x,a,b,c,d,e,f,g,h,n,l,m):
    # print type(x)
    # print a,b,c,d,e,f,g,h,n,l,m
    # print a*x,b*x**2,c*x**3,d*x**4,e*x**5,f*x**6,g*x**7,h*x**8,n*x**9,l*x**10,m*x**11
    return  a*x+b*x**2+c*x**3+d*x**4+e*x**5+f*x**6+g*x**7+h*x**8+n*x**9+l*x**10+m*x**11



if __name__ == '__main__':
    mat=main()
