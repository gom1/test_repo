# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 17:06:13 2017

@author: MalinDahl
"""

from pegasus import Experiment #Custom package for Tribogenics
import matplotlib.pyplot as plt
from pegasus import spectrum_characterization_tools as sct
import numpy as np
from scipy.optimize import curve_fit
from Tkinter import Tk
from tkFileDialog import askopenfilename, askdirectory
import os
import pandas as pd
from scipy.stats import chisquare
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
from lmfit.models import LorentzianModel


Tk().withdraw()
path = askopenfilename(title='Pick an analysis file created by Mustang Aquire')
save_path = askdirectory(title='Pick folder to save figures in')
#peaks_to_remove = [(8.0, 11.0), (21.5, 22.8), (24.5, 25.5)]
#peaks_to_remove = [(21.5, 22.8), (24.5, 25.5)]
peaks_to_remove = [(8,8.8),(9.2,10.2),(21.5, 22.8), (24.5, 25.5)]

PC = False
Pegasus = False


def get_index_of_energy(spectrum_axis, energy):
    index = (np.abs(np.array(spectrum_axis) - energy)).argmin()
    return index

#Remove the peaks from the spectrum
def get_spectrum_without_peaks(spectrum, axis):
    Start = 106
    End = 1200
    peaks = []
    peaks_axis = []

    for j in range(450,1200):
        Sums = 0
        for u in range(0,10):
            Sums += spectrum[j-u]
        if Sums < 5:
             End = j
             break

    spectrum = spectrum[Start:End]
    axis = axis[Start:End]
    for peak in peaks_to_remove:
        start_index = get_index_of_energy(axis, peak[0])
        end_index = get_index_of_energy(axis, peak[1])

        beginning_spectrum = np.array(spectrum[:start_index])
        end_spectrum = np.array(spectrum[end_index:])
        peak = np.array(spectrum[start_index:end_index])
        
        beginning_axis = np.array(axis[:start_index])
        end_axis = np.array(axis[end_index:])
        peak_axis = np.array(axis[start_index:end_index])

        spectrum = np.append(beginning_spectrum, end_spectrum)
        axis = np.append(beginning_axis, end_axis)

        peak_spectrum = np.array(peak)
        peak_axis_spectrum=np.array(peak_axis)
        
        peaks.append(peak_spectrum)
        peaks_axis.append(peak_axis_spectrum)

    return spectrum, axis, peaks, peaks_axis

#fit the peaks of the spectrum with Lorentzian 
def get_peak_fits(peaks_axis, peaks):
    all_fits = []
    params = []
    for j in range(len(peaks)):
        #print 'Peak number: %i'%j
        #mod = LorentzianModel()

        if j == 0:
            best_vals, covar = curve_fit(Lorentzian, peaks_axis[j], peaks[j], [10,20,10])
            fit = Lorentzian(peaks_axis[j], *best_vals)
        else:
            best_vals, covar = curve_fit(Lorentzian, peaks_axis[j], peaks[j], [10, 20, 10])
            fit = Lorentzian(peaks_axis[j], *best_vals)
       
        all_fits.append(fit)
        params.append(best_vals)
    

    return all_fits, params
    
def gaussian(x, amp, cen, wid):
    return (amp/(np.sqrt(2*np.pi)*wid)) * np.exp(-(x-cen)**2 /(2*wid**2))

def double_gaussian(x, amp, wid, cen, amp2, wid2, cen2):
    double_gaussian = (amp/(np.sqrt(2*np.pi)*wid)) * np.exp(-(x-cen)**2 /(2*wid**2)) \
                    + (amp2/(np.sqrt(2*np.pi)*wid2)) * np.exp(-(x-cen2)**2 /(2*wid2**2))
    return double_gaussian

def Lorentzian(x,amp,cen,wid):
    #Lorentzian = amp*(wid)**2/(wid**2+(2*x-2*cen)**2)
    Lorentzian = amp*(wid)/(wid**2+(x-cen)**2)/3.14159
    return Lorentzian
#Function that fits the Bremsstrahlung 
def function(x, A, B):
    return B/x - A 

#Function that fits the Bremsstrahlung (this is the code that is used)
def fcn2min(params,x,data):
    if PC == False:
        A = params['A']
        B = params['B']
        """
        F = params['F']
        G = params['G']
        """
        model = B/x - A 
    else:
        A = params['A']
        B = params['B']
        F = params['F']
        G = params['G']
        model =  B/x - A + F/x/x + G/x/x/x
    """
    #Add F and G
    F = params['F']
    G = params['G']
    """

    return model - data

#helps devide the file so that the spectrum for each burst is retrieved.
def get_spectrum_per_burst(spectra, ontime, offtime):
    spectrums = []
    tmp_spectrum= np.zeros(len(spectra[0]))
    i = 0
    time_to_next_burst = ontime
    while i < len(spectra):
        if i < time_to_next_burst:
            tmp_spectrum += spectra[i]
            i += 1
        else:
            spectrums.append(tmp_spectrum)
            tmp_spectrum = np.zeros(2048)
            time_to_next_burst += ontime + offtime
            i += offtime
    return spectrums
    
def get_experiment_name(path):
    exp_name = path.split('.zip')[0]
    return exp_name[-5:-1] + exp_name[-1]

if Pegasus == True:
    my_exp = Experiment(path)
    OffSet = 10
    expname = get_experiment_name(path)
    ontime = int(my_exp.motor_driver.on_time)
    offtime = int(my_exp.motor_driver.off_time)
    spectrums = my_exp.x_ray_detector.spectrum_over_time
    spectrums_per_burst = get_spectrum_per_burst(spectrums, ontime, offtime)
    axis = my_exp.x_ray_detector.spectrum_axis#*1.14043-0.40005
else:
    

try:
    os.system('mkdir "{}"'.format(str(save_path)).replace('/', '\\'))
except:
    pass

#print 'Function it is fitting according to: B/x + A + F*exp(-Gx)'
#print 'Order of coefficients: A, B, F, G'
i = 0
all_parameters = []
dist = len(spectrums_per_burst)
#dist = 31
#dist = 25

#runs through the bursts
for i in range(0,dist):
    spectrum = spectrums_per_burst[i]
    spectrum_without_peaks_old, axis_without_peaks_old, peaks, peaks_axis= get_spectrum_without_peaks(spectrum, axis)
    #spectrum_without_peaks_old = spectrum
    #axis_without_peaks_old = axis

    axis_without_peaks=axis_without_peaks_old[OffSet:len(axis_without_peaks_old)]
    spectrum_without_peaks=spectrum_without_peaks_old[OffSet:len(spectrum_without_peaks_old)]


    """
    # Measure the total amount of counts in the spectrum
    FluxTotal = 0
    spectrum_flux = spectrum_without_peaks
    for countsTot in spectrum_flux:
        FluxTotal += countsTot
    spectrum_without_peaks = spectrum_flux/FluxTotal
  
    for fluxs in spectrum_flux:
        spectrum_without_peaks.append(fluxs/FluxTotal)
    
    """
    
    air_be_kapton_coefficients=sct.total_transmission_factor(axis_without_peaks) #The attenuation Coefficient (uses custom code for this, data retrieved from NIST)
    
    spectrum_without_filter = spectrum_without_peaks / air_be_kapton_coefficients # reverse filter the spectrum


    energies_you_want_to_fit_to=axis_without_peaks
    spectrum_element = []
    axis_element = []

   # This for loop make the matrix for the fitting
    for energy in energies_you_want_to_fit_to:    
        index = get_index_of_energy(axis_without_peaks, energy) 
        spectrum_element.append(spectrum_without_filter[index])
        axis_element.append(axis_without_peaks[index])
    print 'The fitting parameters for burst %i'%i
    try:
    
        if PC == False:
            params = Parameters()
            params.add('A', value = 10, min =1, max=1000)
            params.add('B', value = 100, min =1, max=50000)
            """
            params.add('F',value = 10, min = 0, max = 10000)
            params.add('G',value = 10, min = 0, max = 10000)
            """
        else:
            params = Parameters()
            params.add('A', value = 10, min =0, max=1000)
            params.add('B', value = 100, min =0, max=5000)
            params.add('F', value =1, min =-100000, max=10000000000)
            params.add('G', value = 1, min =0, max=10)

        # fitting the data
        data = spectrum_element
        
        minner = Minimizer(fcn2min, params, fcn_args=(axis_without_peaks,data), fcn_kws={})
        
       
        #result = minner.minimize(method ='differential_evolution')
        #result = minner.minimize(method ='Nelder-Mead')
        result = minner.minimize(method ='least_squares')
        
      
        final = result.residual  + data
      
        spectrum_fit = final
        spectrum_fit_matrix = final*air_be_kapton_coefficients
        
        # putting the parameters together. PC stands for pressure control on or off. Ignore this, it isn't needed
        if PC == False:
            pars = result.params
            ChiSq = result.redchi
            a_v = pars['A'].value
            b_v = pars['B'].value
            """
            f_v = pars['F'].value
            g_v = pars['G'].value
            """
            v_Energy = b_v/a_v
            #print a_v
            fit_params=[a_v,b_v,v_Energy,ChiSq,i]
        else:
            pars = result.params
            ChiSq = result.redchi
            a_v = pars['A'].value
            b_v = pars['B'].value
            v_Energy = b_v/a_v
            f_v = pars['F'].value
            g_v = pars['G'].value
            fit_params=[a_v,b_v,f_v,g_v,v_Energy,ChiSq,i] 
    
        axis_chopped = []
        for axe in axis:
            if axe > 1:
                if axe < 100:
                    axis_chopped.append(axe)
        
        
        air_be_kapton_coefficients2=sct.total_transmission_factor(axis_chopped)
              
        bigFit = []
        

        for axe in axis_chopped:
            bigFit.append(function(axe,a_v,b_v))
        
        bigFit = bigFit*air_be_kapton_coefficients2

        """
        f_v = pars['F'].value
        g_v = pars['G'].value
        fit_params=[a_v,b_v,f_v,g_v,v_Energy,ChiSq]      
        #fit_params=[a_v,b_v,f_v,v_Energy,ChiSq]      
        """
        #print fit_params
        #all_parameters.append(fit_params)
        
        print v_Energy
        print report_fit(result)

        #print "before all fits"
        ##all_fits, fit_peak_params = get_peak_fits(peaks_axis, peaks)
        #print "after all fits"
        #print len(all_fits)

        
        ##for peaks_fit_2 in fit_peak_params:
        ##    for par2 in peaks_fit_2:
        ##        fit_params.append(par2)
       
        #print len(fit_params)
      ##  all_parameters.append(fit_params)
        #print all_parameters

        #print len(all_parameters[i])
        for k in range(len(peaks_axis)):
            start = get_index_of_energy(axis_chopped, peaks_axis[k][2])
            end = get_index_of_energy(axis_chopped, peaks_axis[k][-3])

            new_bigFit = np.append(bigFit[:start], all_fits[k][2:-3])
            new_bigFit = np.append(new_bigFit, bigFit[end:])
            bigFit = new_bigFit

        plt.plot(axis_chopped, bigFit, 'r', zorder=2 ) 
        
       
    except:
        plt.text(50, 50, 'Could not find parameters to fit spectrum')
        print 'Could not fit the spectra %s'%path
        all_fits = [0,0]
    
    #plt.plot(axis_without_peaks, spectrum_without_filter, 'b', zorder=1)
    plt.plot(axis, spectrum, 'b', zorder=1)
    #axis_chopped=np.array(axis_chopped)
    #bigFit=np.array(bigFit)
     
    #plt.plot(axis_without_peaks,spectrum_fit_matrix, 'r', zorder=2 )  
    #plt.plot(axis_without_peaks,final, 'r', zorder=2 )  

    
    #plt.plot(axis_without_peaks,data, 'b', zorder=1)
    #plt.plot(axis_without_peaks,final, 'r', zorder=2 )    
    

    
    plt.xlabel('Energy [keV]')
    plt.ylabel('Counts')
    plt.axis([0,80,0,600])

    plt.title('Spectra fitting for burst nbr: %i'%i)
    plt.show()
    plt.savefig(save_path + '/'+  expname +'_burst_%i.png'%i)
    plt.close()
    specArray =[]
    specArray2=[]

    for oo in range(0,len(axis)):
        specArray.append([axis[oo],spectrum[oo]])
    for nn in range(0,len(axis_chopped)):
        specArray2.append([axis_chopped[nn], bigFit[nn]])



    specDat = pd.DataFrame(specArray,columns = ['axis','counts'])
    specDat.to_csv(save_path +'/spect %i'%i+ '%s.dat'%expname, sep='\t')

all_parameters = np.array(all_parameters)
#print len(all_parameters[1])
if PC == False:
    df = pd.DataFrame(all_parameters, columns=['A','B','Energy','ChiSq/DOF','i','Amp','cent','W','Amp2','cent2','W2','Amp3','cent3','W3','Amp4','cent4','W4'])
else:
    df = pd.DataFrame(all_parameters, columns=['A','B','Energy','ChiSq/DOF','i','Amp','cent','W','Amp2','cent2','W2','Amp3','cent3','W3','Amp4','cent4','W4'])

#df.to_csv(save_path + '/%s.dat'%expname, sep='\t')
    
    
    




 
    
    
    
    
    
    
    
    