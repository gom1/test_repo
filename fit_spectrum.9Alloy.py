# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 17:06:13 2017

@author: MalinDahl
"""

import os
import csv
import math
import numpy as np
import pandas as pd
#from Tkinter import Tk
from pegasus import parsing
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#from tkFileDialog import askdirectory
from scipy.interpolate import interp1d
from lmfit import Minimizer, Parameters
from tkFileDialog import askopenfilename, askdirectory
import ag_line_energy as ale

def cdte_efficiency(spectrum, offset, gain):
    energies = [1,1.006,1.006001,1.5,2,3,3.5375,3.727,4,4.018,4.3414,4.612,4.612001,4.9392,5,6,8,10,12,15,20,25,26.711,26.71101,30,31.814,31.81401,35,40,45,50,55,60,70,80,90,100,125]
    eff = [0,0,0,0.0002,.0299,.368,.5276,.5839,.6649,.6676,.716,.7564,.7564,.8053,.8144,.888,.9486,.9701,.9779,.9857,.9895,.9895,.9909,.9909,.9916,.9916,.9916,.9916,.9923,.9923,.9927,.9927,.9927,.9879,.9685,.9265,.8646,.6854]
    energy_log = [np.log10(x) for x in energies]
    attenuation_log = [np.log10(x) for x in eff]   
    interpolation_function_log = interp1d(energy_log, attenuation_log, kind='linear', assume_sorted = True)     
    Mat = []
    nChannels = len(spectrum)
    start = 0
    stop = nChannels 
    for i in range(start, stop):
        
        if (i+1)*gain+offset >= 1:
            energy_param = float((i+1)*gain+offset)
            spectrumData = float(spectrum[i])
            
            Mat.append(spectrumData/(10**interpolation_function_log(np.log10(energy_param))))
           
        else:
            Mat.append(0)
    try:
        return Mat
    except:
        raise


def sipin_efficiency(spectrum, offset, gain, fDetThick=0.5):
    """    
    nChannels = len(spectrum)
    start = 0
    stop = nChannels - 1
    ltk = 1.8389
    ltc = 10.6422
    ltn = 2.7345
    fEfficiency = 1.0 
    fAlpha = 1.0
    fDensity = 2.33
    spectrumCorrected = [0] * nChannels
    axisCorrected = [0] * nChannels
    fEfficiencies = [0] * nChannels
    for i in range(start, stop):
        fKev = offset + gain * i
        axisCorrected[i] = fKev
        if fKev > 2.0:
            temp1 = 12.3981 / float(fKev)
            temp2 = float(ltn)
            temp3 = math.pow(temp1, temp2)
            mac = ltk * ltc * temp3
            if(mac > 1e-3):
                temp = mac * fDensity * fDetThick * fAlpha / 10.0
                if temp < 100.0:
                    temp = math.exp(-temp)
                else:
                    temp = 0.0
                fEfficiency = float(1.0 - temp)
            else:
                fEfficiency = 1.0
            fEfficiencies[i] = fEfficiency
            
            spectrumCorrected[i] = float(spectrum[i])
            if fEfficiency > 1e-4:
                counts = float(spectrum[i])/fEfficiency
                if counts < 1e9:
                    spectrumCorrected[i] = counts + 0.5
        """
    #energies = [0.586,0.626,0.669,0.715,0.765,0.818,0.874,0.934,0.999,1.07,1.14,1.22,1.3,1.39,1.49,1.59,1.7,1.8,1.82,1.83,1.84,1.85,1.88,1.95,2.08,2.22,2.38,2.54,2.72,2.9,3.11,3.32,3.55,3.79,4.06,4.33,4.63,4.95,5.3,5.66,6.05,6.47,6.92,7.39,7.9,8.45,9.03,9.65,10.3,11,11.8,12.6,13.5,14.4,15.4,16.5,17.6,18.8,20.1,21.5,23,24.6,26.3,28.1,30,32.1,34.3,36.7,39.2,41.9,44.8,47.9,51.2,54.7,58.5,62.5,66.8,71.5,76.4,81.7,87.3,93.3,99.8,107]
    #eff = [0.00000238,0.000022,0.00014,0.000649,0.00232,0.00668,0.016,.0331,.0602,.1,.152,.214,.283,.356,.429,.5,.567,.619,.629,.632,.636,.576,.588,.621,.675,.725,.768,.805,.837,.864,.887,.907,.923,.937,.948,.957,.965,.971,.976,.98,.984,.987,.989,.99,.991,.99,.985,.973,.952,.92,.875,.92,.756,.686,.615,.545,.477,.415,.358,.308,.263,.224,.191,.162,.138,.117,.0997,.0854,.0735,.0634,.0551,.0483,.0426,.038,.0341,.0309,.0283,.026,.0242,.0226,.0213,.0201,.0191,.0183]                    
    
    energies = [0.512, 0.548, 0.586, 0.626, 0.669, 0.715, 0.765, 0.818, 0.874, 0.934, 0.999, 1.07, 1.14, 1.22, 1.3, 1.39, 1.49, 1.59, 1.7, 1.8, 1.82, 1.83, 1.84, 1.85, 1.88, 1.95, 2.08, 2.22, 2.38, 2.54, 2.72, 2.9, 3.11, 3.32, 3.55, 3.79, 4.06, 4.33, 4.63, 4.95, 5.3, 5.66, 6.05, 6.47, 6.92, 7.39, 7.9, 8.45, 9.03, 9.65, 10.3, 11, 11.8, 12.6, 13.5, 14.4, 15.4, 16.5, 17.6, 18.8, 20.1, 21.5, 23, 24.6, 26.3, 28.1, 30, 32.1, 34.3, 36.7, 39.2, 41.9, 44.8, 47.9, 51.2, 54.7, 58.5, 62.5, 66.8, 71.5, 76.4, 81.7, 87.3, 93.3, 99.8, 107, 114, 122, 130, 139, 149, 159, 170, 182, 194, 208, 222, 237, 254,271, 290, 310, 332, 354, 379, 405]
    #eff = [0.00000000677, 1.65*10**(-7), 2.38*10**(-6), 0.000022, 0.00014, 0.000649, 0.00232, 0.00668, 0.016, 0.0331, 0.0602, 0.1, 0.152, 0.214, 0.283, 0.356, 0.429, 0.5, 0.567, 0.619, 0.629, 0.632, 0.636, 0.576, 0.588, 0.621, 0.675, 0.725, 0.768, 0.805, 0.837, 0.864, 0.887, 0.907, 0.923, 0.937, 0.948, 0.957, 0.965, 0.971, 0.976, 0.98, 0.984, 0.987, 0.989, 0.991, 0.992, 0.992, 0.989, 0.981, 0.965, 0.939, 0.901, 0.851, 0.791, 0.724, 0.654, 0.583, 0.514, 0.449, 0.389, 0.335, 0.288, 0.246, 0.21, 0.179, 0.152, 0.129, 0.11, 0.0944, 0.0813, 0.0702, 0.0611, 0.0535, 0.0473, 0.0421, 0.0378, 0.0343, 0.0313, 0.0289, 0.0268, 0.0251, 0.0236, 0.0223, 0.0212, 0.0203, 0.0195, 0.0187, 0.0181, 0.0175, 0.0169, 0.0164, 0.0159, 0.0155, 0.0151, 0.0147, 0.0143, 0.0139, 0.0135, 0.0132, 0.0128, 0.0125, 0.0121, 0.0118, 0.0115, 0.0111]
    eff = [6.77*10**(-9), 1.65*10**(-7),2.38*10**(-6), 0.000022, 0.00014, 0.000649, 0.00232, 0.00668, 0.016, 0.0331, 0.0602, 0.1, 0.152, 0.214, 0.283, 0.356, 0.429, 0.5, 0.567, 0.619, 0.629, 0.632, 0.636, 0.576, 0.588, 0.621, 0.675, 0.725, 0.768, 0.805, 0.837, 0.864, 0.887, 0.907, 0.923, 0.937, 0.948, 0.957, 0.965,0.971, 0.976, 0.98, 0.984, 0.987, 0.989, 0.991, 0.992, 0.991, 0.988, 0.98, 0.963, 0.934, 0.893, 0.84, 0.777, 0.707, 0.634, 0.56, 0.489, 0.422, 0.361, 0.307, 0.259, 0.217, 0.181, 0.15, 0.124, 0.102, 0.0832, 0.0681, 0.0556, 0.0452, 0.0367, 0.0298, 0.0241, 0.0195, 0.0158, 0.0128, 0.0104, 0.00838, 0.00678, 0.00548, 0.00443, 0.00358, 0.0029, 0.00234, 0.00189, 0.00153, 0.00124, 0.001, 0.000809, 0.000654, 0.000528, 0.000427, 0.000345, 0.000279, 0.000226, 0.000182, 0.000147, 0.000119, 0.0000963, 0.0000779, 0.0000629, 0.0000509, 0.0000411, 0.0000332]
    #eff = [3.83*10**(-5),2.13*10**(-4),0.000896,0.00297,0.00806,0.0185,0.0369,0.0653,0.105,0.156,0.215,0.284,0.356,0.43,0.501,0.567,0.629,0.684,0.732,0.768,0.775,0.777,0.779,0.704,0.713,0.737,0.776,0.812,0.843,0.868,0.89,0.908,0.924,0.937,0.948,0.957,0.965,0.971,0.976,0.98,0.984,0.987,0.989,0.991,0.992,0.994,0.994,0.994,0.991,0.983,0.966,0.94,0.902,0.851,0.792,0.725,0.654,0.583,0.514,0.449,0.389,0.335,0.288,0.246,0.21,0.179,0.152,0.129,0.110,0.0945,0.0813,0.0703,0.0611,0.0535,0.0473,0.0421,0.0379,0.0343,0.0314,0.0289,0.0268,0.0251,0.0236,0.023,0.0212,0.0203,0.0195,0.0187,0.0184,0.0175,0.0169,0.0164,0.0159,0.0155,0.0151,0.0147,0.0143,0.0139,0.0135,0.0132,0.0128,0.0125,0.0121,0.0118,0.0115,0.0111]
    energy_log = [np.log10(x) for x in energies]
    attenuation_log = [np.log10(x) for x in eff]   
    interpolation_function_log = interp1d(energy_log, attenuation_log, kind='linear', assume_sorted = True)     
    Mat = []
    nChannels = len(spectrum)
    start = 0
    stop = nChannels 
    for i in range(start, stop):
        #print i*gain+offset
        if (i+1)*gain+offset >= 0.512:
            energy_param = float((i+1)*gain+offset)
            spectrumData = float(spectrum[i])
            #print spectrumData*10**interpolation_function_log(energy_param)
            Mat.append(spectrumData/(10**interpolation_function_log(np.log10(energy_param))))
            #print 10**interpolation_function_log(energy_param)
            #print energy_param
            #Mat.append(1)
        else:
            Mat.append(0)
    try:
        return Mat
    except:
        raise


#Copied function
def mu_Al(photon_energy):
    energies = [1,1.5,1.5596,2,3,4,5,6,8,10,15,20,30,40,50,60,80,100,150,200]
    mu_m = [1185,402.2,362.1,2263,788,360.5,193.4,115.3,50.33,26.23,7.955,3.441,1.128,0.5685,0.3681,0.2778,0.2018,0.1704,0.1378,0.1223]
    energy_log = [np.log10(x) for x in energies]
    attenuation_log = [np.log10(x) for x in mu_m]   
    interpolation_function_log = interp1d(energy_log, attenuation_log, kind='linear', assume_sorted = True)     
    Mat = []
    for energy in photon_energy:
        Mat.append(10**interpolation_function_log(np.log10(energy))) 
    try:
        return Mat
    except:
        raise

#Copied function
def mu_Kapton(photon_energy): 
    with open('C:\Users\Tribogenics\Desktop\KapArray2.txt','rb') as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        energies = []
        mu_m = []    
        for row in tsvin:
            energies.append(float(row[0])/1000) 
            mu_m.append(10000/float(row[1]))
    energy_log = [np.log10(x) for x in energies]
    attenuation_log = [np.log10(x) for x in mu_m]
    interpolation_function_log = interp1d(energy_log, attenuation_log, kind='linear', assume_sorted = True)
    Mat = []
    for energy in photon_energy:
        if energy > 30.00:
            Mat.append(1)
        else:
            Mat.append(10**interpolation_function_log(np.log10(energy)))
    try:
        return Mat
    except:
        raise

#Copied function 
def mu_air(photon_energy): #valid for 1-200keV
    energies = [1,1.5,2,3,3.2029,4,5,6,8,10,15,20,30,40,50,60,80,100,150,200]
    mu_m = [3606,1191,527.9,162.5,148.5,77.88,40.27,23.41,9.921,5.12,1.614,0.7779,0.3538,0.2485,0.208,0.1875,0.1662,0.1541,0.1356,0.1233]
    energy_log = [np.log10(x) for x in energies]
    attenuation_log = [np.log10(x) for x in mu_m]   
    interpolation_function_log = interp1d(energy_log, attenuation_log, kind='linear', assume_sorted = True)     
    Mat = []
    for energy in photon_energy:
        Mat.append(10**interpolation_function_log(np.log10(energy))) 
    try:
        return Mat
    except:
        raise

def mu_Ag(photon_energy): #valid for 1-200keV
    energies = [1,1.5,2,3,3.3511,3.4363,3.5237,3.66203,3.8058,4,5,6,8,10,15,20,25.514,30,40,50,60,80,100,150,200]
    mu_m = [7039,2790,1401,513.6,388.7,1200,1126,1409,1282,1305,738.5,461,216.4,119.3,39.98,18.36,9.527,36.68,17.2,9.444,5.766,2.651,1.47,0.5426,0.2972]
    energy_log = [np.log10(x) for x in energies]
    attenuation_log = [np.log10(x) for x in mu_m]   
    interpolation_function_log = interp1d(energy_log, attenuation_log, kind='linear', assume_sorted = True)     
    Mat = []
    for energy in photon_energy:
        if energy > 100.00:
            Mat.append(1)
        else:
            Mat.append(10**interpolation_function_log(np.log10(energy))) 
    try:
        return Mat
    except:
        raise


#Copied function 
def mu_Be(photon_energy): #valid for 1-100keV
    energies = [1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0]
    #mu_m = [604.1, 179.7, 74.69, 21.27, 8.685, 4.369, 2.527, 1.124, 0.6466, 0.307, 0.2251, 0.1792, 0.164, 0.1554, 0.1493, 0.1401, 0.1328]
    mu_m = [603.5,179.1,74.22,20.9,8.367,4.081,2.26,.8839,.4255,.1143,0.0478,0.01898,0.01438,0.01401,0.01468,0.01658,0.01836]
    energy_log = [np.log10(x) for x in energies]
    attenuation_log = [np.log10(x) for x in mu_m]
    interpolation_function_log = interp1d(energy_log, attenuation_log, kind='linear', assume_sorted = True)
    Mat = []
    for energy in photon_energy:
        if energy > 100.00:
            Mat.append(1)
        else:
            Mat.append(10**interpolation_function_log(np.log10(energy))) 
    try:
        return Mat
    except:
        print "Error occured for photon energy: " + str(photon_energy)
        raise

#Copied function 
def Al_window_transmission_factor(photon_energy, t=0.0005*2.54, rho=2.699):
    x = rho * t
    mu = mu_Al(photon_energy)
    AlMat = []
    for i in range(0,len(photon_energy)):
        AlMat.append(np.exp(-mu[i]*x))
    return AlMat

#Copied function 
def CdTe_detector(photon_energy):
    with open(parsing.detector_data_path("CdTe-Parameter.dat")) as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        eff_data = [(float(x), float(y)) for (x, y) in reader]

    efficiency_correction = interp1d([x for (x, y) in eff_data], [y for (x, y) in eff_data], kind='linear', assume_sorted = True)
    Mat = []
    for energy in photon_energy:
        Mat.append(efficiency_correction(energy))
    return Mat

#Copied function 
def target_Kapton_transmission_factor(photon_energy):
    t = .001*2.54
    x = t
    mu = mu_Kapton(photon_energy)
    KapMat = []
    for i in range(0,len(photon_energy)):
        KapMat.append(np.exp(-mu[i]*x))
    return KapMat

#Copied function 2.225
#Original air = 2.225
def air_column_transmission_factor(photon_energy, t=2.225, rho=0.0012047):
    #t = 3.7*2.54
    #t = 1.75*2.54
    x = rho * t
    mu = mu_air(photon_energy)
    AirMat = []
    for i in range(0,len(photon_energy)):
        AirMat.append(np.exp(-mu[i]*x))
    return AirMat

def air_column_transmission_factorCdTe(photon_energy, t=23.5, rho=0.0012047):
    x = rho * t
    mu = mu_air(photon_energy)
    AirMat = []
    for i in range(0,len(photon_energy)):
        AirMat.append(np.exp(-mu[i]*x))
    return AirMat

#Copied function 
def Be_window_transmission_factor(photon_energy, t=0.01*2.54, rho=1.848):
    x = rho * t
    mu = mu_Be(photon_energy)
    BeMat = []
    for i in range(0,len(photon_energy)):
        BeMat.append(np.exp(-mu[i]*x))
    return BeMat

def Ag_transmission_factor(photon_energy, t=0.000075, rho=10.5):
    x = rho * t
    mu = mu_Ag(photon_energy)
    AgMat = []
    for i in range(0,len(photon_energy)):
        AgMat.append(np.exp(-mu[i]*x))
    return AgMat

#Copied function 
def total_transmission_factor(photon_energy):
    Be = np.array(Be_window_transmission_factor(photon_energy))
    Air = np.array(air_column_transmission_factor(photon_energy))
    Kapton = np.array(target_Kapton_transmission_factor(photon_energy))
    Agf = np.array(Ag_transmission_factor(photon_energy))
    AlF = np.array(Al_window_transmission_factor(photon_energy))
    #CdTe = np.array(CdTe_detector(photon_energy))

    return Be*Air*Kapton#*AlF
    #return Agf

def get_index_of_energy(spectrum_axis, energy):
    index = (np.abs(np.array(spectrum_axis) - energy)).argmin()
    return index

def get_spectrum_without_peaks(spectrum, axis, peaks_to_remove):
    peaks = []
    peaks_axis = []
    for peak in peaks_to_remove:
        # Indices of where peak starts and end
        start_index = get_index_of_energy(axis, peak[0])
        end_index = get_index_of_energy(axis, peak[1])

        #Spectrum up until peak start
        beginning_spectrum = np.array(spectrum[:start_index])
        #Spectrum from when peak ends
        end_spectrum = np.array(spectrum[end_index:])
        #Spectrum of the peak
        peak = np.array(spectrum[start_index:end_index])
        
        #Axis up until peak starts
        beginning_axis = np.array(axis[:start_index])
        #Axis from when peak ends
        end_axis = np.array(axis[end_index:])
        #Axis of the peak
        peak_axis = np.array(axis[start_index:end_index])

        #Spectrum and axis of all but the peaks
        spectrum = np.append(beginning_spectrum, end_spectrum)
        axis = np.append(beginning_axis, end_axis)

        #Spectrum and axis of peak
        peak_spectrum = np.array(peak)
        peak_axis_spectrum=np.array(peak_axis)
        
        #Add all the peaks togehter
        peaks.append(peak_spectrum)
        peaks_axis.append(peak_axis_spectrum)
    return spectrum, axis, peaks, peaks_axis

def get_peak_fits(peaks_axis, peaks):
    all_fits = []
    params = []
    for j in range(len(peaks)):
        if j == 0:
            best_vals, covar = curve_fit(Lorentzian, peaks_axis[j], peaks[j], [10,20,10])
            fit = Lorentzian(peaks_axis[j], *best_vals)
        else:
            best_vals, covar = curve_fit(Lorentzian, peaks_axis[j], peaks[j], [10, 20, 10])
            fit = Lorentzian(peaks_axis[j], *best_vals)       
        all_fits.append(fit)
        params.append(best_vals)
    return all_fits, params
    
#Not used function
def gaussian(x, amp, cen, wid):
    return (amp/(np.sqrt(2*np.pi)*wid)) * np.exp(-(x-cen)**2 /(2*wid**2))

#Not used function
def double_gaussian(x, amp, wid, cen, amp2, wid2, cen2):
    double_gaussian = (amp/(np.sqrt(2*np.pi)*wid)) * np.exp(-(x-cen)**2 /(2*wid**2)) \
                    + (amp2/(np.sqrt(2*np.pi)*wid2)) * np.exp(-(x-cen2)**2 /(2*wid2**2))
    return double_gaussian

#Not used function
def Lorentzian(x,amp,cen,wid):
    Lorentzian = amp*(wid)/(wid**2+(x-cen)**2)/3.14159
    return Lorentzian

#Function that fits the Bremsstrahlung 
def Kramers_law(x, A, B):
    return B/x - A 

#Function that fits the Bremsstrahlung (this is the code that is used)
def fcn2min(params,x,data):
    A = params['A']
    B = params['B']
    if PC == False:
        model = B/x - A 
    else:
        F = params['F']
        G = params['G']
        model =  B/x - A + F/x/x + G/x/x/x
    return model - data

#Not used
def get_spectrum_per_burst(spectra, ontime, offtime):
    '''
    If you have spectrums from an analysis file, this funtion helps to
    divide all the spectrums per burst into an array.
    '''
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

def binomialSmooth(data, degree):
    for _ in range(degree):
        smoothed = np.ones(len(data))
        smoothed[0] = data[0]
        smoothed[len(data) - 1] = data[-1]
        for k in range(1, len(data) - 2):
            smoothed[k] = (data[k-1] + 2 * data[k] + data[k+1]) / 4
        data = smoothed
    return smoothed

"""
Main code starts here
"""
# Reading the spectrum file
path_to_spectrum_file = "C:\Users\MalinDahl\Desktop\Spectrum4.dat"
#path_to_spectrums = "C:\Users\Tribogenics\Desktop\processed_steel_alloys\Results\\"
#path_to_spectrums = "C:\Users\Tribogenics\Desktop\processed_steel_alloys\\"
path_to_spectrums = "C:\Users\Tribogenics\Desktop\Data Set Test - to delete"
#path_to_spectrums = askdirectory(title='Pick an analysis file created by Mustang Aquire')
path_to_spectrums=path_to_spectrums + "\\"
print path_to_spectrums
#path_to_spectrums = "C:\Users\Tribogenics\Desktop\TestFiles\\"
PC = False #Pressure control
energies = []
AgEnergies=[]
a = []
b = []
alloys = []
for new_path in os.listdir(path_to_spectrums):
    alloy = new_path
    new_path += "\\"
    i = 0
    FileType = ""
    for f in os.listdir(path_to_spectrums + new_path):
        alloys.append(new_path + str(i))
        print f
        spectrum = np.array(pd.read_csv(path_to_spectrums + new_path + f, header=None))
        start_index = np.where(spectrum == "<<DATA>>")[0][0] + 1
        end_index = np.where(spectrum == "<<END>>")[0][0] 
        spectrum = spectrum[start_index: end_index].flatten()
        print f[0]
        if f[0] == "Y":
            axis = np.linspace(0.159497, 90.4052, num=2048)
            offset = 0.159497
            gain = 0.0440653
            #print axis
        elif f[0] == "T":
            #axis = np.linspace(0.0969676, 83.7156, num=2048)
            #offset = 0.0969676
            #gain = 0.0408294
            #offset = -0.0359929
            #gain = 0.0542093 
            offset = -0.122238
            gain =0.0467921
            FileType = "T"
            axis = np.linspace(offset, gain*2048+offset, num=2048)
        elif f[0] == "N":
            axis = np.linspace(0.01813, 75.2806964, num=2048)
            offset = 0.01813
            gain = 0.0367493
        elif f[0] == "F":
            axis = np.linspace(-0.006637, 76.150291, num=2048)
            offset = -0.006637
            gain = 0.037186
        elif f[0] == "B":
            axis = np.linspace(-0.00356152, 87.0662, num=2048)
            offset = -0.00356152
            gain = 0.0425145 
        elif f[0] == "A":
            #offset = 0.21817 
            #gain = 0.0429977
            #offset = -0.96
            #gain = 0.0436818
            offset = 0.0185049  
            gain = 0.0369941
            FileType = "A"
            #offset =-0.0673333
            #gain = 0.0411342
            #offset =0.282667
            #gain = 0.0404848
            axis = np.linspace(offset, gain*2048+offset, num=2048)
        else:
            print f
        sipin_axis = axis
        #Only look at values between 1-100 keV due to attenuation coefficients
        start_energy = 15
        end_energy = 35
        start_energy =4.1
        end_energy = 40
        start_index = np.argmax(axis > start_energy)
        end_index = np.argmax(axis > end_energy) - 1
        start_index_old = np.argmax(axis > 6)
        end_index_old = np.argmax(axis > 80) - 1
        axis_orig=axis
        axis_old = axis[start_index_old: end_index_old]
        axis = axis[start_index: end_index]
        old_spectrum = spectrum
        if FileType == "A":
            sipin_spectrum = sipin_efficiency(spectrum, offset, gain)
        elif FileType == "T":
            sipin_spectrum = cdte_efficiency(spectrum, offset, gain)
        """
        OnesArray = []
        for i in range(0,len(sipin_spectrum)):
            OnesArray.append(100)
        Ones_SiPin = sipin_efficiency(OnesArray, offset, gain)
       
        plt.plot(axis_orig,Ones_SiPin)
        plt.xlim([0, 40])
        plt.ylim([0,20])
        plt.show()
        """
        #print len(axis_orig)
        #print len(spectrum)
        
        #smoothing
        sipin_spectrum=binomialSmooth(sipin_spectrum,4)
        
        spectrum = sipin_spectrum[start_index:end_index]
        spectrum_old = sipin_spectrum[start_index_old:end_index_old]
        #spectrum = spectrum[0::10]
        #axis = axis[0::10]
        # Remove both Ag peaks
        peaks_to_remove = [(20, 26)]
        spectrum_without_peaks, axis_without_peaks, peaks, peaks_axis = get_spectrum_without_peaks(spectrum, axis, peaks_to_remove)
        
        # Reverse the filters
        air_be_kapton_coefficients = total_transmission_factor(axis_without_peaks)
        spectrum_without_filter = spectrum_without_peaks / air_be_kapton_coefficients 

        if PC == False:
            params = Parameters()
            params.add('A', value = 65, min=0, max=10000)
            params.add('B', value = 100, min=0, max=500000)
        else:
            params = Parameters()
            params.add('A', value = 10, min =0, max=1000)
            params.add('B', value = 100, min =0, max=5000)
            params.add('F', value = 1, min =-100000, max=10000000000)
            params.add('G', value = 1, min =0, max=10)
        
        # fitting the data
        data = spectrum_without_filter
        #print axis_old
        AgLineParam = ale.getEnergy(spectrum_old, axis_old)
        #print "AgLineInfo"
        #print AgLineParam
        try:
            #data=binomialSmooth(data,3)
            #Set the initial conditions for the fits
            minner = Minimizer(fcn2min, params, fcn_args=(axis_without_peaks,data), fcn_kws={})
            result = minner.minimize(method ='least_squares')
            #result = minner.minimize(method ='Nelder-Mead',options={'maxiter':10000000000000000000000})
            #The final result after least square minimization.
            final = result.residual + data
            #Adding the filter back to spectrum
            spectrum_fit_matrix = final#*air_be_kapton_coefficients
            # putting the parameters together. PC stands for pressure control on or off. Ignore this, it isn't needed
            if PC == False:
                pars = result.params
                ChiSq = result.redchi
                a_v = pars['A'].value
                b_v = pars['B'].value
                v_Energy = b_v/a_v
                fit_params=[a_v,b_v,v_Energy,ChiSq]
            else:
                pars = result.params
                ChiSq = result.redchi
                a_v = pars['A'].value
                b_v = pars['B'].value
                v_Energy = b_v/a_v
                f_v = pars['F'].value
                g_v = pars['G'].value
                fit_params=[a_v,b_v,f_v,g_v,v_Energy,ChiSq] 
            
            # Kramer's law fit to energies with the calculated parameters
            bigFit = []
            for energy in axis_without_peaks:
                bigFit.append(Kramers_law(energy,a_v,b_v))
            bigFit = bigFit*air_be_kapton_coefficients
            energies.append(v_Energy)
            AgEnergies.append(AgLineParam[0])
            #Print the results to the terminal
            print v_Energy
            print AgLineParam[0]
            print AgLineParam[2]
            print 'params A'
            print a_v
            print 'params B'
            print b_v
            a.append(a_v)
            b.append(b_v)
            """
            Final thing to do is append the peaks back to the spectrum
            """    
        except:
            plt.text(50, 50, 'Could not find parameters to fit spectrum')
            print 'Could not fit the spectra'
            all_fits = [0,0]
        i += 1
     

        fig, ax1 = plt.subplots(figsize=[15, 7])
        ax2 = ax1.twinx()
        #ax1.plot(np.linspace(-0.006637, 76.150291, num=2048), sipin_spectrum, 'r', label="Original Spectrum", linewidth=2)
        #ax1.plot(axis_without_peaks, spectrum_without_peaks, 'g', label="Spectrum No Peaks")
        #ax1.plot(axis_orig, old_spectrum, 'g', label="Original Spectrum pre sipin", linewidth=2)
        ax1.plot(axis_orig, sipin_spectrum, 'r', label="Original Spectrum", linewidth=2)
        ax1.plot(axis_without_peaks,data, 'b', label="Fitting Spectrum", linewidth=2)
        #ax1.plot(axis_orig, Ones_SiPin, 'b', label="Si Pin Eff", linewidth=2)
        #ax1.plot(axis_without_peaks, spectrum_without_filter, 'k', label="Spectrum After 'Removing' Filter")
        #ax1.plot(axis_without_peaks, final, 'y', label="Spectrum Before Adding the Filter Back on")
        ax1.plot(axis_without_peaks, spectrum_fit_matrix, 'm', label="Result of Spectrum Fitting")
        #ax1.plot([v_Energy] * 1000, range(1000), 'k', label="Peak Energy")
        #ax1.set_ylim([0, 2000])
        ax1.set_ylim([0, 2000])
        ax1.set_xlim([0, 80])
        ax1.set_ylabel("Counts")
        ax1.set_xlabel("Energy [keV]")
        plt.title("%s Removed All Ag Peaks. Range [%i:%i] keV. Energy = %.2f. AgEnergy = %.2f"%(alloy, start_energy, end_energy, v_Energy, AgLineParam[0]), fontsize=24)
        ax1.legend()
        plt.show()    


plt.figure(figsize=[12, 7])
plt.plot(energies, 'o')
plt.title("Range [%i, %i]keV analysis"%(start_energy, end_energy), fontsize=24)
plt.xlabel("Different experiments", fontsize=20)
plt.ylabel("Energy [keV]", fontsize=20)
plt.show()

plt.figure(figsize=[12, 7])
plt.plot(AgEnergies, 'o')
plt.title("Range [%i, %i]keV analysis"%(start_energy, end_energy), fontsize=24)
plt.xlabel("Different experiments", fontsize=20)
plt.ylabel("Ag Energy [keV]", fontsize=20)
plt.show()

"""
df = pd.DataFrame([energies]).T
df2 = pd.DataFrame([AgEnergies]).T
df3 = pd.DataFrame([alloys]).T
df.to_csv("C:\Users\Tribogenics\Desktop\New_files\PeakEnergies.csv", index=False, sep="\t")
df2.to_csv("C:\Users\Tribogenics\Desktop\New_files\AgEnergies.csv", index=False, sep="\t")
df3.to_csv("C:\Users\Tribogenics\Desktop\New_files\Alloys.csv", index=False, sep="\t")
"""
    
    
    
    
    
    
    
    