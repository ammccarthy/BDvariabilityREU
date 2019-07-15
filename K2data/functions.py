# Allison McCarthy
# Summer 2019
# This creates functions to use to test the building of the lightcurve and the sin curve fit 

import os
import csv
import lightkurve as lk
from lightkurve import search_targetpixelfile
import numpy as np
from scipy.optimize import leastsq
from tqdm import tqdm
import astropy.units as u
import pylab as plt


### RUN LIGHTCURVE FUNCTION ###

def correct_lightcurve(EPIC=NONE, campaign_num=NONE):
    if campaign_num is not NONE:
        assert EPIC is not NONE, ("You must input an EPIC if you input a campaign number")
        try:
            EPIC = int(EPIC)
            campaign_num = int(campaign_num)
            root_path=s.getcwd()
            dir_path=root_path+str(EPIC)+"_"+str(campaign_num)
            os.mkdir(dir_path)
            os.chdir(dir_path)
        #Actually Running Lightkurve
            tpf = search_targetpixelfile(EPIC, mission='K2', campaign=campaign_num).download()
        #Target Pixel File and Aperture Mask
            user_lc = tpf.to_lightcurve(aperture_mask=tpf.pipeline_mask.astype(bool))
            user_lc = user_lc.remove_nans().remove_outliers()
        #Uncorrected Periodogram
            p = user_lc.to_periodogram(freq_unit=u.microHertz, maximum_frequency=400, minimum_frequency=10)
        #Lomb Scargle Periodogram
            raw_periodogram=lk.periodogram.LombScarglePeriodogram.from_lightcurve(user_lc, minimum_period=0.05, maximum_period =10)
            period_at_max_power_user_lc= raw_periodogram.period_at_max_power
        #Normalize and Correct Lightcurve
            lc = tpf.to_lightcurve().normalize().remove_nans().remove_outliers()
            clc = lc.correct(windows=10).remove_outliers().fill_gaps()
        #Corrected Periodogram
            corrected_periodogram=lk.periodogram.LombScarglePeriodogram.from_lightcurve(clc, minimum_period=0.05, maximum_period =100)
            corrected_period_at_max_power= corrected_periodogram.period_at_max_power
        #Folded Corrected Lightcurve
            corrected_folded_lightcurve = clc.fold(corrected_period_at_max_power.value)
        #Binned Folded Lightcurve
            corrected_bin_folded_lc = corrected_folded_lightcurve.bin(10,method='median')
        #SFF Correcting
            corrector=lk.SFFCorrector(lc)
            new_lc = corrector.correct(lc.centroid_col,lc.centroid_row)
        #Lomb Scargle Periodogram for SFF corrected lightcurve
            SFF_corrected_periodogram=lk.periodogram.LombScarglePeriodogram.from_lightcurve(new_lc, minimum_period=0.05, maximum_period =100)
            SFF_corrected_period=SFF_corrected_periodogram.period_at_max_power
        #SFF Corrected Lightcurve
            SFF_corrected_folded_lightcurve = new_lc.fold(SFF_corrected_period.value)
        #Binned Folded Lightcurve
            SFF_bin_folded_lc = SFF_corrected_folded_lightcurve.bin(10,method='median')
        except:
            print(str(EPIC)+str(campaign_num)+' has a problem')
            os.rename(dir_path,dir_path+'HasAnError')
            pass
    if EPIC is not NONE:
        try:
            EPIC = int(EPIC)
            root_path=s.getcwd()
            dir_path=root_path+str(EPIC)
            os.mkdir(dir_path)
            os.chdir(dir_path)         
        #Actually Running Lightkurve
            tpf = search_targetpixelfile(EPIC, mission='K2').download()
        #Target Pixel File and Aperture Mask
            user_lc = tpf.to_lightcurve(aperture_mask=tpf.pipeline_mask.astype(bool))
            user_lc = user_lc.remove_nans().remove_outliers()
        #Uncorrected Periodogram               
            p = user_lc.to_periodogram(freq_unit=u.microHertz, maximum_frequency=400, minimum_frequency=10)
        #Lomb Scargle Periodogram
            raw_periodogram=lk.periodogram.LombScarglePeriodogram.from_lightcurve(user_lc, minimum_period=0.05, maximum_period =10)
            period_at_max_power_user_lc= raw_periodogram.period_at_max_power
        #Normalize and Correct Lightcurve              
            lc = tpf.to_lightcurve().normalize().remove_nans().remove_outliers()
            clc = lc.correct(windows=10).remove_outliers().fill_gaps()
        #Corrected Periodogram
            corrected_periodogram=lk.periodogram.LombScarglePeriodogram.from_lightcurve(clc, minimum_period=0.05, maximum_period =100)
            corrected_period_at_max_power= corrected_periodogram.period_at_max_power          
        #Folded Corrected Lightcurve             
            corrected_folded_lightcurve = clc.fold(corrected_period_at_max_power.value)
        #Binned Folded Lightcurve
            corrected_bin_folded_lc = corrected_folded_lightcurve.bin(10,method='median')
        #SFF Correcting
            corrector=lk.SFFCorrector(lc)
            new_lc = corrector.correct(lc.centroid_col,lc.centroid_row)
        #Lomb Scargle Periodogram for SFF corrected lightcurve             
            SFF_corrected_periodogram=lk.periodogram.LombScarglePeriodogram.from_lightcurve(new_lc, minimum_period=0.05, maximum_period =100)
            SFF_corrected_period=SFF_corrected_periodogram.period_at_max_power
        #SFF Corrected Lightcurve
            SFF_corrected_folded_lightcurve = new_lc.fold(SFF_corrected_period.value)
        #Binned Folded Lightcurve
            SFF_bin_folded_lc = SFF_corrected_folded_lightcurve.bin(10,method='median')
        except:
            print(str(EPIC)+str(campaign_num)+' has a problem')
            os.rename(dir_path,dir_path+'HasAnError')
            pass


### SIN CURVE FUNCTION ###

def fit_sin_curve(flux):
        
     N=len(bin_folded_lc.flux) # number of data points
     t = np.linspace(-0.5, 0.5, N)

     guess_mean = np.mean(bin_folded_lc.flux)
     guess_std = 3*np.std(bin_folded_lc.flux)/(2**0.5)/(2**0.5)
     guess_phase = 0
     guess_freq = 1
     guess_amp = 1

# we'll use this to plot our first estimate. This might already be good enough for you
     data_first_guess = guess_std*np.sin(t+guess_phase) + guess_mean

# Define the function to optimize, in this case, we want to minimize the difference
# between the actual data and our "guessed" parameters
     optimize_func = lambda x: x[0]*np.sin(x[1]*t+x[2]) + x[3] - bin_folded_lc.flux
     est_amp, est_freq, est_phase, est_mean = leastsq(optimize_func, [guess_amp, guess_freq, guess_phase, guess_mean])[0]

# recreate the fitted curve using the optimized parameters
     data_fit = est_amp*np.sin(est_freq*t+est_phase) + est_mean

# recreate the fitted curve using the optimized parameters
     fine_t = np.arange(-0.5,0.5,0.001)
     data_fit=est_amp*np.sin(est_freq*fine_t+est_phase)+est_mean

### PLOTING FUNCTION ###

def plot_all(EPIC, campaign_num, tpf, user_lc, p, raw_periodogram, clc, corrected_periodogram, corrected_folded_lightcurve, corrected_bin_folded_lc, t, corrected_bin_folded_lc, data_first_guess, fine_t, data_fit):

     tpf.plot(aperture_mask=tpf.pipeline_mask, mask_color='red')
     plt.title(EPIC+'_'+campaign_num+' Pipeline Mask')
     plt.savefig(EPIC+'_'+campaign_num+'PipelineMask')

     user_lc.plot(marker='o', linestyle='None', markersize=2, color='blue')
     plt.title(EPIC+'_'+campaign_num+' Raw Light Curve')
     plt.savefig(EPIC+'_'+campaign_num+'RawLightCurve')

     p.plot(c='k')
     plt.title(EPIC+'_'+campaign_num+' Raw Periodogram')
     plt.savefig(EPIC+'_'+campaign_num+'RawPeriodogram')
 
     raw_periodogram.plot()
     raw_periodogram.period_at_max_power
     plt.text(17,raw_periodogram.period_at_max_power.value,"Per=%s"%(raw_periodogram.period_at_max_power.value))
     plt.title(EPIC+'_'+campaign_num+' Raw Lomb Scargel Periodogram')
     plt.savefig(EPIC+'_'+campaign_um+'RawLombScargelPeriodogram')

     clc.plot()
     plt.title(EPIC+'_'+campaign_num+' Corrected Light Curve')
     plt.savefig(EPIC+'_'+campaign_num+'CorrectedLightCurve')

     corrected_periodogram.plot()
     corrected_periodogram.period_at_max_power
     plt.text(17,corrected_periodogram.period_at_max_power.value,"Per=%s"%(corrected_periodogram.period_at_max_power.value))
     plt.title(EPIC+'_'+campaign_num+' Raw Lomb Scargel Periodogram')
     plt.title(EPIC+'_'+campaign_num+' Corrected Lomb Scargel Periodogram')
     plt.savefig(EPIC+'_'+campaign_num+'CorrectedLombScargelPeriodogram')

     corrected_folded_lightcurve.plot(marker='o',linestyle='None')
     plt.text(17,corrected_periodogram.period_at_max_power.value,"Per=%s"%(corrected_periodogram.period_at_max_power.value))
     plt.title(EPIC+'_'+campaign_num+' Corrected Folded Light Curve')
     plt.savefig(EPIC+'_'+campaign_num+'CorrectedFoldedLightCurve')

     corrected_bin_folded_lc.plot(marker='o', linestyle='None', markersize=2, color='blue')
     plt.text(17,corrected_periodogram.period_at_max_power.value,"Per=%s"%(corrected_periodogram.period_at_max_power.value))
     plt.title(EPIC+'_'+campaign_num+' Corrected Binned Folded Light Curve')
     plt.savefig(EPIC+'_'+campaign_num+'CorrectedBinnedFoldedLightCurve')

     plt.plot(t, corrected_bin_folded_lc_.flux, marker='.', linestyle='none', color='black')
     plt.plot(t, data_first_guess, label='first guess', color='yellow')
     plt.plot(fine_t, data_fit_first, label='after fitting', color='red')
     plt.text(-0.525,1.045,"Per=%s"%(corrected_periodogram.period_at_max_power.value))
     plt.text(-0.525,1.04,"Amp=%s"%(est_amp_first))
     plt.title(EPIC+'_'+campaign_num+' Sine Fit Curve to Binned Folded Corrected Periodogram')
     plt.savefig(EPIC+'_'+campaign_num+'SineFit')


### SEPARATION OF LIGHTCURVE INTO SECTIONS ###

def sep_light_curve(clc):
    
     first_time=[]
     second_time=[]
     third_time=[]
     fourth_time=[]
     fifth_time=[]
     sixth_time=[]
     seventh_time=[]
     eighth_time=[]
     ninth_time=[]
     tenth_time=[]
     first_flux=[]
     second_flux=[]
     third_flux=[]
     fourth_flux=[]
     fifth_flux=[]
     sixth_flux=[]
     seventh_flux=[]
     eighth_flux=[]
     ninth_flux=[]
     tenth_flux=[]

     row=0
     first_day=clc.time[0]
     num_data=clc.time.size
     last_day=clc.time[num_data-1]
     num_days=last_day-first_day
     sep_length=num_days/10

     for i in user_lc.time:
         if (i<(first_day+sep_length)):
             first_time.append(i)
             first_flux.append(user_lc.flux[row])
         if (i>(first_day+sep_length)) and (i<(first_day+2*sep_length)):
             second_time.append(i)
             second_flux.append(user_lc.flux[row])
         if (i>(first_day+2*sep_length)) and (i<(first_day+3*sep_length)):
             third_time.append(i)
             third_flux.append(user_lc.flux[row])
         if (i>(first_day+3*sep_length)) and (i<(first_day+4*sep_length)):
             fourth_time.append(i)
             fourth_flux.append(user_lc.flux[row])
         if (i>(first_day+4*sep_length)) and (i<(first_day+5*sep_length)):
             fifth_time.append(i)
             fifth_flux.append(user_lc.flux[row])
         if (i>(first_day+5*sep_length)) and (i<(first_day+6*sep_length)):
             sixth_time.append(i)
             sixth_flux.append(user_lc.flux[row])
         if (i>(first_day+6*sep_length)) and (i<(first_day+7*sep_length)):
             seventh_time.append(i)
             seventh_flux.append(user_lc.flux[row])
         if (i>(first_day+7*sep_length)) and (i<(first_day+8*sep_length)):
             eighth_time.append(i)
             eighth_flux.append(user_lc.flux[row])
         if (i>(first_day+8*sep_length)) and (i<(first_day+9*sep_length)):
             ninth_time.append(i)
             ninth_flux.append(user_lc.flux[row])
         if (i>(first_day+9*sep_length)) and (i<(first_day+10*sep_length)):
             tenth_time.append(i)
             tenth_flux.append(user_lc.flux[row])
         row=row+1

     list_of_time_lengths=[len(first_time),len(second_time),len(third_time),len(fourth_time),len(fifth_time),len(sixth_time),len(seventh_time),len(eighth_time),len(ninth_time),len(tenth_time)]

### PRODUCING PHASE FOLDED BINNED LIGHT CURVES FOR EACH SEPARATION OF VARIABLE OBJECT ###

def separation_lightcruves(EPIC, campaign_num, clc, first_time, first_flux, second_time, second_flux, third_time, third_flux, fourth_time, fourth_flux, fifth_time, fifth_flux, sixth_time, sixth_flux, seventh_time, seventh_flux, eighth_time, eighth_flux, ninth_time, ninth_flux, tenth_time, tenth_flux):

     x=1
     y=list_of_time_lengths[0]
     for i in range (1,11):       
          root_path=s.getcwd()
          os.chdir(root_path+str(EPIC)+"_"+str(campaign_num))
          dir_name=root_path+str(EPIC)+"_"+str(campaign_num)+"SeparatedTime"+i
          os.mkdir(dir_name)
          os.chdir(dir_name)

          periodogram=lk.periodogram.LombScarglePeriodogram.from_lightcurve(clc[x:y], minimum_period=0.05, maximum_period =10)
          periodogram.plot()
          periodogram.period_at_max_power

          folded_lightcurve = clc[x:y].fold(periodogram.period_at_max_power.value)
          folded_lightcurve.plot(marker='o',linestyle='none')

          bin_folded_lc = folded_lightcurve_first.bin(5,method='median')
          bin_folded_lc.plot(marker='o',linestyle='None',markersize=4,color='blue')

          N = len(bin_folded_lc.flux) # number of data points
          t = np.linspace(-0.5, 0.5, N)
          guess_mean = np.mean(bin_folded_lc.flux)
          guess_std = 3*np.std(bin_folded_lc.flux)/(2**0.5)/(2**0.5)
          guess_phase = 0
          guess_freq = 1
          guess_amp = 1
          data_first_guess = guess_std*np.sin(t+guess_phase) + guess_mean
          optimize_func = lambda x: x[0]*np.sin(x[1]*t+x[2]) + x[3] - bin_folded_lc.flux
          est_amp, est_freq, est_phase, est_mean = leastsq(optimize_func, [guess_amp, guess_freq, guess_phase, guess_mean])[0]
          data_fit = est_amp*np.sin(est_freq*t+est_phase) + est_mean
          fine_t = np.arange(-0.5,0.5,0.001)
          data_fit=est_amp*np.sin(est_freq*fine_t+est_phase)+est_mean
          plt.plot(t, bin_folded_lc.flux, marker='.', linestyle='none', color='black')
          plt.plot(t, data_first_guess, label='first guess', color='yellow')
          plt.plot(fine_t, data_fit, label='after fitting', color='red')
          plt.text(-0.525,1.045,"Per=%s"%(periodogram.period_at_max_power.value))
          plt.text(-0.525,1.04,"Amp=%s"%(est_amp))

          x=x+y
          y=y+list_of_time_lengths[i]
      
          
