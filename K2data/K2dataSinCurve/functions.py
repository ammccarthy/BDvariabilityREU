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

def correct_lightcurve(csv=NONE, EPIC=NONE, campaign_num=NONE):
    if csv is not NONE:
        f=open(csv)
        csv_f=csv.reader(f)
        header=next(csv_f)
        for row in csv_f:
            try:
                EPIC = int(row[0])
                ob_name = row[1]
                campaign_num = int(row[2])
                root_path=s.getcwd()
                dir_path=root_path+str(EPIC)+"_"+str(campaign_num)+"_"+ob_name
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

     counter=0
     first_day=clc.time[0]
     num_data=clc.time.size
     last_day=clc.time[num_data-1]
     num_days=last_day-first_day
     sep_length=num_days/10

     for i in user_lc.time:
         if (i<(first_day+sep_length)):
             first_time.append(i)
             row=counter
             first_flux.append(user_lc.flux[row])
         if (i>(first_day+sep_length)) and (i<(first_day+2*sep_length)):
             second_time.append(i)
             row=counter
             second_flux.append(user_lc.flux[row])
         if (i>(first_day+2*sep_length)) and (i<(first_day+3*sep_length)):
             third_time.append(i)
             row=counter
             third_flux.append(user_lc.flux[row])
         if (i>(first_day+3*sep_length)) and (i<(first_day+4*sep_length)):
             fourth_time.append(i)
             row=counter
             fourth_flux.append(user_lc.flux[row])
         if (i>(first_day+4*sep_length)) and (i<(first_day+5*sep_length)):
             fifth_time.append(i)
             row=counter
             fifth_flux.append(user_lc.flux[row])
         if (i>(first_day+5*sep_length)) and (i<(first_day+6*sep_length)):
             sixth_time.append(i)
             row=counter
             sixth_flux.append(user_lc.flux[row])
         if (i>(first_day+6*sep_length)) and (i<(first_day+7*sep_length)):
             seventh_time.append(i)
             row=counter
             seventh_flux.append(user_lc.flux[row])
         if (i>(first_day+7*sep_length)) and (i<(first_day+8*sep_length)):
             eighth_time.append(i)
             row=counter
             eighth_flux.append(user_lc.flux[row])
         if (i>(first_day+8*sep_length)) and (i<(first_day+9*sep_length)):
             ninth_time.append(i)
             row=counter
             ninth_flux.append(user_lc.flux[row])
         if (i>(first_day+9*sep_length)) and (i<(first_day+10*sep_length)):
             tenth_time.append(i)
             row=counter
             tenth_flux.append(user_lc.flux[row])
         counter=counter+1

