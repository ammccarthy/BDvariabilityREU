# Allison McCarthy
# Summer 2019
# This code will run through a list of EPICs and create a lightcurve for each. 

import os
import csv
import lightkurve as lk
from lightkurve import search_targetpixelfile
import matplotlib.pyplot as plt
import astropy.units as u
import sys, string, calendar, datetime, traceback

#Need to create a list maybe to read from with all EPICs, maybe read from file? Prompt user?
path_to_csv_file = input('Please type path to csv file: \n')
f=open(path_to_csv_file)
csv_f=csv.reader(f)
header=next(csv_f)
for row in csv_f:
  try:
     EPIC = int(row[0])
     ob_name = row[1]
     campaign_num = int(row[2])
     root_path="/Users/AllieMcCarthy/REU/K2data/K2dataJohannasNotes/"
     dir_path=root_path+str(EPIC)+"_"+str(campaign_num)+"_"+ob_name
     os.mkdir(dir_path)
     os.chdir(dir_path)

#Code from jupyter notebook

     #Aperture Mask
     tpf = search_targetpixelfile(EPIC, mission='K2', campaign=campaign_num, cadence='short').download()
     tpf.plot(aperture_mask=tpf.pipeline_mask, mask_color='red')
     plt.title(ob_name+' Pipeline Mask')
     plt.savefig(ob_name+'PipelineMask')

     #first plot
     user_lc = tpf.to_lightcurve(aperture_mask=tpf.pipeline_mask.astype(bool))
     # Clean the light curve
     user_lc = user_lc.remove_nans().remove_outliers()
     user_lc.plot(marker='o',linestyle='None',markersize=4,color='blue')
     #plt.ylim(0.9,1.1)
     plt.title(ob_name+' Raw LC')
     plt.savefig(ob_name+'RAWLC.png')

     p = user_lc.to_periodogram(freq_unit=u.microHertz, maximum_frequency=400, minimum_frequency=10)
     ax = p.plot(c='k');
     plt.title(ob_name+' Raw Periodogram')
     plt.savefig(ob_name+'RawPeriodogram.png')
     
     periodogram=lk.periodogram.LombScarglePeriodogram.from_lightcurve(user_lc, minimum_period=0.05, maximum_period =10)
     periodogram.plot()
     periodogram.period_at_max_power
     plt.title(ob_name+'Lomb Scargel Periodogram')
     plt.savefig(ob_name+'LombScarglePeriodogram.png')

     lc = tpf.to_lightcurve().normalize().remove_nans().remove_outliers()
     clc = lc.correct(windows=10).remove_outliers().fill_gaps()
     clc.plot() 
     plt.title(ob_name+' Corrected Light Curve')
     plt.savefig(ob_name+'CorrectedLightCurve.png')

     periodogram=lk.periodogram.LombScarglePeriodogram.from_lightcurve(clc, minimum_period=0.05, maximum_period =100)
     periodogram.plot()
     periodogram.period_at_max_power
     plt.title(ob_name+'Lomb Scargel Corrected Periodogram')
     plt.savefig(ob_name+'LombScargelCorrectedPeriodogram.png')    

     folded_lightcurve = clc.fold(periodogram.period_at_max_power.value)
     folded_lightcurve.plot(marker='o',linestyle='none')
     plt.title(ob_name+' Folded Light Curve')
     plt.savefig(ob_name+'FoldedLightCurve.png')

     bin_folded_lc = folded_lightcurve.bin(50,method='median')
     bin_folded_lc.plot(marker='o',linestyle='None',markersize=4,color='blue')
     plt.title(ob_name+' Bin Folded Light Curve')
     plt.savefig(ob_name+'BinFoldedLightCurve.png')

     corrector=lk.SFFCorrector(lc)
     new_lc = corrector.correct(lc.centroid_col,lc.centroid_row)

     periodogram=lk.periodogram.LombScarglePeriodogram.from_lightcurve(new_lc, minimum_period=0.05, maximum_period =100)
     periodogram.plot()
     periodogram.period_at_max_power
     plt.title(ob_name+' SFF Corrected Lomb Scargle Periodogram')
     plt.savefig(ob_name+'SFFCorrectedLombScarglePeriodogram')

     folded_lightcurve = new_lc.fold(periodogram.period_at_max_power.value)
     bin_folded_lc = folded_lightcurve.bin(50,method='median')
     bin_folded_lc.plot(marker='o',linestyle='None',markersize=4,color='blue')
     plt.title(ob_name+' SFF Corrected Bin Folded Light Curve')
     plt.savefig(ob_name+'SFFCorrectedBinFoldedLightCurve')

     plt.close("all")
     os.chdir("/Users/AllieMcCarthy/REU/K2data")

  except:
     print(str(EPIC)+str(campaign_num)+' has a problem')
     os.rename(dir_path,dir_path+'HasAnError')
     pass
