# Allison McCarthy
# Summer 2019
# This script calls all functions to do analysis for all EPICs in a CSV. Individual EPICs can also be called.


from functions import correct_lightcurve
from functions import fit_sin_curve
from functions import plot_all
from functions import sep_light_curve
from functions import separation_lightcurves
import csv

#Eventually allow user to select CSV or individual EPIC
path_to_csv_file = input('Please type path to csv file: \n')
f=open(path_to_csv_file)
csv_f=csv.reader(f)
header=next(csv_f)
for row in csv_f:
  try:
     EPIC = int(row[0])
     ob_name = row[1]
     campaign_num = int(row[2])
     variability=str(row[6])

     if variability=='Yes':

          tpf, user_lc, p, raw_periodogram, period_at_max_power_user_lc, lc, clc, corrected_periodogram, corrected_periodogram_at_max_power, corrected_folded_lightcurve, corrected_bin_folded_lc, corrector, new_lc, SFF_corrected_periodogram, SFF_corrected_period, SFF_corrected_folded_lightcurve, SFF_bin_folded_lc = correct_lightcurve(EPIC,campaign_num,ob_name)

          N, t, guess_mean, guess_std, guess_phase, guess_freq, guess_amp, data_first_guess, optimize_func, est_amp, est_freq, est_phase, est_mean, data_fit, fine_t = fit_sin_curve(corrected_bin_folded_lc)

          plot_all(EPIC, campaign_num, tpf, user_lc, p, raw_periodogram, clc, corrected_periodogram, corrected_folded_lightcurve, t, corrected_bin_folded_lc, data_first_guess, fine_t, data_fit, est_amp)

          list_of_time_lengths, first_time, first_flux, second_time, second_flux, third_time, third_flux, fourth_time, fourth_flux, fifth_time, fifth_flux, sixth_time, sixth_flux, seventh_time, seventh_flux, eighth_time, eighth_flux, ninth_time, ninth_flux, tenth_time, tenth_flux = sep_light_curve(clc)

          separation_lightcurves(EPIC, ob_name, campaign_num, clc, list_of_time_lengths)


  except:
     print(str(EPIC)+"_"+str(campaign_num)+' has a problem')
     pass 
