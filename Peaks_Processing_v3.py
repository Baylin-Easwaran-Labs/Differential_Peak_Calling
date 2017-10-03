import pandas as pd
import numpy as np
import re
import os
import glob
from math import log
import multiprocessing
#import Differential_Peak_Calling
#from K4_Final_Differential_Peak_Calling_v3 import MultiprocessingModuleTwoInputsOneConstant, MultiProcessingPreparationsTwoInputs, FinalDataFrameReconstruction

#####################################################################################################
def DesisionMaking(start_peak, end_peak, start_bin, end_bin):
    size = 1000
    if (int(start_peak) - size) <= int(start_bin) and int(start_bin) <= (int(end_peak) + size):
        return 1
    elif int(end_bin) >= (int(start_peak) - size) and int(end_bin) <= (int(end_peak) + size):
        return 1
    elif int(start_bin) <= int(start_peak) and int(end_bin) >= int(end_peak):
        return 1
    else:
        return 0
#####################################################################################################
def PeakAnnotation(delta_delta_df, peak_df, peak_name):
    
    ######################################################################
    def SecondRound(start_bin, end_bin, df_peak_c = peak_df):
        df_peak_c['Pos'] = np.vectorize(DesisionMaking)(df_peak_c['Start'], df_peak_c['End'], start_bin, end_bin)
        df_peak_c_pos = df_peak_c[df_peak_c['Pos'] == 1]
        if len(df_peak_c_pos) > 1:
            array_t = list(df_peak_c_pos['Name'])
            myString = ",".join(array_t)
            return myString
        elif len(df_peak_c_pos) == 1:
            array_t = list(df_peak_c_pos['Name'])
            return array_t[0]
        else:
            return 'None'
    ######################################################################

    delta_delta_df[peak_name] = np.vectorize(SecondRound)(delta_delta_df['Start'],delta_delta_df['End'])
    return delta_delta_df
#####################################################################################################
def PeakNameCreator(chrom, start, end):
    string = chrom+'_'+str(start)+'_'+str(end)

    return string
#####################################################################################################


    


