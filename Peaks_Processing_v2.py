import pandas as pd
import numpy as np
import re
import os
import glob
from math import log
import multiprocessing
#import Differential_Peak_Calling

#########################################################################################
def PeakFilesCompinations(peak1, peak2):
    df_p = pd.concat([peak1,peak2])
    df_p.sort_values(by = ['Start'], ascending=[True], inplace=True)
    return df_p
#########################################################################################
def MakingDF(array1, array2,array3, array4):
    df1 = pd.DataFrame(array1)
    df2 = pd.DataFrame(array2)
    df3 = pd.DataFrame(array3)
    df4 = pd.DataFrame(array4)
    
    df = pd.concat([df1,df2,df3,df4], axis = 1)
    return df
#########################################################################################
def DimessionReductionPeak(condition, chrom, df, array_column_names, array_number_rows):
    bin_size = condition

    start_temp = 0
    end_temp = 0
    counter = 0
    number_add = 0
    array_chrom = []
    array_start = []
    array_end = []
    array_number = []
    
    for index, row in df.iterrows():
        #print (row[1], row[2])
        start = row[array_number_rows[0]]
        end = row[array_number_rows[1]]
        number = row[array_number_rows[2]]

        if end < start:
            print ('Error, The peak end is smaller than the start')
            exit()

        if counter == 0:
            start_temp = int(start) - bin_size
            end_temp = int(end)
            number_add = float(number)
            counter = 1

        elif (int(end_temp) + bin_size) >= int(start):
            end_temp = int(end)
            number_add = number_add + float(number)
        
        elif (int(end_temp) + bin_size) < int(start):
            array_chrom.append(chrom)
            array_start.append(start_temp)
            array_end.append(end_temp)
            #counts_add = float(counts_add) + float(ads_delta_counts)
            if number_add == 0:
                number_add = float(number)
            array_number.append(number_add)
            number_add = float(number)
            #intesity_add = 0
            start_temp = int(start) - bin_size
            end_temp = int(end)

        else:
            print('Never go here')
            
    array_chrom.append(chrom)
    array_start.append(start)
    array_end.append(end)
    number_add = number_add +  float(number)
    array_number.append(number_add)
    df_f = MakingDF(array_chrom, array_start, array_end, array_number)
    df_f.columns = array_column_names
    return df_f
#########################################################################################
def HigherWorkerOneForPeaks(chrom, df_peak_1, df_peak_2):
    results = []
    df1 = PeakFilesCompinations(df_peak_1, df_peak_2)
    df1['length'] = df1['End'] - df1['Start']
    number_of_peaks = len(df1)
    total_length = df1['length'].sum(axis=0)
    total_counts = df1['Intensity'].sum(axis = 0)
    results.append(chrom)
    results.append(number_of_peaks)
    results.append(total_length)
    results.append(total_counts)
    
    df_t1 = pd.DataFrame(results)
    df_t2 = df_t1.T
    df_t2.columns = ['CHR','Number_of_Peaks','Total_Length', 'Total_Count']

    return df_t2
#########################################################################################
def DimessionReductionPeakWorker(chrom, df_peak_1, df_peak_2):
    df = PeakFilesCompinations(df_peak_1, df_peak_2)
    array_column_names = ['CHR','Start', 'End','Intensity']
    array_number_rows = [1,2,3]
    df1 = DimessionReductionPeak(2000, chrom, df, array_column_names, array_number_rows)
    return df1
#########################################################################################
def MultiProcessingPreparationByChromosomesTwoInputs(df1, df2):
    chrom_list1 = list(set(df1['CHR']))
    chrom_list2 = list(set(df2['CHR']))
    chrom_list = list(set(chrom_list1+chrom_list2))

    df1_array = []
    df2_array = []
    df_chrom_array = []
   
    for chrom in chrom_list:
        df1_chrom = df1[df1['CHR']==chrom]
        df2_chrom = df2[df2['CHR'] == chrom]
        print('\nThe number of peaks in {0} is {1}\n'.format(chrom, len(df1_chrom)))
        print('\nThe number of peaks at the {0} is {1}\n'.format(chrom, len(df2_chrom)))
        
        if len(df1_chrom) > 0 and len(df2_chrom) > 0:
            df1_array.append(df1_chrom)
            df2_array.append(df2_chrom)
            df_chrom_array.append(chrom)
 
    return df_chrom_array, df1_array, df2_array

######################################################################
def AveragingPeaksReads(df):
    average_read_per_bean = df['DeltaDelta'].mean()
    df['Expected_reads'] = average_read_per_bean
    return df




    


