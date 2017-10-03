import pandas as pd
import numpy as np
import re
import os
import glob
import multiprocessing
import logging
import argparse
import sys
sys.path.insert(0, "/starter/starter-02/ikagiamp/data/Filtered_BAM_files/")
from Peaks_Processing_v3 import *
from Chi_Square_v3 import *
######################################################################

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.FileHandler('K27_log')
handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
    # add the handlers to the logger
logger.addHandler(handler)
######################################################################
def parserMain():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c','--cpus', help = "Number of CPUs", type = int)  
    args = parser.parse_args()
        
    return args.cpus
######################################################################
def MakeDirectory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
######################################################################
def SaveDFtoPickle(directory, output_name, df):
    MakeDirectory(directory)
    path = os.path.join(directory, str(output_name))
    df.to_pickle(path)
######################################################################
def SaveDFtoCSV(directory,output_name, df):
    MakeDirectory(directory)
    path = os.path.join(directory, str(output_name))
    df.to_csv(path, sep=',', index = False)
######################################################################
def FilesRead(folder_name, file_name):
    path = os.path.join(folder_name, str(file_name))
    df = pd.read_csv(path, sep = '\t', header = None)
    df.columns = ['CHR','Start','End','Counts']
    return df
######################################################################
def DifferenceAminusB(counts1, counts2):

    c = float(counts1) - float(counts2)
    if c >= 0:
        return c
    elif c < 0:
        c = 0
        return c
    else:
        return 0

######################################################################
def SignalMinusContol(folder_name, signal,control):

    print('\nReading File: {0}\n'.format(signal))
    df_con1 = FilesRead(folder_name, signal)
    df_con1 = df_con1[df_con1['CHR'] != 'M']
    
    #TO BE REMOVED
    #df_con1 = df_con1[df_con1['CHR'] == 'chr5']
    #print(df_con1.head())
    #TO BE REMOVED-END

    print('\nReading File: {0}\n'.format(control))
    df_input1 = FilesRead(folder_name, control)
    df_input1 = df_input1[df_input1['CHR'] != 'M']

    #TO BE REMOVED
    #df_input1 = df_input1[df_input1['CHR'] == 'chr5']
    #print(df_input1.head())
    #TO BE REMOVED-END

    print('\nMerging....\n')
    df_con_full = pd.merge(df_con1,df_input1, on = ['CHR','Start','End'], how='outer')
    df_con_full.fillna(0, inplace = True)
    print(df_con_full.head())
    print('\nFile Name: {0}, length: {1}\n'.format(signal, len(df_con1)))
    print('\nFile Name: {0}, length: {1}\n'.format(control, len(df_input1)))
    print('\nMerged table length {0}.\nColumns Name: {1}\n'.\
        format(len(df_con_full),list(df_con_full.columns)))
    print('\nDifference 1: {0}'.format((len(df_con_full)-len(df_con1))))
    print('\nDifference 2: {0}'.format((len(df_con_full) - len(df_input1))))

    print('\nSubstructing....\n')
    df_con_full['Diff_counts'] = np.vectorize(DifferenceAminusB)(df_con_full['Counts_x'], \
        df_con_full['Counts_y'])

    #SaveDFtoCSV('Test_DAC_Mock_chrom5', signal, df_con_full)
    return df_con_full

######################################################################
def MergeDeltaDF(output1, output2, con1_name,con2_name):

    def DifferenceCalculator(n1,n2):
        c = float(n1) - float(n2)
        return float(c), abs(float(c))
        
    try:
        output1.drop(['Counts_x','Counts_y'], axis = 1, inplace = True)
        output2.drop(['Counts_x','Counts_y'], axis = 1, inplace = True)
    except:
        logger.info('Error at dropping')
        logger.info(con1_name,con2_name)
        logger.info(output1.head())
        logger.info(output2.head())
        pass
    left_suffix = '_x_'+con2_name
    right_suffix = '_y_'+con1_name
    df_delta_merged = pd.merge(output2,output1, on = ['CHR','Start','End'], suffixes=(left_suffix, right_suffix), how='outer')
    df_delta_merged.fillna(0, inplace = True)
    #SaveDFtoCSV('Test_DAC_Mock_chrom5', 'df_delta_merged.csv', df_delta_merged)
    print(df_delta_merged.head())
    print('\nDataFrame1 length: {0}\n'.format(len(output1)))
    print('\nDataFrame2 length: {0}\n'.format(len(output2)))
    print('\nDataFrame2 + DataFrame1 merge length: {0}\n'.format(len(df_delta_merged)))


    print('\nSubstructing.....\n')
    left_name = 'Diff_counts'+left_suffix
    right_name  = 'Diff_counts'+right_suffix
    df_delta_merged['DeltaDelta'], df_delta_merged['Abs_DeltaDelta'] =\
    np.vectorize(DifferenceCalculator)(df_delta_merged[left_name], df_delta_merged[right_name])

    print(df_delta_merged.head())
    
    return df_delta_merged
######################################################################
def FinalDataFrameReconstruction(df_array):
    df_final = pd.DataFrame()
    for df1 in df_array:
        df_final = pd.concat([df_final, df1])
    print('Final DataFrame:\n{0}\n'.format(df_final.head()))
    return df_final
######################################################################
def CallPeakImport(folder_name, file_name):
    path = os.path.join(folder_name, file_name)
    print('\nReading File: {0}\n'.format(file_name))
    df = pd.read_csv(path, sep = '\t', header = None)
    df1 = df[[0,1,2,4,9]]
    df1.columns = ['CHR','Start','End','Intensity','Top']
    #print(df1.head())
    return df1

######################################################################
def MultiProcessingPreparationsTwoInputs(df1, df2):
    chrom_list1 = list(set(df1['CHR']))
    #chrom_list2 = list(set(df2['CHR']))
    #chrom_list = list(set(chrom_list1+chrom_list2))

    df1_array = []
    df2_array = []
    df_chrom_array = []
    division_n = 0

    for chrom in chrom_list1:
        df1_chrom = df1[df1['CHR']==chrom]
        df2_chrom = df2[df2['CHR'] == chrom]
        print('\nThe number of DeltaDelta bins in {0} is {1}\n'.format(chrom, len(df1_chrom)))
        print('\nThe number of peaks at the {0} is {1}\n'.format(chrom, len(df2_chrom)))
        
        if len(df1_chrom) > 0:
            division_n = (len(df1_chrom)//1000)
        
        if len(df2_chrom) == 0 or len(df1_chrom) == 0:
            pass
        elif division_n == 0:
            if len(df1_chrom) > 0:
                df1_array.append(df1_chrom)
                df2_array.append(df2_chrom)
                df_chrom_array.append(chrom)
            #print(division_n, len(df_450_chrom), len(df_en_chrom), chrom)
        elif division_n > 0:
            df1_temp_list = np.array_split(df1_chrom, division_n)
            for df in df1_temp_list:
                df1_array.append(df)
                df2_array.append(df2_chrom)
                df_chrom_array.append(chrom)
                #print(division_n, len(df), len(df_en_chrom), chrom)
        else:
            print('Why here?')

    return df_chrom_array, df1_array, df2_array

######################################################################
def MultiprocessingModuleThreeInputs(function, df_chrom_array, df_array1, df_array2 ):
    #MultiprocessingModuleThreeInputs(deltadelta_array, peak_array, df_chrom_array)
    results = []
    #results_async = [MultiprocessingChrom(chrom, df_annot_c, df_en_c) for chrom, df_annot_c, df_en_c in zip(df_chrom_array, df_annot_array, df_en_array)]
    if len(df_array1) > 0 and len(df_array2) > 0 and len(df_chrom_array) > 0:
        results_async = [pool.apply_async(function, (chrom, df1_c, df2_c)) \
        for chrom, df1_c, df2_c in zip( df_chrom_array, df_array1, df_array2)]

        results=[r.get() for r in results_async]

    return results
######################################################################
def MultiprocessingModuleTwoInputs(function, df_chrom_array, df_array1):
    #MultiprocessingModuleThreeInputs(deltadelta_array, peak_array, df_chrom_array)
    results = []
    #results_async = [MultiprocessingChrom(chrom, df_annot_c, df_en_c) for chrom, df_annot_c, df_en_c in zip(df_chrom_array, df_annot_array, df_en_array)]
    if len(df_array1) > 0  and len(df_chrom_array) > 0:
        results_async = [pool.apply_async(function, (chrom, df1_c)) \
        for chrom, df1_c in zip( df_chrom_array, df_array1)]

        results=[r.get() for r in results_async]

    return results
######################################################################
def MultiprocessingModuleTwoInputsOneConstant(function, df_chrom_array, df_array1, OneConstant):
    #MultiprocessingModuleThreeInputs(deltadelta_array, peak_array, df_chrom_array)
    constant_array = []
    for i in range(0,len(df_chrom_array)):
        constant_array.append(OneConstant)

    results = []
    #results_async = [MultiprocessingChrom(chrom, df_annot_c, df_en_c) for chrom, df_annot_c, df_en_c in zip(df_chrom_array, df_annot_array, df_en_array)]
    if len(df_array1) > 0  and len(df_chrom_array) > 0:
        results_async = [pool.apply_async(function, (chrom, df1_c, con)) \
        for chrom, df1_c, con in zip( df_chrom_array, df_array1, constant_array)]

        results=[r.get() for r in results_async]

    return results
######################################################################
def SplitAndSaveBigDataFrame(df,cutoff,foldername,filename):
    array_df = []
    pickle_folder_name = 'pickle_'+ foldername
    if len(df) > int(cutoff):
        number = len(df)//int(cutoff)
        array_df = np.array_split(df, (number + 1))
        for index, df1 in enumerate(array_df):
            new_name = 'part'+str(index)+'_'+filename
            print('The length of the {0} DataFrame is:{1}'.format(new_name, len(df1)))
            print('\nSaving {0}'.format(new_name))
            logger.info('The length of the {0} DataFrame is:{1}'.format(new_name, len(df1)))
            SaveDFtoCSV(foldername,new_name, df1)
            #SaveDFtoPickle(pickle_folder_name,new_name, df1)
            print('Done Saving {0}'.format(new_name))
    elif len(df)<= int(cutoff):
        print('The length of the {0} DataFrame is:{1}'.format(filename, len(df)))
        print('\nSaving {0}'.format(filename))
        logger.info('The length of the {0} DataFrame is:{1}'.format(filename, len(df)))
        SaveDFtoCSV(foldername,filename, df)
        #SaveDFtoPickle(pickle_folder_name,filename, df)
        print('Done Saving {0}'.format(filename))
    else:
        logger.info('The length of the {0} DataFrame is:{1}'.format(filename, len(df)))

######################################################################
def MainWorkerPeakAnnotation(df_Pvalues,df_peaks, peak_name):
    df_peaks['Name'] = np.vectorize(PeakNameCreator)(df_peaks['CHR'],df_peaks['Start'], df_peaks['End'])
    
    df_chrom_array, deltadelta_array, peak_array = MultiProcessingPreparationsTwoInputs(df_Pvalues, df_peaks)
    print('\nAnnotating Peaks from condition: {0}'.format(peak_name))
    results_peaks = MultiprocessingModuleTwoInputsOneConstant(PeakAnnotation, deltadelta_array, peak_array, peak_name)
    df_with_peaks = FinalDataFrameReconstruction(results_peaks)
    print('\nAnnotation Finished')
    #SaveDFtoCSV('Test_peaks', 'peaks.csv', df_with_peaks)
    return df_with_peaks
######################################################################
def CleanSelection(con1, con2):
    if con1 != 'None' or con2 != 'None':
        return 1
    else:
        return 0
######################################################################
def MainStatistics(ac, peak_output1, peak_output2,con1_name, con2_name, output_file_name, foldername, expected_number_for_chi_square):
    
    
    df_chrom_array_for_Chisquare, df_array_data_for_Chisquare \
    = MultiProcessingPreparationByChromosomesOneInput(ac)
    df_Pvalues = MultiprocessingModuleTwoInputsOneConstant(GroupByPeaks, df_chrom_array_for_Chisquare, df_array_data_for_Chisquare, expected_number_for_chi_square)
    final_df_Pvalues = FinalDataFrameReconstruction(df_Pvalues)
    
    
    #SaveDFtoCSV('K27_intermediate_data',output_file_name, final_df_Pvalues)
    #path = 'K27_intermediate_data/'+ output_file_name
    #final_df_Pvalues = pd.read_csv(path, sep = ',')
    print(final_df_Pvalues.head())

    print('\nStart Annotating Peaks from condition: {0}'.format(con1_name))
    final_df_Pvalues_with_peak1 = MainWorkerPeakAnnotation(final_df_Pvalues,peak_output1, con1_name)
    print('\nStart Annotating Peaks from condition: {0}'.format(con2_name))
    final_df_Pvalues_with_peaks = MainWorkerPeakAnnotation(final_df_Pvalues_with_peak1,peak_output2, con2_name)

    print('\nCleaning Peaks...')
    final_df_Pvalues_with_peaks['Clean'] = np.vectorize(CleanSelection)(final_df_Pvalues_with_peaks[con1_name],final_df_Pvalues_with_peaks[con2_name])
    final_df_Pvalues_with_peaks_clean = final_df_Pvalues_with_peaks.loc[final_df_Pvalues_with_peaks['Clean'] == 1]
    print(final_df_Pvalues_with_peaks_clean.head())
    output_file_name1 = 'clean_'+ output_file_name

    cutoff = 2000000
    print('\nSaving...')
    SplitAndSaveBigDataFrame(final_df_Pvalues_with_peaks,cutoff,foldername,output_file_name)
    SplitAndSaveBigDataFrame(final_df_Pvalues_with_peaks_clean,cutoff,foldername,output_file_name1)
######################################################################
def MakingNoiseZero(noise, abs1, deltadelta):
    if abs1 < noise:
        abs1 = 0.0
    
    if deltadelta < 0:
        if deltadelta > -noise:
            deltadelta = 0.0
    elif deltadelta > 0:
        if deltadelta < noise:
            deltadelta = 0.0

    return abs1, deltadelta
       
######################################################################
def ContinuesPickAnnotation(df):
    print('Working on Peaks')
    df_chrom_array, df_array_data = MultiProcessingPreparationByChromosomesOneInput(df)
    array_dfs = MultiprocessingModuleTwoInputs(ProcessingArrayofZeroandOne, df_chrom_array, df_array_data)
    df1 = FinalDataFrameReconstruction(array_dfs)

    return df1
######################################################################
def AfterInputsProcessing(peak_output1, peak_output2, output1, output2,con1_name,con2_name, output_file_name, foldername):
    
    #Compining Peaks
    #To be removed
    #peak_output1 = peak_output1[peak_output1['CHR'] == 'chr5']
    #peak_output2 = peak_output2[peak_output2['CHR'] == 'chr5']
    #End To be removed

    peak_output1 = peak_output1[peak_output1['CHR'] != 'M']
    peak_output2 = peak_output2[peak_output2['CHR'] != 'M']
    
    print('\nSubstructing...\n')
    out2minusout1 = MergeDeltaDF(output1, output2, con1_name,con2_name)
    ##SaveDFtoCSV('Test_DAC_Mock_chrom5', 'out2minusout1.csv', out2minusout1)

    print('Calculating the median of all differences...')
    possitive_temp = out2minusout1[out2minusout1['Abs_DeltaDelta'] > 0]
    expected_number_for_chi_square = possitive_temp['Abs_DeltaDelta'].median()
    print('The Media is: {0}'.format(expected_number_for_chi_square))
    print('Remove Noise....')
    out2minusout1['Abs_DeltaDelta'], out2minusout1['DeltaDelta'] = np.vectorize(MakingNoiseZero)(expected_number_for_chi_square, out2minusout1['Abs_DeltaDelta'], out2minusout1['DeltaDelta'])
    
    print('\nStart the binning process...')
    df_chrom_array, df_array_data_out2minusout1 = MultiProcessingPreparationByChromosomesOneInput(out2minusout1)
    out2minusout1_temp = MultiprocessingModuleTwoInputs(ChisquearePreparationWorker, df_chrom_array, df_array_data_out2minusout1)
    df_with_bin_numbers = FinalDataFrameReconstruction(out2minusout1_temp)
    print(df_with_bin_numbers.head())
    
    #SaveDFtoCSV('Testing', 'df_with_bin_numbers.csv', df_with_bin_numbers)
    #SaveDFtoPickle('Testing', 'df_with_bin_numbers', df_with_bin_numbers)
    #df_with_bin_numbers = pd.read_pickle('Testing/df_with_bin_numbers')
    
    out2minusout1_pos = df_with_bin_numbers[df_with_bin_numbers['Peak'] > 0]
    out2minusout1_neg = df_with_bin_numbers[df_with_bin_numbers['Peak'] < 0]

    print('\nNumbering the neaby bins and removing the small number of aggregated bins...')
    out2minusout1_pos_n = ContinuesPickAnnotation(out2minusout1_pos)
    out2minusout1_neg_n = ContinuesPickAnnotation(out2minusout1_neg)
    
    #SaveDFtoCSV('Testing', 'out2minusout1_pos_n.csv', out2minusout1_pos_n)
    #SaveDFtoPickle('Testing', 'out2minusout1_pos_n', out2minusout1_pos_n)

    #out2minusout1_pos_n = pd.read_pickle('Testing/out2minusout1_pos_n')
    #print(out2minusout1_pos_n.head())
    
    #out2minusout1_pos_n = []
    #out2minusout1_neg_n = []
    #expected_number_for_chi_square = 1

    print('Start Statistics....')
    MainStatistics(out2minusout1_pos_n,peak_output1, peak_output2, con1_name,con2_name, output_file_name, foldername, expected_number_for_chi_square)
    deaccetylation_output = 'NegativePeaks_'+output_file_name
    MainStatistics(out2minusout1_neg_n,peak_output1, peak_output2, con1_name,con2_name, deaccetylation_output, foldername, expected_number_for_chi_square)

    
######################################################################
def ReadingInputsPeaks(folder_name_peakcall, fpeak_1_1):
    peak_output = CallPeakImport(folder_name_peakcall, fpeak_1_1)
    return peak_output
######################################################################
def ReadingInputBEDFiles(folder_name,file1_1,file2_1):
    #SignalMinusContol(folder_name, signal,control)
    output =  SignalMinusContol(folder_name,file1_1,file2_1)
    print(output.head())
    
    return output
######################################################################
if __name__ == '__main__':

    try:
        number_of_cpus = parserMain()
        print('\nNumber of CPUs: {0}\n'.format(number_of_cpus))
        logger.info('Number of CPUs: {0}'.format(number_of_cpus))
        pool = multiprocessing.Pool(number_of_cpus)
    except:
        pool = multiprocessing.Pool()

    ##############################################################
    folder_name_peakcall = '/starter/starter-02/ikagiamp/data/Filtered_BAM_files/MACS2_bdgpeakcall_files'

    fpeak_out_1 = 'Galaxy1743-[Mock_K27_Qnorm_bdgpeakcall_1500_300].txt'
    fpeak_out_2 = 'Galaxy1744-[455_K27_Qnorm_bdgpeakcall_1500_300].txt'
    fpeak_out_3 = 'Galaxy1745-[DAC_K27_Qnorm_bdgpeakcall_1500_300].txt'
    fpeak_out_4 = 'Galaxy1746-[DAC455_K27_Qnorm_bdgpeakcall_1500_300].txt'
    
    ##############################################################
    folder_name = '/starter/starter-02/ikagiamp/data/Filtered_BAM_files/MACS2_bdgdiff_bedgraph_input_files'

    #file1_1 = 'Galaxy12-[Mock_K27ac_Qnorm_Final].bedgraph'
    file1_2 = 'Galaxy9-[455_K27ac_Qnorm_Final].bedgraph'
    #file1_3 = 'Galaxy11-[DAC_K27ac_Qnorm_Final].bedgraph'
    file1_4 = 'Galaxy10-[DAC455_K27ac_Qnorm_Final].bedgraph'
    

    #file2_1 = 'Galaxy8-[Mock_INP_qNorm_Final].bedgraph'
    file2_2 = 'Galaxy5-[455_INP_qNorm_Final].bedgraph'
    #file2_3 = 'Galaxy7-[DAC_INP_qNorm_Final].bedgraph'
    file2_4 = 'Galaxy6-[DAC455_INP_qNorm_Final].bedgraph'

    ##############################################################

   
    #peak_output1 = ReadingInputsPeaks(folder_name_peakcall, fpeak_out_1)
    peak_output2 = ReadingInputsPeaks(folder_name_peakcall, fpeak_out_2)
    #peak_output3 = ReadingInputsPeaks(folder_name_peakcall, fpeak_out_3)
    peak_output4 = ReadingInputsPeaks(folder_name_peakcall, fpeak_out_4)
    #SaveDFtoPickle('outputK27', 'peak_output1', peak_output1)
    #SaveDFtoPickle('outputK27', 'peak_output2', peak_output2)
    #SaveDFtoPickle('outputK27', 'peak_output3', peak_output3)
    #SaveDFtoPickle('outputK27', 'peak_output4', peak_output4)

    #output1 = ReadingInputBEDFiles(folder_name,file1_1,file2_1)
    output2 = ReadingInputBEDFiles(folder_name,file1_2,file2_2)
    #output3 = ReadingInputBEDFiles(folder_name,file1_3,file2_3)
    output4 = ReadingInputBEDFiles(folder_name,file1_4,file2_4)
    #SaveDFtoPickle('outputK27', 'output1', output1)
    #SaveDFtoPickle('outputK27', 'output2', output2)
    #SaveDFtoPickle('outputK27', 'output3', output3)
    #SaveDFtoPickle('outputK27', 'output4', output4)

    #peak_output1 = pd.read_pickle('outputK27/peak_output1')
    #peak_output2 = pd.read_pickle('outputK27/peak_output2')
    #peak_output3 = pd.read_pickle('outputK27/peak_output3')
    #peak_output4 = pd.read_pickle('outputK27/peak_output4')
    #output1 = pd.read_pickle('outputK27/output1')
    #output2 = pd.read_pickle('outputK27/output2')
    #output3 = pd.read_pickle('outputK27/output3')
    #output4 = pd.read_pickle('outputK27/output4')
    #output1 = []
    #output2 = []
    #output3 = []
    #output4 = []

    out_path = '/starter/starter-02/ikagiamp/data/Filtered_BAM_files/'
    #print('Start Processing....')
    #logger.info('Start 1')
    #print('Start_1')
    #AfterInputsProcessing(peak_output1, peak_output2, output1, output2,'mock_K27','455_K27', '455_K27_vs_mock_K27.csv', out_path+'455_K27_vs_mock_K27')
    #logger.info('Start 2')
    #print('Start_2')
    #AfterInputsProcessing(peak_output1, peak_output3, output1, output3,'mock_K27','DAC_K27', 'DAC_K27_vs_mock_K27.csv', out_path+'DAC_K27_vs_mock_K27')
    #logger.info('Start 3')
    #print('Start_3')
    #AfterInputsProcessing(peak_output1, peak_output4, output1, output4,'mock_K27', 'DAC455_K27', 'DAC455_K27_vs_mock_K27.csv', out_path+'DAC455_K27_vs_mock_K27')
    #logger.info('Start 4')
    #print('Start_4')
    #AfterInputsProcessing(peak_output3, peak_output4, output3, output4,'DACK27','DAC455_K27', 'DAC455_K27_vs_DACK27.csv', out_path+'DAC455_K27_vs_DACK27')
    #logger.info('Start 5')
    print('Start_5')
    AfterInputsProcessing(peak_output2, peak_output4, output2, output4,'455K27', 'DAC455_K27', 'DAC455_K27_vs_455K27.csv', out_path+'DAC455_K27_vs_455K27')
    logger.info('Finished')
    print('Finished')


   
