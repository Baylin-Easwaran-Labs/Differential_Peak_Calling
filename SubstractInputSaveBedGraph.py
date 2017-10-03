import pandas as pd
import numpy as np
import re
import os
import glob

######################################################################
def MakeDirectory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
######################################################################
def SaveDFtoCSV(directory,output_name, df):
    MakeDirectory(directory)
    path = os.path.join(directory, str(output_name))
    df.to_csv(path, sep='\t', index = False, header = False)
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

    df_small = df_con_full[['CHR','Start','End','Diff_counts']]
    df_small['Start'] = df_small['Start'].astype(int)
    df_small['End'] = df_small['End'].astype(int)

    signal_name = 'Clean_'+signal
    SaveDFtoCSV('Result_Homer',  signal_name, df_small)
    return df_small

######################################################################
def ReadingInputBEDFiles(folder_name,file1_1,file2_1):
    #SignalMinusContol(folder_name, signal,control)
    output =  SignalMinusContol(folder_name,file1_1,file2_1)
    print(output.head())
    
    return output
######################################################################
if __name__ == '__main__':

  
    ##############################################################
    folder_name = 'MACS2_bdgdiff_bedgraph_input_files'

    file1_1 = 'Galaxy12-[Mock_K27ac_Qnorm_Final].bedgraph'
    file1_2 = 'Galaxy1280-[455_K27ac_Qnorm_Final].bedgraph'
    file1_3 = 'Galaxy11-[DAC_K27ac_Qnorm_Final].bedgraph'
    file1_4 = 'Galaxy1281-[DAC455_K27ac_Qnorm_Final].bedgraph'
    

    file2_1 = 'Galaxy8-[Mock_INP_qNorm_Final].bedgraph'
    file2_2 = 'Galaxy5-[455_INP_qNorm_Final].bedgraph'
    file2_3 = 'Galaxy7-[DAC_INP_qNorm_Final].bedgraph'
    file2_4 = 'Galaxy6-[DAC455_INP_qNorm_Final].bedgraph'

    ##############################################################

    output1 = ReadingInputBEDFiles(folder_name,file1_1,file2_1)
    print(output1.head())
    output2 = ReadingInputBEDFiles(folder_name,file1_2,file2_2)
    print(output2.head())
    output3 = ReadingInputBEDFiles(folder_name,file1_3,file2_3)
    print(output3.head())
    output4 = ReadingInputBEDFiles(folder_name,file1_4,file2_4)
    print(output4.head())
