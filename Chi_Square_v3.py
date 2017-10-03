import pandas as pd
import numpy as np
import scipy
from scipy.stats import chisquare


############################################################################
def MultiProcessingPreparationByChromosomesOneInput(df):
    chrom_list = list(set(df['CHR']))
    array_chrom = []
    array_chrom_data = []
    for chrom in chrom_list:
        df_c = df[df['CHR'] == chrom]
        array_chrom_data.append(df_c)
        array_chrom.append(chrom)
    return array_chrom, array_chrom_data
############################################################################
def GroupByPeaks(chrom, df, expected_number_for_chi_square):
    print('GroupByPeaks: {0}'.format(chrom))
    print(df.head())
    if len(df) == 0:
        return df
    elif len(df) > 0:
        min_bin = df['Start'].min()
        max_bin = df['End'].min()
        dif_bin = int(max_bin) - min_bin
        number_bins = dif_bin/25
        #expected_read = expected_number_for_chi_square
         #########################################
        def ChiSquareTest(g, expected_read = expected_number_for_chi_square, number_bins = number_bins):
            array_expected = []
            for i in range(0,len(g)):
                array_expected.append(abs(expected_read))
            
            array_observed = g.tolist()
            if number_bins > len(array_observed):
                dif = number_bins - len(array_observed)
                for k in range(0,dif):
                    array_observed.append(0)

            #The first value in the returned tuple is the chi^2 value itself, 
            #while the second value is the p-value computed using nu=k-1 1 where k is the number of values in each array.
            observed_values=scipy.array(array_observed)
            expected_values=scipy.array(array_expected)
            ChisquareScore, Pvalues = chisquare(observed_values, f_exp=expected_values)
            return Pvalues
        #########################################
    
        df_n = df[df['Group'] > 0]
        
        col_array = df.columns
        dif_x = 'Diff_counts_x'
        dif_y = 'Diff_counts_y'
        for col_name in col_array:
            if col_name[:13] == 'Diff_counts_x':
                dif_x = col_name
            if col_name[:13] == 'Diff_counts_y':
                dif_y = col_name


        f = {'Start':{'Start':'min'},\
            'End': {'End':'max'},\
            dif_x:{dif_x:'sum'},\
            dif_y:{dif_y:'sum'},\
            'DeltaDelta':{'DeltaDelta':'sum'},\
            'Abs_DeltaDelta': {'Pvalues': ChiSquareTest}}
        test = df_n.groupby('Group').agg(f)
        df_out = pd.DataFrame(test)
        df_out.columns = df_out.columns.droplevel()
        df_out['CHR'] = chrom
        return df_out
############################################################################

def ProcessingArrayofZeroandOne(chrom, df):
    print('Processing Peaks at: {0}'.format(chrom))
    df.reset_index(inplace = True)
    start_array = list(df['Start'])
    end_array = list(df['End'])
    peak_number_array = []

    i = 0
    counter = 1
    for start in start_array:

        if i == 0:
            peak_number_array.append(counter)
        elif int(start) == int(end_array[i-1]):
            peak_number_array.append(counter)
        elif int(start) != int(end_array[i-1]):
            counter = counter + 1
            peak_number_array.append(counter)
            
        i = i + 1

    df1 = pd.DataFrame(peak_number_array)
    df1.columns = ['Group']

    df2 = pd.concat([df,df1], axis = 1)
    print('Grouping..')
    g = df2.groupby('Group')
    print('Filtering...')
    df3 = g.filter(lambda x: len(x) >3)
    df3.reset_index(drop = True, inplace = True)
    df3.drop(['index'], axis = 1, inplace = True)

    return df3
############################################################################
def ChisquearePreparationWorker(chrom, df):
    print(chrom)
    print('Making the Bins as 0, 1 and -1')
    number_bins = 8
    df.reset_index(drop = True, inplace = True)
    df['index'] = df.index
    diff_array = df['DeltaDelta']
    
    def MainWorker(input1, index1,  diff_array = diff_array):
        
        if input1 > 0:
            return 1
        elif input1 < 0:
            return -1
        elif input1 == 0 :
            back = index1 - number_bins
            front = index1 + number_bins
            ar_back = np.array(diff_array[back:index1])
            ar_front = np.array(diff_array[index1:front])

            zero_back = np.count_nonzero(ar_back)
            if zero_back >= number_bins-1:
                return 0
            zero_front = np.count_nonzero(ar_front)
            if zero_front >= number_bins - 1:
                return 0

            pos_back = (ar_back > 0).sum()
            neg_back = (ar_back < 0).sum()

            pos_front = (ar_front > 0).sum()
            neg_front = (ar_front < 0).sum()
            
            if pos_back > number_bins - 3 and pos_front > number_bins - 3 and neg_front == 0 and neg_back == 0:
                return 1
            elif neg_back > number_bins - 3 and neg_front > number_bins - 3 and pos_front == 0 and pos_back == 0:
                return -1
            elif pos_back > number_bins - 2 and pos_front > number_bins - 2:
                return 1
            elif neg_back > number_bins - 2 and neg_front > number_bins - 2:
                return -1

            ar_back_small = np.array(diff_array[(index1-3):index1])
            ar_front_small = np.array(diff_array[index1:(index1+3)])
            
            pos_back_small = (ar_back > 0).sum()
            pos_front_small = (ar_front > 0).sum()
            
            if pos_back_small == 3 and pos_front_small == 3:
                return 1

            neg_back_small = (ar_back < 0).sum()
            neg_front_small = (ar_front < 0).sum()
            
            if neg_back_small == 3 and neg_front_small == 3:
                return -1
            else:
                return 0

    df['Peak'] = np.vectorize(MainWorker)(df['DeltaDelta'], df['index'])

    print(df.head())
    df.drop(['index'], axis =1, inplace = True)
  
    return df

############################################################################






            











