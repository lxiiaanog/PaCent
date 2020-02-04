#!/usr/bin/env python
# -*- coding:utf-8 -*-
'''
Created by Jian Song on 20200203
Modified by Xiao Liang on 20200203
'''
import numpy as np
import pandas as pd
import argparse

def main_fun(fpath, rt_shift_tol, out_path):
    # read OpenSWATH result
    df_osw = pd.read_csv(fpath, sep=None, engine='python')

    # order via peptide seqences
    df_osw['public_seq'] = df_osw['FullPeptideName'].str.replace('(UniMod:188)','', regex=False)
    df_osw.loc[:, 'L_H'] = 'L'
    df_osw.loc[df_osw.FullPeptideName.str.contains('188'), 'L_H'] = 'H'
    df_osw['pg_id'] = df_osw.FullPeptideName + '_' + df_osw.Charge.astype(str) + '_' + df_osw.RT.round(2).astype(str)

    df_osw.sort_values(by=['Charge', 'public_seq', 'L_H'],
                       ascending=[True, True, False],
                       inplace=True)
    df_osw.reset_index(inplace=True, drop=True)

    def getPairFullName(name):
        if name.find('188)') != -1:
            return name.replace('R(UniMod:188)', 'R').replace('K(UniMod:188)', 'K') # L for H
        else:
            return name.replace('R', 'R(UniMod:188)').replace('K', 'K(UniMod:188)') # H for L

    # pair-centric peptide process
    result_rt_shift = []
    result_pair_idx = []
    result_pair_id = []
    result_pair_LH_ratio = []

    row_idx = 0
    df_osw_rows = len(df_osw)

    start_print_idx = 0
    while row_idx < df_osw_rows:
        if (row_idx == 0) or ((row_idx - start_print_idx) // 10000 == 1):
            start_print_idx = row_idx
            print('{:8d}/{:d}'.format(row_idx//10000*10000, df_osw_rows))

        # get L/H index of seqences
        bias = 0
        seq_L = df_osw.iloc[row_idx]['FullPeptideName']

        seq_H = getPairFullName(seq_L)
        charge = df_osw.iloc[row_idx]['Charge']
        bias_L, bias_H = [], []
        while True:
            if row_idx+bias < df_osw_rows:
                row = df_osw.iloc[row_idx+bias]
                if (row['FullPeptideName'] == seq_L) and (row['Charge'] == charge):
                    bias_L.append(bias)
                    bias += 1
                elif (row['FullPeptideName'] == seq_H) and (row['Charge'] == charge):
                    bias_H.append(bias)
                    bias += 1
                else: break
            else: break

        # for singular peptides, complement with nan
        if len(bias_H) == 0 or len(bias_L) == 0:

            result_rt_shift.extend(np.array([np.nan] * (len(bias_L) + len(bias_H))))
            result_pair_idx.extend(np.array([np.nan] * (len(bias_L) + len(bias_H))))
            result_pair_id.extend(np.array([np.nan] * (len(bias_L) + len(bias_H))))
            result_pair_LH_ratio.extend(np.array([np.nan] * (len(bias_L) + len(bias_H))))
            row_idx = row_idx + len(bias_H) + len(bias_L)
            continue

        # calculate rt_shift of peptide pairs
        matrix_rts = np.zeros((len(bias_L), len(bias_H)))

        for i in bias_L:

            i_row = df_osw.iloc[row_idx+i]
            i_rts = float(i_row['RT'])

            for j in bias_H:
                j_row = df_osw.iloc[row_idx+j]
                j_rts = float(j_row['RT'])
                matrix_rts[i][j-len(bias_L)] = np.abs(i_rts-j_rts)

        # generate 2D matrix with rows annotating light and columns annotating H
        for i in bias_L:
            j = len(bias_L) + np.argmin(matrix_rts[i])
            L_intensity = df_osw.loc[row_idx + i, 'Intensity']
            H_intensity = df_osw.loc[row_idx + j, 'Intensity']
            L_rt = '%.2f' % df_osw.loc[row_idx + i, 'RT']
            H_rt = '%.2f' % df_osw.loc[row_idx + j, 'RT']
            public_seq = df_osw.loc[row_idx + i, 'public_seq']
            charge = str(df_osw.loc[row_idx + i, 'Charge'])
            pair_id = public_seq + '_' + charge + '_' + L_rt + '_' + H_rt

            result_pair_id.append(pair_id)
            result_pair_LH_ratio.append(L_intensity/H_intensity)
            result_rt_shift.append(np.min(matrix_rts[i]))
            result_pair_idx.append(row_idx+j)

        for j in bias_H:
            i = np.argmin(matrix_rts[:,j-len(bias_L)])
            L_intensity = df_osw.loc[row_idx + i, 'Intensity']
            H_intensity = df_osw.loc[row_idx + j, 'Intensity']
            L_rt = '%.2f' % df_osw.loc[row_idx + i, 'RT']
            H_rt = '%.2f' % df_osw.loc[row_idx + j, 'RT']
            public_seq = df_osw.loc[row_idx + i, 'public_seq']
            charge = str(df_osw.loc[row_idx + i, 'Charge'])
            pair_id = public_seq + '_' + charge + '_' + L_rt + '_' + H_rt

            result_pair_id.append(pair_id)
            result_pair_LH_ratio.append(L_intensity / H_intensity)
            result_rt_shift.append(np.min(matrix_rts[:, j-len(bias_L)]))
            result_pair_idx.append(row_idx + i)

        row_idx = row_idx + len(bias_H) + len(bias_L)

    df_osw['pair_idx'] = result_pair_idx
    df_osw['pair_rt_shift'] = result_rt_shift
    df_osw['pair_id'] = result_pair_id
    df_osw['pair_LH_ratio'] = result_pair_LH_ratio

    # delete unpaired singulars
    df_osw_without_na = df_osw.loc[~np.isnan(df_osw.pair_idx)]

    # delete pairs with over tolerance rt_shift 
    df_filter = df_osw_without_na.loc[df_osw_without_na.pair_rt_shift < rt_shift_tol, :].copy()

    # combine discrimination scores
    columns = df_osw.columns.values
    var_idx_start, var_idx_end = (np.where(df_osw.columns.str.startswith('var'))[0][0],
                                  np.where(df_osw.columns.str.startswith('var'))[0][-1])
    df_scores_all = pd.DataFrame(np.empty((len(df_filter), (var_idx_end - var_idx_start + 1) * 2)))
    df_scores_all[:] = np.nan
    df_scores_all.columns = list(columns[var_idx_start:(var_idx_end + 1)] + '_L') + \
                            list(columns[var_idx_start:(var_idx_end + 1):] + '_H')

    df_L = df_filter.loc[df_filter.L_H == 'L', :]
    df_H = df_filter.loc[df_filter.L_H == 'H', :]

    df_scores_all.iloc[np.where(df_filter.L_H == 'L')[0], 0:(var_idx_end - var_idx_start + 1)] = \
        df_filter.loc[df_filter.L_H == 'L', columns[var_idx_start:(var_idx_end+1)]].values
    df_scores_all.iloc[np.where(df_filter.L_H == 'H')[0], (var_idx_end - var_idx_start + 1):] = \
        df_filter.loc[df_filter.L_H == 'H', columns[var_idx_start:(var_idx_end+1)]].values
    df_scores_all.iloc[np.where(df_filter.L_H == 'H')[0], 0:(var_idx_end - var_idx_start + 1)] = \
        df_osw.loc[df_H.pair_idx.astype(int), columns[var_idx_start:(var_idx_end+1)]].values
    df_scores_all.iloc[np.where(df_filter.L_H == 'L')[0], (var_idx_end - var_idx_start + 1):] = \
        df_osw.loc[df_L.pair_idx, columns[var_idx_start:(var_idx_end+1)]].values

    # average "main" score for combined pairs 
    df_filter.loc[:, 'main_var_xx_swath_prelim_score'] = \
        (df_filter.main_var_xx_swath_prelim_score.values + df_osw.loc[df_filter.pair_idx, 'main_var_xx_swath_prelim_score'].values)/2

    # clean remains of previous information
    df_filter.drop(columns[var_idx_start:(var_idx_end+1)], axis=1, inplace=True)
    df_filter.reset_index(drop=True, inplace=True)
    df_scores_all.reset_index(drop=True, inplace=True)
    result = pd.concat([df_filter, df_scores_all], axis=1)

    # output file
    print('{:8d}/{:d}'.format(df_osw_rows, df_osw_rows))
    result.to_csv(out_path, sep='\t', index=False)



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description= 'PaCent_Pairing - Combine light/heavy paired peptides.')
    parser.add_argument('-rt_shift_tolerance', type=float,
                        help='maximum rt_shift time (seconds) of light/heavy pairs, example:5', required=True)
    parser.add_argument('-input',
                        help='tab separated openswath result, example:openswath.tsv', required=True)
    parser.add_argument('-output', 
                        help='tab separated pacent_pairing result, example:pacent_pairing.tsv', required=True)

    args = parser.parse_args()

    main_fun(args.input, args.rt_shift_tolerance, args.output)

    print('Done.')

