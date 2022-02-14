#!/usr/bin/env python
# encoding:utf-8
"""
Author : Yuanqing Mei
Date : 2022/1/7
Time: 19:40
File: 5_auc_precision_recall_testing.py
HomePage : https://github.com/meiyuanqing
Email : dg1533019@smail.nju.edu.cn

Performance metric: AUC,precision,recall
Thresholds derived from five supervised and their meta-analysis will be tested on 66 versions.

"""


def predictive_performance(metric, metric_p, metric_t, df_test):
    from sklearn.metrics import recall_score, precision_score, f1_score, roc_curve, auc, roc_auc_score, confusion_matrix

    if metric_p > 0:
        df_test['predictBinary'] = df_test[metric].apply(lambda x: 1 if float(x) >= metric_t else 0)
    else:
        df_test['predictBinary'] = df_test[metric].apply(lambda x: 1 if float(x) <= metric_t else 0)

    # confusion_matrix()函数中需要给出label, 0和1，否则该函数算不出TP,因为不知道哪个标签是poistive.
    c_matrix = confusion_matrix(df_test["bugBinary"], df_test['predictBinary'], labels=[0, 1])
    tn, fp, fn, tp = c_matrix.ravel()

    if (tn + fp) == 0:
        tnr_value = 0
    else:
        tnr_value = tn / (tn + fp)

    if (fp + tn) == 0:
        fpr = 0
    else:
        fpr = fp / (fp + tn)

    auc_value = roc_auc_score(df_test['bugBinary'], df_test['predictBinary'], labels=[0, 1])
    recall_value = recall_score(df_test['bugBinary'], df_test['predictBinary'], labels=[0, 1])
    precision_value = precision_score(df_test['bugBinary'], df_test['predictBinary'], labels=[0, 1])
    f1_value = f1_score(df_test['bugBinary'], df_test['predictBinary'], labels=[0, 1])
    gm_value = (recall_value * tnr_value) ** 0.5
    pfr = recall_value
    pdr = fpr  # fp / (fp + tn)
    bpp_value = 1 - (((0 - pfr) ** 2 + (1 - pdr) ** 2) * 0.5) ** 0.5

    valueOfbugBinary = df_test["predictBinary"].value_counts()  # 0 和 1 的各自的个数

    if len(valueOfbugBinary) <= 1:
        if valueOfbugBinary.keys()[0] == 0:
            value_0 = valueOfbugBinary[0]
            value_1 = 0
        else:
            value_0 = 0
            value_1 = valueOfbugBinary[1]
    else:
        value_0 = valueOfbugBinary[0]
        value_1 = valueOfbugBinary[1]

    if auc_value > 1 or auc_value < 0:
        auc_value = 0.5
    elif auc_value < 0.5:
        auc_value = 1 - auc_value
    Q1 = auc_value / (2 - auc_value)
    Q2 = 2 * auc_value * auc_value / (1 + auc_value)
    auc_value_variance = auc_value * (1 - auc_value) + (value_1 - 1) * (Q1 - auc_value * auc_value) \
                         + (value_0 - 1) * (Q2 - auc_value * auc_value)
    auc_value_variance = auc_value_variance / (value_0 * value_1)

    return precision_value, recall_value, auc_value, auc_value_variance, gm_value, f1_value, bpp_value


def thresholds_test(working_dir, result_dir):
    import os
    import csv
    import pandas as pd
    from sklearn.metrics import recall_score, precision_score, f1_score, roc_curve, auc, roc_auc_score, confusion_matrix

    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.width', 5000)

    threshold_dir = working_dir + "staticThresholdsMeta\\"
    testing_dir = working_dir + "pyData\\data_defects_java_testing\\"
    pearson_dir = working_dir + "pearsonMeta\\"

    # read_csv(path, keep_default_na=False, na_values=[""])  只有一个空字段将被识别为NaN
    df_metric_names = pd.read_csv(threshold_dir + "universal_thresholds_meta.csv", keep_default_na=False,
                                  na_values=[""])
    df_pearson = pd.read_csv(pearson_dir + "Pearson_effects_meta.csv", keep_default_na=False, na_values=[""])
    metric_names = sorted(set(df_metric_names.metric.values.tolist()))

    with open(testing_dir + 'List.txt') as l:
        lines = l.readlines()

    with open(threshold_dir + 'List.txt') as l_t:
        lines_t = l_t.readlines()

    # stores predictive performance of the universal threshold of each metric on the testing project
    performnce_columns = ['fileName', 'metric', 'Sample_size', 'threshold', 'threshold_stdError', 'precision', 'recall',
                          'auc', 'auc_var', 'gm', 'f1', 'bpp']
    universal_t_performance = pd.DataFrame(columns=performnce_columns, dtype=object)
    auc_t_performance = pd.DataFrame(columns=performnce_columns, dtype=object)
    bender_t_performance = pd.DataFrame(columns=performnce_columns, dtype=object)
    bpp_t_performance = pd.DataFrame(columns=performnce_columns, dtype=object)
    gm_t_performance = pd.DataFrame(columns=performnce_columns, dtype=object)
    mfm_t_performance = pd.DataFrame(columns=performnce_columns, dtype=object)

    for line in lines:
        project = line.replace("\n", "")

        for root, dirs, files in os.walk(testing_dir + project):
            for name in files:
                print(project, name)

                df_test = pd.read_csv(testing_dir + project + '\\' + name)
                # drop all rows that have any NaN values,删除表中含有任何NaN的行,并重新设置行号
                df_test = df_test.dropna(axis=0, how='any', inplace=False).reset_index(drop=True)

                for metric in metric_names:
                    print(name, metric)

                    df_test = df_test[~df_test[metric].isin(['undef', 'undefined'])]
                    # 由于bug中存储的是缺陷个数,转化为二进制存储,若x>2,则可预测bug为3个以上的阈值,其他类推
                    df_test['bugBinary'] = df_test.bug.apply(lambda x: 1 if x > 0 else 0)
                    # print(df_test.head())
                    # print(df_test.shape)
                    # print(df_test['bugBinary'].sum())
                    # print(df_test['bugBinary'].values)
                    # print(len(df_test[df_test[metric] != 0]))
                    # print(len(df_test[df_test[metric] != 0]) < 6)
                    if len(set(df_test['bugBinary'].values)) == 1:
                        continue

                    if len(df_test[df_test[metric] != 0]) < 6:
                        continue

                    for line_t in lines_t:
                        file_t = line_t.replace("\n", "")
                        # print('the file is ', file_t, file_t.split('_')[0])

                        df_threshold = pd.read_csv(threshold_dir + file_t, keep_default_na=False, na_values=[""])
                        # print(df_threshold.head())
                        # print(df_threshold[df_threshold['metric'] == metric])
                        # print(metric not in df_threshold.metric.values.tolist())
                        # print(df_threshold.metric.values.tolist())

                        if metric not in df_threshold.metric.values.tolist():
                            # append a nan value to universal_performance
                            continue

                        metric_t = df_threshold[df_threshold['metric'] == metric].Pooled_meta_threshold.values[0]
                        metric_t_stdError = \
                            df_threshold[df_threshold['metric'] == metric].Pooled_meta_threshold_stdError.values[0]

                        metric_p = df_pearson[df_pearson['metric'] == metric].Pearson_effects_meta.values[0]

                        precision, recall, auc, auc_var, gm, f1, bpp = \
                            predictive_performance(metric, metric_p, metric_t, df_test)
                        print(name, metric, file_t.split('_')[0], precision, recall, auc, auc_var, gm, f1, bpp)
                        if file_t.split('_')[0] == 'auc':
                            auc_t_performance = auc_t_performance.append(
                                {'fileName': name[:-4], 'metric': metric, 'Sample_size': len(df_test),
                                 'threshold': metric_t, 'threshold_stdError': metric_t_stdError, 'precision': precision,
                                 'recall': recall, 'auc': auc, 'auc_var': auc_var, 'gm': gm, 'f1': f1, 'bpp': bpp},
                                ignore_index=True)
                            auc_t_performance.to_csv(result_dir + 'auc_t_performance.csv', index=False)

                        if file_t.split('_')[0] == 'bender':
                            bender_t_performance = bender_t_performance.append(
                                {'fileName': name[:-4], 'metric': metric, 'Sample_size': len(df_test),
                                 'threshold': metric_t, 'threshold_stdError': metric_t_stdError, 'precision': precision,
                                 'recall': recall, 'auc': auc, 'auc_var': auc_var, 'gm': gm, 'f1': f1, 'bpp': bpp},
                                ignore_index=True)
                            bender_t_performance.to_csv(result_dir + 'bender_t_performance.csv', index=False)

                        if file_t.split('_')[0] == 'bpp':
                            bpp_t_performance = bpp_t_performance.append(
                                {'fileName': name[:-4], 'metric': metric, 'Sample_size': len(df_test),
                                 'threshold': metric_t, 'threshold_stdError': metric_t_stdError, 'precision': precision,
                                 'recall': recall, 'auc': auc, 'auc_var': auc_var, 'gm': gm, 'f1': f1, 'bpp': bpp},
                                ignore_index=True)
                            bpp_t_performance.to_csv(result_dir + 'bpp_t_performance.csv', index=False)

                        if file_t.split('_')[0] == 'gm':
                            gm_t_performance = gm_t_performance.append(
                                {'fileName': name[:-4], 'metric': metric, 'Sample_size': len(df_test),
                                 'threshold': metric_t, 'threshold_stdError': metric_t_stdError, 'precision': precision,
                                 'recall': recall, 'auc': auc, 'auc_var': auc_var, 'gm': gm, 'f1': f1, 'bpp': bpp},
                                ignore_index=True)
                            gm_t_performance.to_csv(result_dir + 'gm_t_performance.csv', index=False)

                        if file_t.split('_')[0] == 'mfm':
                            mfm_t_performance = mfm_t_performance.append(
                                {'fileName': name[:-4], 'metric': metric, 'Sample_size': len(df_test),
                                 'threshold': metric_t, 'threshold_stdError': metric_t_stdError, 'precision': precision,
                                 'recall': recall, 'auc': auc, 'auc_var': auc_var, 'gm': gm, 'f1': f1, 'bpp': bpp},
                                ignore_index=True)
                            mfm_t_performance.to_csv(result_dir + 'mfm_t_performance.csv', index=False)

                        if file_t.split('_')[0] == 'universal':
                            universal_t_performance = universal_t_performance.append(
                                {'fileName': name[:-4], 'metric': metric, 'Sample_size': len(df_test),
                                 'threshold': metric_t, 'threshold_stdError': metric_t_stdError, 'precision': precision,
                                 'recall': recall, 'auc': auc, 'auc_var': auc_var, 'gm': gm, 'f1': f1, 'bpp': bpp},
                                ignore_index=True)
                            universal_t_performance.to_csv(result_dir + 'universal_t_performance.csv', index=False)

                    # break

        # break


if __name__ == '__main__':
    import os
    import sys
    import csv
    import math
    import time
    import random
    import shutil
    from datetime import datetime
    import pandas as pd
    import numpy as np

    s_time = time.time()

    working_Directory = "F:\\MTmeta\\"
    result_Directory = "F:\\MTmeta\\staticThresholdsTesting\\"
    os.chdir(working_Directory)

    thresholds_test(working_Directory, result_Directory)

    e_time = time.time()
    execution_time = e_time - s_time

    print("The __name__ is ", __name__, ".\nFrom ", time.asctime(time.localtime(s_time)), " to ",
          time.asctime(time.localtime(e_time)), ",\nThis", os.path.basename(sys.argv[0]), "ended within",
          execution_time, "(s), or ", (execution_time / 60), " (m).")
