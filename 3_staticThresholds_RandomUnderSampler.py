#!/usr/bin/env python
# encoding:utf-8
"""
Author : Yuanqing Mei
Date : 2021/3/16
Time: 17:59
File: 3_staticThresholds.py
HomePage : http://github.com/yuanqingmei
Email : dg1533019@smail.nju.edu.cn

Use Bender method to derive metric threshold in each system.


Reference:
[1]  Bender, R. Quantitative risk assessment in epidemiological studies investigating threshold effects.
     Biometrical Journal, 41 (1999), 305-319.（计算VARL的SE（标准误）的参考文献P310）
"""


# input: a dataframe including the 'bug', 'bugBinary', 'intercept' and the metric column
# output: 5 five kinds of threshold
# note that the dataframe should input astype(float), i.e., bender_auc_threshold(df.astype(float))
def bender_auc_threshold(df):
    import statsmodels.api as sm
    from imblearn.under_sampling import RandomUnderSampler
    from sklearn.metrics import recall_score, precision_score, f1_score, roc_curve, auc, roc_auc_score, confusion_matrix

    metric_name = ''
    for col in df.columns:
        if col not in ['bug', 'bugBinary', 'intercept']:
            metric_name = col

    try:
        rus = RandomUnderSampler(random_state=0)
        X_resampled, y_resampled = rus.fit_resample(df.loc[:, [metric_name, 'intercept']], df['bugBinary'])
        print(y_resampled)
        print(X_resampled)
        logit = sm.Logit(y_resampled, X_resampled)   # RandomUnderSampler
        # logit = sm.Logit(df['bugBinary'], df.loc[:, [metric_name, 'intercept']])   # 原始数据
        result_logit = logit.fit()
    except Exception as err1:
        print(err1)
        return 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

    tau = result_logit.params[0]
    beta = result_logit.params[1]
    tau_pvalue = result_logit.pvalues[0]
    beta_removed_pvalue = result_logit.pvalues[1]

    # 1. bender method to derive threshold for metric_inverted_comma
    if tau == 0 or tau_pvalue > 0.05:
        bender_t = 0
    else:
        valueOfbugBinary = df['bugBinary'].value_counts()  # 0 和 1 的各自的个数

        # 用缺陷为大于0的模块数占所有模块之比
        BaseProbability_1 = valueOfbugBinary[1] / (valueOfbugBinary[0] + valueOfbugBinary[1])
        # 计算VARL阈值
        bender_t = (np.log(BaseProbability_1 / (1 - BaseProbability_1)) - beta) / tau

    if bender_t > float(df[metric_name].max()) or bender_t < float(df[metric_name].min()):
        bender_t = 0

    # 2. 依次用该度量每一个值作为阈值计算出各自预测性能值,然后选择预测性能值最大的作为阈值,分别定义存入五个值list,取最大值和最大值的下标值
    AUCs = []
    GMs = []
    BPPs = []
    MFMs = []

    auc_max_value = 0
    gm_max_value = 0
    bpp_max_value = 0
    f1_max_value = 0

    i_auc_max = 0
    i_gm_max = 0
    i_bpp_max = 0
    i_f1_max = 0

    # 判断每个度量与bug之间的关系,用于阈值判断正反例
    Corr_metric_bug = df.loc[:, [metric_name, 'bug']].corr('spearman')

    Spearman_value = Corr_metric_bug[metric_name][1]
    Pearson_value = 2 * np.sin(np.pi * Spearman_value / 6)

    # the i value in this loop, is the subscript value in the list of AUCs, GMs etc.
    for i in range(len(df)):

        t = df.loc[i, metric_name]
        # if Corr_metric_bug[metric_name][1] < 0:
        if Pearson_value < 0:
            df['predictBinary'] = df[metric_name].apply(lambda x: 1 if x <= t else 0)
        else:
            df['predictBinary'] = df[metric_name].apply(lambda x: 1 if x >= t else 0)

        # confusion_matrix()函数中需要给出label, 0和1，否则该函数算不出TP,因为不知道哪个标签是poistive.
        c_matrix = confusion_matrix(df["bugBinary"], df['predictBinary'], labels=[0, 1])
        tn, fp, fn, tp = c_matrix.ravel()

        if (tn + fp) == 0:
            tnr_value = 0
        else:
            tnr_value = tn / (tn + fp)

        if (fp + tn) == 0:
            fpr = 0
        else:
            fpr = fp / (fp + tn)

        auc_value = roc_auc_score(df['bugBinary'], df['predictBinary'])
        recall_value = recall_score(df['bugBinary'], df['predictBinary'], labels=[0, 1])
        # precision_value = precision_score(df['bugBinary'], df['predictBinary'], labels=[0, 1])
        f1_value = f1_score(df['bugBinary'], df['predictBinary'], labels=[0, 1])

        gm_value = (recall_value * tnr_value) ** 0.5
        pfr = recall_value
        pdr = fpr  # fp / (fp + tn)
        bpp_value = 1 - (((0 - pfr) ** 2 + (1 - pdr) ** 2) * 0.5) ** 0.5

        AUCs.append(auc_value)
        GMs.append(gm_value)
        BPPs.append(bpp_value)
        MFMs.append(f1_value)

        # 求出上述五个list中最大值，及对应的i值，可能会有几个值相同，且为最大值，则取第一次找到那个值(i)为阈值
        if auc_value > auc_max_value:
            auc_max_value = auc_value
            i_auc_max = i

        if gm_value > gm_max_value:
            gm_max_value = gm_value
            i_gm_max = i

        if bpp_value > bpp_max_value:
            bpp_max_value = bpp_value
            i_bpp_max = i

        if f1_value > f1_max_value:
            f1_max_value = f1_value
            i_f1_max = i

    # 计算auc阈值,包括其他四个类型阈值
    auc_t = df.loc[i_auc_max, metric_name]
    gm_t = df.loc[i_gm_max, metric_name]
    bpp_t = df.loc[i_bpp_max, metric_name]
    mfm_t = df.loc[i_f1_max, metric_name]

    list_t = []
    if bender_t != 0:
        list_t.append(bender_t)
    if auc_t != 0:
        list_t.append(auc_t)
    if gm_t != 0:
        list_t.append(gm_t)
    if bpp_t != 0:
        list_t.append(bpp_t)
    if mfm_t != 0:
        list_t.append(mfm_t)
    avg_t = np.mean(list_t)
    var_t = np.var(list_t)

    return Pearson_value, tau, tau_pvalue, bender_t, auc_t, gm_t, bpp_t, mfm_t, list_t, avg_t, var_t


def static_threshold(working_dir, result_dir):
    # display all columns and rows, and set the item of row of dataframe
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.width', 5000)

    working_directory = working_dir
    result_directory = result_dir
    os.chdir(working_directory)

    with open(working_directory + "List.txt") as l_all:
        lines_all = l_all.readlines()

    for line in lines_all:
        project = line.replace("\n", "")
        print("the file is ", project)

        for root, dirs, files in os.walk(working_directory + project):
            for name in files:

                if not os.path.exists(result_directory + project):
                    os.mkdir(result_directory + project)

                staticThreshold = pd.DataFrame(
                    columns=['version', 'metric', 'pearson', 'tau', 'tau_pvalue', 'bender_t', 'auc_t', 'gm_t',
                             'bpp_t', 'mfm_t', 'avg_t', 'var_t'], dtype=object)

                if os.path.exists(result_directory + project + '\\' + 'staticThreshold_' + name):
                    print("the csv file is created in last execution, so it will not be created this time.")
                    continue

                df_name = pd.read_csv(working_directory + project + '\\' + name)

                # exclude the non metric fields and 31 size metrics
                non_metric = ["relName", "className", 'prevsloc', 'currsloc', 'addedsloc', 'deletedsloc', 'changedsloc',
                              'totalChangedsloc', "bug"]

                # metric_data stores the metric fields (102 items)
                def fun_1(m):
                    return m if m not in non_metric else None

                metric_data = filter(fun_1, df_name.columns)

                for metric in metric_data:
                    print("the current file is ", name, "the current metric is ", metric)

                    # if metric != 'SLOC':
                    #     continue

                    # 由于bug中存储的是缺陷个数,转化为二进制存储,若x>2,则可预测bug为3个以上的阈值,其他类推
                    df_name['bugBinary'] = df_name.bug.apply(lambda x: 1 if x > 0 else 0)

                    # 删除度量中空值和undef值
                    df_metric = df_name[~df_name[metric].isin(['undef', 'undefined'])].loc[:, ['bug', 'bugBinary',
                                                                                               metric]]

                    df_metric = df_metric.dropna(subset=[metric]).reset_index(drop=True)

                    # exclude those data sets in which the metric m has fewer than six non-zero data points
                    # (each corresponding to a class).
                    # if len(df_metric) - len(df_metric[df_metric[metric] == 0]) < 6:
                    if len(df_metric[df_metric[metric] != 0]) < 6:
                        continue

                    df_metric['intercept'] = 1.0

                    pearson_t0, tau_t0, tau_p_t0, bender_t0, auc_t0, gm_t0, bpp_t0, mfm_t0, list_t0, avg_t0, var_t0 = \
                        bender_auc_threshold(df_metric.loc[:, ['bug', 'bugBinary', metric, 'intercept']].astype(float))
                    print(pearson_t0, tau_t0, tau_p_t0, bender_t0, auc_t0, gm_t0, bpp_t0, mfm_t0, list_t0, avg_t0,
                          var_t0)

                    # if pearson correlation is zero, continue
                    if pearson_t0 == 0:
                        continue

                    staticThreshold = staticThreshold.append(
                        {'version': name[:-4], 'metric': metric, 'pearson': pearson_t0, 'tau': tau_t0,
                         'tau_pvalue': tau_p_t0, 'bender_t': bender_t0, 'auc_t': auc_t0, 'gm_t': gm_t0, 'bpp_t': bpp_t0,
                         'mfm_t': mfm_t0, 'avg_t': avg_t0, 'var_t': var_t0}, ignore_index=True)

                    # break
                staticThreshold.to_csv(result_directory + project + '\\' + 'staticThreshold_' + name, index=False)
                # break
        # break


if __name__ == "__main__":
    import os
    import csv
    import sys
    import time
    import numpy as np
    from sklearn.model_selection import StratifiedKFold
    import pandas as pd
    import statsmodels.api as sm
    from sklearn.metrics import roc_curve, auc
    from sklearn.metrics import confusion_matrix

    s_time = time.time()

    working_Directory = "F:\\MTmeta\\pyData\\data_defects_java_training\\"
    result_Directory = "F:\\MTmeta\\staticThresholds\\"
    os.chdir(working_Directory)
    print(os.getcwd())

    static_threshold(working_Directory, result_Directory)

    # bender_threshold()
    e_time = time.time()
    execution_time = e_time - s_time
    print("The __name__ is ", __name__, ".\nFrom ", time.asctime(time.localtime(s_time)), " to ",
          time.asctime(time.localtime(e_time)), ",\nThis", os.path.basename(sys.argv[0]), "ended within",
          execution_time, "(s), or ", (execution_time / 60), " (m).")
