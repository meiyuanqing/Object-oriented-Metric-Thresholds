#!/usr/bin/env python
# encoding:utf-8
"""
Author : Yuanqing Mei
Date : 2022/1/6
Time: 23:06
File: 4_staticThreshold_meta.py
HomePage : https://github.com/meiyuanqing
Email : dg1533019@smail.nju.edu.cn

Deriving Pooled (methods) thresholds by meta-analysis (including sensitivity analysis).

Five supervised methods: Bender, ROC, BPP, MFM, GM.

For each OO metric on 66 projects, five supervised methods derived five thresholds.

So, there are 66*5=330 samples to perform meta-analysis.

"""

import os
import csv
from scipy.stats import norm  # norm.cdf() the cumulative normal distribution function in Python
from scipy import stats  # 根据卡方分布计算p值: p_value=1.0-stats.chi2.cdf(chisquare,freedom_degree)
import numpy as np
import pandas as pd


# 输入：两个匿名数组，effect_size中存放每个study的effect size，variance存放对应的方差
# 输出：fixed effect model固定效应元分析后的结果，包括
#      (1)fixedMean：固定效应元分析后得到的效应平均值；(2) fixedStdError：固定效应元分析的效应平均值对应的标准错
def fixed_effect_meta_analysis(effect_size, variance):
    fixed_weight = []
    sum_Wi = 0
    sum_WiYi = 0
    d = {}  # return a dict
    study_number = len(variance)
    for i in range(study_number):
        if variance[i] == 0:
            continue
        fixed_weight.append(1 / variance[i])
        sum_Wi = sum_Wi + fixed_weight[i]
        sum_WiYi = sum_WiYi + effect_size[i] * fixed_weight[i]
    fixedMean = sum_WiYi / sum_Wi  # 固定模型元分析后得到的效应平均值
    fixedStdError = (1 / sum_Wi) ** 0.5  # 固定模型元分析的效应平均值对应的标准错
    d['fixedMean'] = fixedMean
    d['fixedStdError'] = fixedStdError
    return d


# 输入：两个匿名数组，effect_size中存放每个study的effect size，variance存放对应的方差
# 输出：random effect model随机效应元分析后的结果，包括：
#      (1) randomMean：随机模型元分析后得到的效应平均值； (2) randomStdError：随机模型元分析的效应平均值对应的标准错
def random_effect_meta_analysis(effect_size, variance):
    sum_Wi = 0
    sum_WiWi = 0
    sum_WiYi = 0  # Sum(Wi*Yi), where i ranges from 1 to k, and k is the number of studies
    sum_WiYiYi = 0  # Sum(Wi*Yi*Yi), where i ranges from 1 to k, and k is the number of studies

    sum_Wistar = 0
    sum_WistarYi = 0
    d = {}  # return a dict

    study_number = len(variance)
    fixed_weight = [0 for i in range(study_number)]  # 固定模型对应的权值
    random_weight = [0 for i in range(study_number)]  # 随机模型对应的权值

    for i in range(study_number):
        if variance[i] == 0:
            continue
        fixed_weight[i] = 1 / variance[i]
        sum_Wi = sum_Wi + fixed_weight[i]
        sum_WiWi = sum_WiWi + fixed_weight[i] * fixed_weight[i]
        sum_WiYi = sum_WiYi + effect_size[i] * fixed_weight[i]
        sum_WiYiYi = sum_WiYiYi + fixed_weight[i] * effect_size[i] * effect_size[i]

    Q = sum_WiYiYi - sum_WiYi * sum_WiYi / sum_Wi
    df = study_number - 1
    C = sum_Wi - sum_WiWi / sum_Wi

    # for PII metric C = 0 20210423
    # 当元分析过程中只有一个study研究时，没有研究间效应，故研究间的方差为零
    # if study_number == 1 or C == 0:
    if study_number == 1:
        T2 = 0
    else:
        T2 = (Q - df) / C  # sample estimate of tau square

    if T2 < 0:
        T2 = 0  # 20210411，若T2小于0，取0,   M.Borenstein[2009] P114

    for i in range(study_number):
        random_weight[i] = 1 / (variance[i] + T2)  # random_weight 随机模型对应的权值

    for i in range(study_number):
        sum_Wistar = sum_Wistar + random_weight[i]
        sum_WistarYi = sum_WistarYi + random_weight[i] * effect_size[i]

    randomMean = sum_WistarYi / sum_Wistar  # 随机模型元分析后得到的效应平均值
    randomStdError = (1 / sum_Wistar) ** 0.5  # 随机模型元分析的效应平均值对应的标准错
    # 当元分析过程中只有一个study研究时，没有研究间异质性，故异质性为零
    if study_number == 1:
        I2 = 0
    else:
        I2 = ((Q - df) / Q) * 100  # Higgins et al. (2003) proposed using a statistic, I2,
        # the proportion of the observed variance reflects real differences in effect size

    if I2 < 0:
        I2 = 0  # 20210420，若I2小于0，取0,   M.Borenstein[2009] P110

    pValue_Q = 1.0 - stats.chi2.cdf(Q, df)  # pValue_Q = 1.0 - stats.chi2.cdf(chisquare, freedom_degree)

    d["C"] = C
    d["mean"] = randomMean
    d["stdError"] = randomStdError
    d["LL_CI"] = randomMean - 1.96 * randomStdError  # The 95% lower limits for the summary effect
    d["UL_CI"] = randomMean + 1.96 * randomStdError  # The 95% upper limits for the summary effect
    # 20210719 adds the 84% CI for the summary effect
    d["LL_CI_84"] = randomMean - 1.4051 * randomStdError  # The 84% lower limits for the summary effect
    d["UL_CI_84"] = randomMean + 1.4051 * randomStdError  # The 84% upper limits for the summary effect

    d["ZValue"] = randomMean / randomStdError  # a Z-value to test the null hypothesis that the mean effect is zero
    # 20210414 双侧检验时需要增加绝对值符号np.abs
    d["pValue_Z"] = 2 * (1 - norm.cdf(np.abs(randomMean / randomStdError)))  # norm.cdf() 返回标准正态累积分布函数值
    d["Q"] = Q
    d["df"] = df
    d["pValue_Q"] = pValue_Q
    d["I2"] = I2
    d["tau"] = T2 ** 0.5
    d["LL_ndPred"] = randomMean - 1.96 * (T2 ** 0.5)  # tau、randomMean 已知情况下的新出现的study的effctsize所落的区间
    d["UL_ndPred"] = randomMean + 1.96 * (T2 ** 0.5)  # tau、randomMean 已知情况下的新出现的study的effctsize所落的区间
    # tau、randomMean 未知情况（估计）下的新出现的study的effctsize所落的区间
    # stats.t.ppf(0.975,df)返回学生t分布单尾alpha=0.025区间点(双尾是alpha=0.05)的函数，它是stats.t.cdf()累积分布函数的逆函数
    d["LL_tdPred"] = randomMean - stats.t.ppf(0.975, df) * ((T2 + randomStdError * randomStdError) ** 0.5)
    # tau、randomMean 未知情况（估计）下的新出现的study的effctsize所落的区间
    d["UL_tdPred"] = randomMean + stats.t.ppf(0.975, df) * ((T2 + randomStdError * randomStdError) ** 0.5)

    fixedMean = sum_WiYi / sum_Wi  # 固定模型元分析后得到的效应平均值
    fixedStdError = (1 / sum_Wi) ** 0.5  # 固定模型元分析的效应平均值对应的标准错
    d['fixedMean'] = fixedMean
    d['fixedStdError'] = fixedStdError
    return d


def getEstimatedK0(effectSizeArray, mean):
    centeredEffectSizeArray = []
    absoluteCenteredEffectSizeArray = []
    size = len(effectSizeArray)
    for i in range(size):
        centeredEffectSizeArray.append(effectSizeArray[i] - mean)
        absoluteCenteredEffectSizeArray.append(np.abs(effectSizeArray[i] - mean))
    sortedArray = sorted(absoluteCenteredEffectSizeArray)
    rank = {sortedArray[0]: 1}  # return a dict
    initialRankValue = 1
    predValue = sortedArray[0]
    for i in range(size):
        if sortedArray[i] > predValue:
            predValue = sortedArray[i]
            initialRankValue += 1
        rank[sortedArray[i]] = initialRankValue
    finalRank = []
    for i in range(size):
        if centeredEffectSizeArray[i] < 0:
            finalRank.append((-1) * rank[absoluteCenteredEffectSizeArray[i]])
        else:
            finalRank.append(rank[absoluteCenteredEffectSizeArray[i]])
    gamma = finalRank[size - 1] + finalRank[0]
    SumPositiveRank = 0
    for i in range(size):
        if finalRank[i] < 0:
            continue
        SumPositiveRank = SumPositiveRank + finalRank[i]
    R0 = int(gamma + 0.5) - 1
    temp = (4 * SumPositiveRank - size * (size + 1)) / (2 * size - 1)
    L0 = int(temp + 0.5)
    if R0 < 0:
        R0 = 0
    if L0 < 0:
        L0 = 0
    return R0, L0


# Duval and Tweedie's trim and fill method
def trimAndFill(effect_size, variance, isAUC):
    effectSizeArray = effect_size
    varianceArray = variance
    size = len(effect_size)
    # 检查是否需要切换方向，因为trim and fill方法假设miss most negative的研究
    flipFunnel = 0
    metaAnalysisForFlip = fixed_effect_meta_analysis(effectSizeArray, varianceArray)
    meanForFlip = metaAnalysisForFlip["fixedMean"]

    tempSorted = sorted(effectSizeArray)
    min = tempSorted[0] - meanForFlip
    max = tempSorted[-1] - meanForFlip

    if np.abs(min) > np.abs(max):
        flipFunnel = 1
        for i in range(size):
            effectSizeArray[i] = (-1) * effectSizeArray[i]

    # 按effect size排序
    merge = []
    for i in range(size):
        merge.append([effect_size[i], variance[i]])
    sortedMerge = sorted(merge)
    OrignalEffectSizeArray = []
    OrignalVarianceArray = []
    for i in range(len(sortedMerge)):
        OrignalEffectSizeArray.append(sortedMerge[i][0])
        OrignalVarianceArray.append(sortedMerge[i][1])
    # 迭代算法，估算k0
    metaAnalysisResult = fixed_effect_meta_analysis(OrignalEffectSizeArray, OrignalVarianceArray)
    mean = metaAnalysisResult["fixedMean"]
    RL = getEstimatedK0(OrignalEffectSizeArray, mean)
    R0 = RL[0]
    L0 = RL[1]
    k0 = L0  # 默认的情况利用L0来估算k0
    if (k0 == 0) or (k0 > size):
        result = random_effect_meta_analysis(effect_size, variance)
        result["k0"] = k0
        return result
    trimmedMean = mean
    change = 1
    count = 0
    while change and (size - k0) > 2 and (count < 1000):
        count += 1
        upperBound = size - k0 - 1
        trimmedEffectSizeArray = []
        trimmedVarianceArray = []
        for i in range(upperBound):
            trimmedEffectSizeArray.append(OrignalEffectSizeArray[i])
            trimmedVarianceArray.append(OrignalVarianceArray[i])
        trimmedMetaAnalysisResult = fixed_effect_meta_analysis(trimmedEffectSizeArray, trimmedVarianceArray)
        trimmedMean = trimmedMetaAnalysisResult["fixedMean"]
        trimmedR0_L0 = getEstimatedK0(OrignalEffectSizeArray, trimmedMean)
        trimmedR0 = trimmedR0_L0[0]
        trimmedL0 = trimmedR0_L0[1]
        k1 = trimmedL0
        if k1 == k0:
            change = 0
        k0 = k1
    filledEffectSizeArray = []
    filledVarianceArray = []

    for j in range(k0):
        imputedEffectSize = 2 * trimmedMean - OrignalEffectSizeArray[size - j - 1]
        imputedVariance = OrignalVarianceArray[size - j - 1]
        filledEffectSizeArray.append(imputedEffectSize)
        filledVarianceArray.append(imputedVariance)
    fullEffectSizeArray = filledEffectSizeArray
    fullVarianceArray = filledVarianceArray
    fullEffectSizeArray.extend(OrignalEffectSizeArray)
    fullVarianceArray.extend(OrignalVarianceArray)
    if flipFunnel:
        newSize = len(fullEffectSizeArray)
        for i in range(newSize):
            fullEffectSizeArray[i] = -1 * fullEffectSizeArray[i]

    if isAUC:
        # AUC应该在0到1之间，否则有错
        filteredFullEffectSizeArray = []
        filteredFullVarianceArray = []
        for i in range(len(fullEffectSizeArray)):
            if fullEffectSizeArray[i] < 0:
                continue
            if fullEffectSizeArray[i] > 1:
                continue
            filteredFullEffectSizeArray.append(fullEffectSizeArray[i])
            filteredFullVarianceArray.append(fullVarianceArray[i])
        result = random_effect_meta_analysis(filteredFullEffectSizeArray, filteredFullVarianceArray)
        finalk0 = len(filteredFullEffectSizeArray) - len(OrignalEffectSizeArray)
    else:
        result = random_effect_meta_analysis(fullEffectSizeArray, fullVarianceArray)
        finalk0 = len(fullEffectSizeArray) - len(OrignalEffectSizeArray)
    result["k0"] = finalk0
    result["flipFunnel"] = flipFunnel
    return result


def do_random_meta(meta_dir, file_name, metric_name, effect_size, variance):
    try:
        resultMetaAnalysis = random_effect_meta_analysis(np.array(effect_size), np.array(variance))

        adjusted_result = trimAndFill(np.array(effect_size), np.array(variance), 0)

        with open(meta_dir + file_name, 'a+', encoding="utf-8", newline='') as f:
            writer_f = csv.writer(f)
            if os.path.getsize(meta_dir + file_name) == 0:
                writer_f.writerow(
                    ["metric", "Pooled_meta_threshold", "Pooled_meta_threshold_stdError", "LL_CI", "UL_CI",
                     "LL_CI_84", "UL_CI_84",
                     "ZValue", "pValue_Z", "Q", "df", "pValue_Q", "I2", "tau", "LL_ndPred", "UL_ndPred",
                     "number_of_effect_size",
                     "k_0", "Pooled_meta_threshold_adjusted", "Pooled_meta_threshold_stdError_adjusted",
                     "LL_CI_adjusted", "UL_CI_adjusted", "LL_CI_84_adjusted", "UL_CI_84_adjusted",
                     "pValue_Z_adjusted", "Q_adjusted", "df_adjusted",
                     "pValue_Q_adjusted", "I2_adjusted", "tau_adjusted", "LL_ndPred_adjusted",
                     "UL_ndPred_adjusted"])
            writer_f.writerow([metric_name, resultMetaAnalysis["mean"], resultMetaAnalysis["stdError"],
                               resultMetaAnalysis["LL_CI"], resultMetaAnalysis["UL_CI"],
                               resultMetaAnalysis["LL_CI_84"], resultMetaAnalysis["UL_CI_84"],
                               resultMetaAnalysis["ZValue"], resultMetaAnalysis["pValue_Z"],
                               resultMetaAnalysis["Q"], resultMetaAnalysis["df"], resultMetaAnalysis["pValue_Q"],
                               resultMetaAnalysis["I2"], resultMetaAnalysis["tau"], resultMetaAnalysis["LL_ndPred"],
                               resultMetaAnalysis["UL_ndPred"], len(effect_size),
                               adjusted_result["k0"], adjusted_result["mean"], adjusted_result["stdError"],
                               adjusted_result["LL_CI"], adjusted_result["UL_CI"],
                               adjusted_result["LL_CI_84"], adjusted_result["UL_CI_84"],
                               adjusted_result["pValue_Z"],
                               adjusted_result["Q"], adjusted_result["df"], adjusted_result["pValue_Q"],
                               adjusted_result["I2"], adjusted_result["tau"], adjusted_result["LL_ndPred"],
                               adjusted_result["UL_ndPred"]])

    except Exception as err1:
        print(err1)


def pool_meta(working_dir, result_dir):

    # display all columns and rows, and set the item of row of dataframe
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.width', 5000)

    metric_dir = working_dir
    meta_dir = result_dir
    os.chdir(metric_dir)
    print(os.getcwd())

    cohesion = ['LCOM1', 'LCOM2', 'LCOM3', 'LCOM4', 'Co', 'NewCo', 'LCOM5', 'NewLCOM5', 'TCC', 'LCC', 'ICH', 'OCC',
                'PCC', 'DCd', 'DCi', 'CAMC', 'NHD', 'SNHD']

    coupling = ['ACAIC', 'ACMIC', 'AMMIC', 'DMMEC', 'OCAEC', 'OCAIC', 'OCMEC', 'OCMIC', 'OMMEC', 'OMMIC', 'DCAEC',
                'DCMEC', 'CBI', 'CBO', 'DAC', 'ICP', 'IHICP', 'MPC', 'NIHICP', 'RFC']

    inheritance = ["AID", "CLD", "DIT", "DP", "DPA", "DPD", "NMA", "NMI", "NMO", "NOA", "NOC", "NOD", "NOP", "SIX",
                   "SP", "SPA", "SPD"]

    size = ["NA", "NAIMP", "NM", "NMIMP", "NumPara", "SLOC", "stms"]

    metrics = cohesion + coupling + inheritance + size

    # stores the threshold of each metric on the training project deriving from 9 methods
    universal_thresholds = pd.DataFrame()

    with open(metric_dir + "List.txt") as l_all:
        lines_all = l_all.readlines()

    # for one metric
    for metric in metrics:

        print("the current metric is ", metric)

        # appends nine method's thresholds and their variances of each metric in training date in turn
        threshold_effects, threshold_variances = [], []
        bender_t_effects, bender_t_variances = [], []
        auc_t_effects, auc_t_variances = [], []
        gm_t_effects, gm_t_variances = [], []
        bpp_t_effects, bpp_t_variances = [], []
        mfm_t_effects, mfm_t_variances = [], []


        for line in lines_all:
            project = line.replace("\n", "")
            print("the file is ", project)

            bender_t, auc_t, gm_t, bpp_t, mfm_t, = [], [], [], [], []

            for root, dirs, files in os.walk(metric_dir + project):
                for name in files:
                    print(metric, name)

                    df_name = pd.read_csv(metric_dir + project + '\\' + name, index_col=False, keep_default_na=False,
                                          na_values=[""])

                    if metric not in df_name.metric.values:
                        continue
                    if df_name[df_name['metric'] == metric]['bender_t'].values[0] != 0:
                        bender_t.append(df_name[df_name['metric'] == metric]['bender_t'].values[0])
                    auc_t.append(df_name[df_name['metric'] == metric]['auc_t'].values[0])
                    gm_t.append(df_name[df_name['metric'] == metric]['gm_t'].values[0])
                    bpp_t.append(df_name[df_name['metric'] == metric]['bpp_t'].values[0])
                    mfm_t.append(df_name[df_name['metric'] == metric]['mfm_t'].values[0])
                    # break

            if len(bender_t) > 0 and np.mean(bender_t) != 0 and np.var(bender_t) != 0:
                bender_t_effects.append(np.mean(bender_t))
                bender_t_variances.append(np.var(bender_t))
                threshold_effects.append(np.mean(bender_t))
                threshold_variances.append(np.var(bender_t))
            if len(auc_t) > 0 and np.mean(auc_t) != 0 and np.var(auc_t) != 0:
                auc_t_effects.append(np.mean(auc_t))
                auc_t_variances.append(np.var(auc_t))
                threshold_effects.append(np.mean(auc_t))
                threshold_variances.append(np.var(auc_t))
            if len(gm_t) > 0 and np.mean(gm_t) != 0 and np.var(gm_t) != 0:
                gm_t_effects.append(np.mean(gm_t))
                gm_t_variances.append(np.var(gm_t))
                threshold_effects.append(np.mean(gm_t))
                threshold_variances.append(np.var(gm_t))
            if len(bpp_t) > 0 and np.mean(bpp_t) != 0 and np.var(bpp_t) != 0:
                bpp_t_effects.append(np.mean(bpp_t))
                bpp_t_variances.append(np.var(bpp_t))
                threshold_effects.append(np.mean(bpp_t))
                threshold_variances.append(np.var(bpp_t))
            if len(mfm_t) > 0 and np.mean(mfm_t) != 0 and np.var(mfm_t) != 0:
                mfm_t_effects.append(np.mean(mfm_t))
                mfm_t_variances.append(np.var(mfm_t))
                threshold_effects.append(np.mean(mfm_t))
                threshold_variances.append(np.var(mfm_t))
            # break

        # do_random_meta(meta_dir, file_name, metric_name, effect_size, variance)
        if len(bender_t_effects) > 0:
            do_random_meta(meta_dir, 'bender_thresholds_meta.csv', metric, bender_t_effects, bender_t_variances)
        if len(auc_t_effects) > 0:
            do_random_meta(meta_dir, 'auc_thresholds_meta.csv', metric, auc_t_effects, auc_t_variances)
        if len(gm_t_effects) > 0:
            do_random_meta(meta_dir, 'gm_thresholds_meta.csv', metric, gm_t_effects, gm_t_variances)
        if len(bpp_t_effects) > 0:
            do_random_meta(meta_dir, 'bpp_thresholds_meta.csv', metric, bpp_t_effects, bpp_t_variances)
        if len(mfm_t_effects) > 0:
            do_random_meta(meta_dir, 'mfm_thresholds_meta.csv', metric, mfm_t_effects, mfm_t_variances)
        if len(threshold_effects) > 0:
            do_random_meta(meta_dir, 'universal_thresholds_meta.csv', metric, threshold_effects, threshold_variances)

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

    working_Directory = "F:\\MTmeta\\staticThresholds\\"
    result_Directory = "F:\\MTmeta\\staticThresholdsMeta\\"

    os.chdir(working_Directory)

    pool_meta(working_Directory, result_Directory)

    e_time = time.time()
    execution_time = e_time - s_time

    print("The __name__ is ", __name__, ".\nFrom ", time.asctime(time.localtime(s_time)), " to ",
          time.asctime(time.localtime(e_time)), ",\nThis", os.path.basename(sys.argv[0]), "ended within",
          execution_time, "(s), or ", (execution_time / 60), " (m).")
