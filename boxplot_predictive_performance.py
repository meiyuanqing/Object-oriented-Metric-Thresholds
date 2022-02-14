#!/usr/bin/env python
# encoding:utf-8
"""
Author : Yuanqing Mei
Date : 2022/1/9
Time: 19:39
File: boxplot_predictive_performance.py
HomePage : https://github.com/meiyuanqing
Email : dg1533019@smail.nju.edu.cn

Compute the avg std max min values and draw the box plot of precision and recall.

"""


def predictvie_plot(working_dir, plot_dir):

    import os
    import csv
    import pandas as pd
    import matplotlib.pyplot as plt

    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.width', 5000)

    plt.rcParams['savefig.dpi'] = 900  # 图片像素
    plt.rcParams['figure.figsize'] = (8.0, 4.0)

    os.chdir(working_dir)

    df = pd.read_csv(working_dir + "universal_t_performance.csv", keep_default_na=False, na_values=[""])

    # df = df.drop(axis=0, how='any', inplace=False).reset_index(drop=True)

    metric_names = sorted(set(df.metric.values.tolist()))
    print("the metric_names are ", df.columns.values.tolist())
    print("the metric_names are ", metric_names)
    print("the len metric_names are ", len(metric_names))

    with open(working_dir + "precision_table.csv", 'a+', encoding="utf-8", newline='') as table:

        writer = csv.writer(table)
        if os.path.getsize(working_dir + "precision_table.csv") == 0:
            writer.writerow(["metric", "Sample_size", "recall_max", "recall_min", "recall_median", "recall_mean",
                             "recall_variance", "precision_max", "precision_min", "precision_median", "precision_mean",
                             "precision_variance", "f1_max", "f1_min", "f1_median", "f1_mean", "f1_variance"])

        # 需要把同类型所有的度量的性能指标画在一张图上，定义一个空数据框，用于存入行数相同的度量性能结果
        size = ["NA", "NAIMP", "NM", "NMIMP", "NumPara", "SLOC", "stms"]

        cohesion = ['LCOM1', 'LCOM2', 'LCOM3', 'LCOM4', 'Co', 'NewCo', 'LCOM5', 'NewLCOM5', 'TCC', 'LCC', 'ICH', 'OCC',
                    'PCC', 'DCd', 'DCi', 'CAMC', 'NHD', 'SNHD']

        coupling = ['ACAIC', 'ACMIC', 'AMMIC', 'DMMEC', 'OCAEC', 'OCAIC', 'OCMEC', 'OCMIC', 'OMMEC', 'OMMIC', 'DCAEC',
                    'DCMEC', 'CBI', 'CBO', 'DAC', 'ICP', 'IHICP', 'MPC', 'NIHICP', 'RFC']

        inheritance = ["AID", "CLD", "DIT", "DP", "DPA", "DPD", "NMA", "NMI", "NMO", "NOA", "NOC", "NOD", "NOP", "SIX",
                       "SP", "SPA", "SPD"]

        size_recall_df = pd.DataFrame(columns=size)
        size_precision_df = pd.DataFrame(columns=size)
        size_auc_df = pd.DataFrame(columns=size)
        cohesion_recall_df = pd.DataFrame(columns=cohesion)
        cohesion_precision_df = pd.DataFrame(columns=cohesion)
        cohesion_auc_df = pd.DataFrame(columns=cohesion)
        coupling_recall_df = pd.DataFrame(columns=coupling)
        coupling_precision_df = pd.DataFrame(columns=coupling)
        coupling_auc_df = pd.DataFrame(columns=coupling)
        inheritance_recall_df = pd.DataFrame(columns=inheritance)
        inheritance_precision_df = pd.DataFrame(columns=inheritance)
        inheritance_auc_df = pd.DataFrame(columns=inheritance)

        for metric in metric_names:

            print("the current metric is ", metric)
            metric_row = [metric, len(df[df["metric"] == metric].loc[:, "recall"])]

            metric_row.append(df[df["metric"] == metric].loc[:, "recall"].max())
            metric_row.append(df[df["metric"] == metric].loc[:, "recall"].min())
            metric_row.append(df[df["metric"] == metric].loc[:, "recall"].median())
            metric_row.append(df[df["metric"] == metric].loc[:, "recall"].mean())
            metric_row.append(df[df["metric"] == metric].loc[:, "recall"].var())

            metric_row.append(df[df["metric"] == metric].loc[:, "precision"].max())
            metric_row.append(df[df["metric"] == metric].loc[:, "precision"].min())
            metric_row.append(df[df["metric"] == metric].loc[:, "precision"].median())
            metric_row.append(df[df["metric"] == metric].loc[:, "precision"].mean())
            metric_row.append(df[df["metric"] == metric].loc[:, "precision"].var())

            metric_row.append(df[df["metric"] == metric].loc[:, "f1"].max())
            metric_row.append(df[df["metric"] == metric].loc[:, "f1"].min())
            metric_row.append(df[df["metric"] == metric].loc[:, "f1"].median())
            metric_row.append(df[df["metric"] == metric].loc[:, "f1"].mean())
            metric_row.append(df[df["metric"] == metric].loc[:, "f1"].var())

            print("the mean value of recall is ", df[df["metric"] == metric].loc[:, "recall"].mean())
            writer.writerow(metric_row)

            if metric in size:
                size_recall_df[metric] = df[df["metric"] == metric].loc[:, "recall"].dropna(axis=0,
                                  how='any', inplace=False).reset_index(drop=True)
                size_precision_df[metric] = df[df["metric"] == metric].loc[:, "precision"].dropna(axis=0,
                                  how='any', inplace=False).reset_index(drop=True)
                size_auc_df[metric] = df[df["metric"] == metric].loc[:, "auc"].dropna(axis=0,
                                  how='any', inplace=False).reset_index(drop=True)
            if metric in cohesion:
                cohesion_recall_df[metric] = df[df["metric"] == metric].loc[:, "recall"].dropna(axis=0,
                                  how='any', inplace=False).reset_index(drop=True)
                cohesion_precision_df[metric] = df[df["metric"] == metric].loc[:, "precision"].dropna(axis=0,
                                  how='any', inplace=False).reset_index(drop=True)
                cohesion_auc_df[metric] = df[df["metric"] == metric].loc[:, "auc"].dropna(axis=0,
                                  how='any', inplace=False).reset_index(drop=True)
            if metric in coupling:
                coupling_recall_df[metric] = df[df["metric"] == metric].loc[:, "recall"].dropna(axis=0,
                                  how='any', inplace=False).reset_index(drop=True)
                coupling_precision_df[metric] = df[df["metric"] == metric].loc[:, "precision"].dropna(axis=0,
                                  how='any', inplace=False).reset_index(drop=True)
                coupling_auc_df[metric] = df[df["metric"] == metric].loc[:, "auc"].dropna(axis=0,
                                  how='any', inplace=False).reset_index(drop=True)
            if metric in inheritance:
                inheritance_recall_df[metric] = df[df["metric"] == metric].loc[:, "recall"].dropna(axis=0,
                                  how='any', inplace=False).reset_index(drop=True)
                inheritance_precision_df[metric] = df[df["metric"] == metric].loc[:, "precision"].dropna(axis=0,
                                  how='any', inplace=False).reset_index(drop=True)
                inheritance_auc_df[metric] = df[df["metric"] == metric].loc[:, "auc"].dropna(axis=0,
                                  how='any', inplace=False).reset_index(drop=True)

        print(size_recall_df)
        plt.rcParams['savefig.dpi'] = 900  # 图片像素
        plt.rcParams['figure.figsize'] = (20.0, 6.0)
        size_recall_df.plot.box(title="(d) Box Plot of Recall values for Size Metrics")
        plt.grid(linestyle="--", alpha=0.3)
        # plt.xticks(rotation=0, fontsize=9.0)
        plt.xticks(rotation=0)
        plt.gcf().subplots_adjust(bottom=0.34) #这是增加下面标签的显示宽度，20210807找了下午5个小时
        plt.savefig(plot_dir + 'SizeRecallMetrics.png')
        plt.close()

        size_precision_df.plot.box(title="(d) Box Plot of Precision values for Size Metrics")
        plt.grid(linestyle="--", alpha=0.3)
        plt.xticks(rotation=0)
        plt.gcf().subplots_adjust(bottom=0.34)
        plt.savefig(plot_dir + 'SizePrecisionMetrics.png')
        plt.close()

        size_auc_df.plot.box(title="(d) Box Plot of AUC values for Size Metrics")
        plt.grid(linestyle="--", alpha=0.3)
        plt.xticks(rotation=0)
        plt.gcf().subplots_adjust(bottom=0.34)
        plt.savefig(plot_dir + 'SizeAUCMetrics.png')
        plt.close()

        cohesion_recall_df.plot.box(title="(a) Box Plot of Recall values for Cohesion Metrics")
        plt.grid(linestyle="--", alpha=0.3)
        # plt.xticks(rotation=0, fontsize=9.0)
        plt.xticks(rotation=0)
        plt.gcf().subplots_adjust(bottom=0.34) #这是增加下面标签的显示宽度，20210807找了下午5个小时
        plt.savefig(plot_dir + 'cohesionRecallMetrics.png')
        plt.close()

        cohesion_precision_df.plot.box(title="(a) Box Plot of Precision values for Cohesion Metrics")
        plt.grid(linestyle="--", alpha=0.3)
        plt.xticks(rotation=0)
        plt.gcf().subplots_adjust(bottom=0.34)
        plt.savefig(plot_dir + 'cohesionPrecisionMetrics.png')
        plt.close()

        cohesion_auc_df.plot.box(title="(a) Box Plot of AUC values for Cohesion Metrics")
        plt.grid(linestyle="--", alpha=0.3)
        plt.xticks(rotation=0)
        plt.gcf().subplots_adjust(bottom=0.34)
        plt.savefig(plot_dir + 'cohesionAUCMetrics.png')
        plt.close()

        coupling_recall_df.plot.box(title="(b) Box Plot of Recall values for Coupling Metrics")
        plt.grid(linestyle="--", alpha=0.3)
        # plt.xticks(rotation=0, fontsize=9.0)
        plt.xticks(rotation=0)
        plt.gcf().subplots_adjust(bottom=0.34) #这是增加下面标签的显示宽度，20210807找了下午5个小时
        plt.savefig(plot_dir + 'couplingRecallMetrics.png')
        plt.close()

        coupling_precision_df.plot.box(title="(b) Box Plot of Precision values for Coupling Metrics")
        plt.grid(linestyle="--", alpha=0.3)
        plt.xticks(rotation=0)
        plt.gcf().subplots_adjust(bottom=0.34)
        plt.savefig(plot_dir + 'couplingPrecisionMetrics.png')
        plt.close()

        coupling_auc_df.plot.box(title="(b) Box Plot of AUC values for Coupling Metrics")
        plt.grid(linestyle="--", alpha=0.3)
        plt.xticks(rotation=0)
        plt.gcf().subplots_adjust(bottom=0.34)
        plt.savefig(plot_dir + 'couplingAUCMetrics.png')
        plt.close()

        inheritance_recall_df.plot.box(title="(c) Box Plot of Recall values for Inheritance Metrics")
        plt.grid(linestyle="--", alpha=0.3)
        # plt.xticks(rotation=0, fontsize=9.0)
        plt.xticks(rotation=0)
        plt.gcf().subplots_adjust(bottom=0.34) #这是增加下面标签的显示宽度，20210807找了下午5个小时
        plt.savefig(plot_dir + 'inheritanceRecallMetrics.png')
        plt.close()

        inheritance_precision_df.plot.box(title="(c) Box Plot of Precision values for Inheritance Metrics")
        plt.grid(linestyle="--", alpha=0.3)
        plt.xticks(rotation=0)
        plt.gcf().subplots_adjust(bottom=0.34)
        plt.savefig(plot_dir + 'inheritancePrecisionMetrics.png')
        plt.close()

        inheritance_auc_df.plot.box(title="(c) Box Plot of AUC values for Inheritance Metrics")
        plt.grid(linestyle="--", alpha=0.3)
        plt.xticks(rotation=0)
        plt.gcf().subplots_adjust(bottom=0.34)
        plt.savefig(plot_dir + 'inheritanceAUCMetrics.png')
        plt.close()


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

    working_Directory = "F:\\MTmeta\\staticThresholdsTesting\\"
    plot_Directory = "F:\\MTmeta\\staticThresholdsTesting\\boxPlot\\"
    os.chdir(working_Directory)

    predictvie_plot(working_Directory, plot_Directory)

    e_time = time.time()
    execution_time = e_time - s_time

    print("The __name__ is ", __name__, ".\nFrom ", time.asctime(time.localtime(s_time)), " to ",
          time.asctime(time.localtime(e_time)), ",\nThis", os.path.basename(sys.argv[0]), "ended within",
          execution_time, "(s), or ", (execution_time / 60), " (m).")
