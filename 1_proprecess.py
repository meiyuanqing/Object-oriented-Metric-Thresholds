#!/usr/bin/env python
# encoding:utf-8
"""
Author : Yuanqing Mei
Date : 2021/12/13
Time: 11:44
File: 1_proprecess.py
HomePage : https://github.com/meiyuanqing
Email : dg1533019@smail.nju.edu.cn

OO度量阈值元分析实验数据预处理：
（1）根据JIRA中版本时间2020年12月31日（versionDate.csv）筛选出符合要求的版本。虽然JIRA缺陷数据收集截止时间为2021年9月30日，
    但由于在9月份发布的版本还没有充足的时间发现已发布版本上的缺陷，【正常的窗口期是6个月】故设定该截止日期。
（2）在剩下的版本中，若该版本上没有缺陷的被剔除。因为需要在每个版本上计算阈值，若二分类变量（因变量）只有一个类型变量模型不能学习
    并区分两种状态下的OO度量值差异，故舍去。
（3）在符合（1）（2）的版本中，对每个OO度量，若非零数值小于6个，则被舍去，原因也是不能区分不同的模块。
（4）由于个别内聚度量收集过程中出现“undef”值，也舍去。
（5）在每个项目中8种语言中，只剩下JAVA为扩展名的类。
（6）目前有124个understand度量和PERL度量，需要检查在每个项目是否都，即需要统计出每个OO度量满足上述5个预处理要求的版本数。

"""


def pre_process(directory):
    # display all columns and rows, and set the item of row of dataframe
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.width', 5000)

    version_date_20201231_before = pd.DataFrame(columns=['project', 'name', 'project-name', 'releaseDate'],
                                                dtype=object)
    version_date_20201231_after = pd.DataFrame(columns=['project', 'name', 'project-name', 'releaseDate'], dtype=object)

    with open(directory + "List.txt") as l_all:
        lines_all = l_all.readlines()

    for line in lines_all:

        project = line.replace("\n", "")
        print("the file is ", project)
        # read versionDate, excluding the date of version after 2020.12.31

        version_date = pd.read_csv("F:\\MTmeta\\pyData\\versionDate\\" + project.upper() + '_versionDate.csv',
                                   index_col=False)
        # print(version_date.head())
        # print(version_date.columns)
        for i in range(len(version_date)):
            print(version_date.loc[i, 'releaseDate'])
            i_version_date = datetime.strptime(version_date.loc[i, 'Date'], "%Y-%m-%dT%H:%M:%S")
            expire_time = datetime.strptime('2020-12-31', "%Y-%m-%d")
            print(i_version_date, expire_time)
            if i_version_date > expire_time:
                print("Current version date exceeds the deadline!")
                version_date_20201231_after = version_date_20201231_after.append({'project': project,
                                                                                  'name': version_date.loc[i, 'name'],
                                                                                  'project-name': project + '-' +
                                                                                                  version_date.loc[
                                                                                                      i, 'name'],
                                                                                  'releaseDate': version_date.loc[
                                                                                      i, 'releaseDate']},
                                                                                 ignore_index=True)
            else:
                version_date_20201231_before = version_date_20201231_before.append({'project': project,
                                                                                    'name': version_date.loc[i, 'name'],
                                                                                    'project-name': project + '-' +
                                                                                                    version_date.loc[
                                                                                                        i, 'name'],
                                                                                    'releaseDate': version_date.loc[
                                                                                        i, 'releaseDate']},
                                                                                   ignore_index=True)
            # break
        for root, dirs, files in os.walk(directory + project):
            for name in files:

                if not os.path.exists('F:\\MTmeta\\pyData\\data_20201231_before\\' + project):
                    os.mkdir('F:\\MTmeta\\pyData\\data_20201231_before\\' + project)

                print(name)
                print(name[:-4])
                print(version_date_20201231_before['project-name'].values)
                if name[0:-4] in version_date_20201231_before['project-name'].values:
                    print("Current version date does not exceed the deadline!")

                    if os.path.exists('F:\\MTmeta\\pyData\\data_20201231_before\\' + project + '\\' + name):
                        print("the csv file is created in last execution, so it will not be created this time.")
                        continue

                    try:
                        shutil.copy(directory + project + '\\' + name,
                                    'F:\\MTmeta\\pyData\\data_20201231_before\\' + project + '\\' + name)
                    except IOError as e:
                        print("Unable to copy file. %s" % e)

                # break

        # break

    version_date_20201231_after.to_csv("F:\\MTmeta\\pyData\\version_date_20201231_after.csv", index=False)
    version_date_20201231_before.to_csv("F:\\MTmeta\\pyData\\version_date_20201231_before.csv", index=False)


def pre_process_defect(directory):
    # display all columns and rows, and set the item of row of dataframe
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.width', 5000)

    statistics_languages = pd.DataFrame(columns=['project', 'versions', 'bugs', 'bug_counts', 'h', 'java', 'cpp', 'cs',
                                                 'py', 'js', 'hpp', 'cc'], dtype=object)

    with open(directory + "List.txt") as l_all:
        lines_all = l_all.readlines()

    for line in lines_all:

        project = line.replace("\n", "")
        print("the file is ", project)

        versions, bugs, bug_counts, h, java, cpp, cs, py, js, hpp, cc = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

        for root, dirs, files in os.walk(directory + 'data_20201231_before\\' + project):
            for name in files:
                print(name)

                if not os.path.exists(directory + 'data_defects\\' + project):
                    os.mkdir(directory + 'data_defects\\' + project)

                df_name = pd.read_csv(directory + 'data_20201231_before\\' + project + '\\' + name, index_col=False)

                if df_name['bug'].sum() > 0:
                    print("The bugs in current version is not null!")

                    if os.path.exists(directory + 'data_defects\\' + project + '\\' + name):
                        print("the csv file is created in last execution, so it will not be created this time.")
                        continue

                    try:
                        shutil.copy(directory + 'data_20201231_before\\' + project + '\\' + name,
                                    directory + 'data_defects\\' + project + '\\' + name)
                    except IOError as e:
                        print("Unable to copy file. %s" % e)

                    versions += 1
                    bugs = bugs + df_name['bug'].sum()
                    bug_counts = bug_counts + (len(df_name) - len(df_name[df_name['bug'] == 0]))
                    for i in range(len(df_name)):
                        item = df_name.loc[i, 'relName']
                        if '\\' in item:
                            item = item.replace('\\', '/')
                        language = item.split('/')[-1].split('.')[-1]
                        # if item.split('/')[-1].split('.')[-1] not in language:
                        #     language.append(item.split('/')[-1].split('.')[-1])
                        if language == 'h':
                            h += 1
                        if language == 'java':
                            java += 1
                        if language == 'cpp':
                            cpp += 1
                        if language == 'cs':
                            cs += 1
                        if language == 'py':
                            py += 1
                        if language == 'js':
                            js += 1
                        if language == 'hpp':
                            hpp += 1
                        if language == 'cc':
                            cc += 1

                # break
        print(versions, bugs, bug_counts)
        statistics_languages = statistics_languages.append({'project': project, 'versions': versions, 'bugs': bugs,
                                                            'bug_counts': bug_counts, 'h': h, 'java': java, 'cpp': cpp,
                                                            'cs': cs, 'py': py, 'js': js, 'hpp': hpp, 'cc': cc},
                                                           ignore_index=True)
        # break
    statistics_languages = statistics_languages.append({'project': 'total',
                                                        'versions': statistics_languages['versions'].sum(),
                                                        'bugs': statistics_languages['bugs'].sum(),
                                                        'bug_counts': statistics_languages['bug_counts'].sum(),
                                                        'h': statistics_languages['h'].sum(),
                                                        'java': statistics_languages['java'].sum(),
                                                        'cpp': statistics_languages['cpp'].sum(),
                                                        'cs': statistics_languages['cs'].sum(),
                                                        'py': statistics_languages['py'].sum(),
                                                        'js': statistics_languages['js'].sum(),
                                                        'hpp': statistics_languages['hpp'].sum(),
                                                        'cc': statistics_languages['cc'].sum()}, ignore_index=True)
    statistics_languages.to_csv(directory + "languages_statistics.csv", index=False)


def pre_process_defect_java(directory):
    # display all columns and rows, and set the item of row of dataframe
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.width', 5000)
    statistics_languages_java = pd.DataFrame(columns=['project', 'name', 'project-name', 'bugs', 'bug_counts',
                                                      'bugs_java', 'bug_counts_java', 'h', 'java', 'cpp', 'cs', 'py',
                                                      'js', 'hpp', 'cc'], dtype=object)

    statistics_languages_total = pd.DataFrame(columns=['project', 'versions', 'bugs', 'bug_counts', 'java'],
                                              dtype=object)

    with open(directory + "List.txt") as l_all:
        lines_all = l_all.readlines()

    for line in lines_all:

        project = line.replace("\n", "")
        # if project != 'xercesj':
        #     continue
        print("the file is ", project)

        if not os.path.exists(directory + 'data_defects\\' + project):
            continue

        df_last_name = pd.read_csv("F:\\MTmeta\\pyData\\version_date_20201231_before.csv", index_col=False)
        df_last_name = df_last_name[df_last_name['project'] == project].reset_index()
        last_name = df_last_name.loc[len(df_last_name) - 1, 'project-name']

        versions, bugs_versions, bug_counts_versions, java_versions = 0, 0, 0, 0

        for root, dirs, files in os.walk(directory + 'data_defects\\' + project):
            for name in files:
                print(name)
                print(name[:-4])
                bugs, bug_counts, h, java, cpp, cs, py, js, hpp, cc = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                bugs_java, bug_counts_java = 0, 0
                versions += 1
                # if not os.path.exists(directory + 'data_defects_java\\' + project):
                # os.mkdir(directory + 'data_defects_java\\' + project)

                df_name = pd.read_csv(directory + 'data_defects\\' + project + '\\' + name, index_col=False)

                bugs = bugs + df_name['bug'].sum()
                bug_counts = bug_counts + (len(df_name) - len(df_name[df_name['bug'] == 0]))
                for i in range(len(df_name)):
                    item = df_name.loc[i, 'relName']
                    if '\\' in item:
                        item = item.replace('\\', '/')
                    language = item.split('/')[-1].split('.')[-1]
                    # if item.split('/')[-1].split('.')[-1] not in language:
                    #     language.append(item.split('/')[-1].split('.')[-1])
                    if language == 'h':
                        h += 1
                    if language == 'java':
                        java += 1
                        bugs_java = bugs_java + df_name.loc[i, 'bug']
                        if df_name.loc[i, 'bug'] > 0:
                            bug_counts_java += 1
                    if language == 'cpp':
                        cpp += 1
                    if language == 'cs':
                        cs += 1
                    if language == 'py':
                        py += 1
                    if language == 'js':
                        js += 1
                    if language == 'hpp':
                        hpp += 1
                    if language == 'cc':
                        cc += 1

                if java < 6:
                    continue

                if bugs_java == 0:
                    continue

                java_versions = java_versions + java
                bugs_versions = bugs_versions + bugs_java
                bug_counts_versions = bug_counts_versions + bug_counts_java

                if not os.path.exists(directory + 'data_defects_java\\' + project):
                    os.mkdir(directory + 'data_defects_java\\' + project)
                if not os.path.exists(directory + 'data_defects_java_training\\' + project):
                    os.mkdir(directory + 'data_defects_java_training\\' + project)
                if not os.path.exists(directory + 'data_defects_java_testing\\' + project):
                    os.mkdir(directory + 'data_defects_java_testing\\' + project)

                if os.path.exists(directory + 'data_defects_java_training\\' + project + '\\' + name):
                    print("the csv file is created in last execution, so it will not be created this time.")
                    continue

                try:
                    shutil.copy(directory + 'data_defects\\' + project + '\\' + name,
                                directory + 'data_defects_java\\' + project + '\\' + name)
                    if last_name == name[:-4]:
                        shutil.copy(directory + 'data_defects\\' + project + '\\' + name,
                                    directory + 'data_defects_java_testing\\' + project + '\\' + name)
                    else:
                        shutil.copy(directory + 'data_defects\\' + project + '\\' + name,
                                    directory + 'data_defects_java_training\\' + project + '\\' + name)
                except IOError as e:
                    print("Unable to copy file. %s" % e)

                statistics_languages_java = statistics_languages_java.append(
                    {'project': project, 'name': name, 'project-name': project + '-' + name, 'bugs': bugs,
                     'bug_counts': bug_counts, 'bugs_java': bugs_java, 'bug_counts_java': bug_counts_java, 'h': h,
                     'java': java, 'cpp': cpp, 'cs': cs, 'py': py, 'js': js, 'hpp': hpp, 'cc': cc}, ignore_index=True)
                # break

        if not os.path.exists(directory + 'data_defects_java_testing\\' + project + '\\' + last_name + '.csv'):
            for j in range(len(df_last_name)):
                last_name = df_last_name.loc[len(df_last_name) - 1 - j, 'project-name']
                if not os.path.exists(directory + 'data_defects_java_training\\' + project + '\\' + last_name + '.csv'):
                    continue
                else:
                    try:
                        shutil.copy(directory + 'data_defects_java_training\\' + project + '\\' + last_name + '.csv',
                                    directory + 'data_defects_java_testing\\' + project + '\\' + last_name + '.csv')
                    except IOError as e:
                        print("Unable to copy file. %s" % e)
                    os.remove(directory + 'data_defects_java_training\\' + project + '\\' + last_name + '.csv')
                    break

        statistics_languages_total = statistics_languages_total.append(
            {'project': project, 'versions': versions, 'bugs': bugs_versions, 'bug_counts': bug_counts_versions,
             'java': java_versions}, ignore_index=True)
        # break
    statistics_languages_java = statistics_languages_java.append({'project': 'total', 'name': '/', 'project-name': '/',
                                                                  'bugs': statistics_languages_java['bugs'].sum(),
                                                                  'bug_counts': statistics_languages_java[
                                                                      'bug_counts'].sum(),
                                                                  'bugs_java': statistics_languages_java[
                                                                      'bugs_java'].sum(),
                                                                  'bug_counts_java': statistics_languages_java[
                                                                      'bug_counts_java'].sum(),
                                                                  'h': statistics_languages_java['h'].sum(),
                                                                  'java': statistics_languages_java['java'].sum(),
                                                                  'cpp': statistics_languages_java['cpp'].sum(),
                                                                  'cs': statistics_languages_java['cs'].sum(),
                                                                  'py': statistics_languages_java['py'].sum(),
                                                                  'js': statistics_languages_java['js'].sum(),
                                                                  'hpp': statistics_languages_java['hpp'].sum(),
                                                                  'cc': statistics_languages_java['cc'].sum()},
                                                                 ignore_index=True)
    statistics_languages_java.to_csv(directory + "statistics_languages_java.csv", index=False)

    statistics_languages_total = statistics_languages_total.append(
        {'project': 'total', 'versions': statistics_languages_total['versions'].sum(),
         'bugs': statistics_languages_total['bugs'].sum(), 'bug_counts': statistics_languages_total['bug_counts'].sum(),
         'java': statistics_languages_total['java'].sum()}, ignore_index=True)
    statistics_languages_total.to_csv(directory + "statistics_languages_java_versions.csv", index=False)


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

    Directory = "F:\\metricData\\perlProcessed\\"
    result_Directory = "F:\\MTmeta\\pyData\\"
    os.chdir(Directory)

    # 1. 从3903个版本中筛选出3720个版本（日期在2020.12.31之前）
    # pre_process(Directory)

    # 2.从3626个版本中筛选去没有缺陷的版本。3409版本（72个项目）满足要求
    # pre_process_defect(result_Directory)

    # 3.从3409版本（72个项目）满足要求版本上，筛选出只有JAVA的版本,3276 (3188)个版本符合要求
    pre_process_defect_java(result_Directory)

    e_time = time.time()
    execution_time = e_time - s_time

    print("The __name__ is ", __name__, ".\nFrom ", time.asctime(time.localtime(s_time)), " to ",
          time.asctime(time.localtime(e_time)), ",\nThis", os.path.basename(sys.argv[0]), "ended within",
          execution_time, "(s), or ", (execution_time / 60), " (m).")
