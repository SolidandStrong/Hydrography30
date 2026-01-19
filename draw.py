# -*- coding: utf-8 -*-
"""
@Time ： 2025/4/1 19:00
@Auth ：
@File ：draw.py
@IDE ：PyCharm
"""
import csv
import os
from matplotlib import cm, colors
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import geopandas as gpd
# plt.rcParams["font.family"] = "Times New Roman"  # 设置全局字体
from pylab import mpl
from scipy.stats import gaussian_kde


import special_process

mpl.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams["font.family"] = "Times New Roman"  # 设置全局字体
def draw_NHD_relative_error():
    plt.rcParams["font.family"] = "Times New Roman"  # 设置全局字体
    NHD_FAV_CSI_result = r'F:\青藏高原水体数据集\DATA\valid\NHD\NHD_FAB_CSI_result.csv'
    with open(NHD_FAV_CSI_result,'r') as f:
        FAB_NHD_area_0_005 = [[],[]]
        FAB_NHD_area_005_01 = [[], []]
        FAB_NHD_area_01_02 = [[], []]
        FAB_NHD_area_02_99 = [[], []]
        reader = csv.reader(f)
        n = 0
        for i in reader:
            if n==0:
                n+=1
                continue
            infos = i[0].split('_')
            if float(infos[4]) <= 0.05:
                FAB_NHD_area_0_005[0].append(float((infos[2]))/1000)
                FAB_NHD_area_0_005[1].append(float((infos[1]))/1000)
            if 0.05 < float(infos[4]) <= 0.1:
                FAB_NHD_area_005_01[0].append(float((infos[2]))/1000)
                FAB_NHD_area_005_01[1].append(float((infos[1]))/1000)
            if 0.1 < float(infos[4]) <= 0.2:
                FAB_NHD_area_01_02[0].append(float((infos[2]))/1000)
                FAB_NHD_area_01_02[1].append(float((infos[1]))/1000)
            if 0.2 < float(infos[4]) :
                FAB_NHD_area_02_99[0].append(float((infos[2]))/1000)
                FAB_NHD_area_02_99[1].append(float((infos[1]))/1000)
        f.close()

    plt.scatter(FAB_NHD_area_02_99[1], FAB_NHD_area_02_99[0], alpha=0.5, c='green')
    plt.scatter(FAB_NHD_area_01_02[1],FAB_NHD_area_01_02[0],alpha=0.5,c = 'blue')
    plt.scatter(FAB_NHD_area_005_01[1], FAB_NHD_area_005_01[0], alpha=0.5, c='orange')
    plt.scatter(FAB_NHD_area_0_005[1], FAB_NHD_area_0_005[0], alpha=0.5, c='red')
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.plot([0,750],[0,750],color = 'grey')
    plt.show()


def binned_stat(x, y, bins=20):
    bins = pd.qcut(x, bins, duplicates='drop')
    g = pd.DataFrame({'x': x, 'y': y}).groupby(bins)
    return pd.DataFrame({
        'x_center': g['x'].median(),
        'y_median': g['y'].median(),
        'y_q25': g['y'].quantile(0.25),
        'y_q75': g['y'].quantile(0.75),
    })

def draw_Global_relative_error():
    plt.rcParams["font.family"] = "Times New Roman"  # 设置全局字体
    NHD_FAV_CSI_result = r'F:\青藏高原水体数据集\New_DATA\GRDC\制图\各大洲GRDC验证csv\Global_CSI.csv'#r'F:\青藏高原水体数据集\description\制图\全球CSI验证\CSI文件\Global_CSI.csv'

    with open(NHD_FAV_CSI_result,'r') as f:
        FAB_NHD_area_0_005 = [[],[]]
        FAB_NHD_area_005_01 = [[], []]
        FAB_NHD_area_01_02 = [[], []]
        FAB_NHD_area_02_99 = [[], []]
        reader = csv.reader(f)
        n = 0
        all_re = []
        all_CSI = []
        CSI_60 = []
        area = []
        for i in reader:
            if n==0:
                n+=1
                # continue
            infos = i.copy()
            # print(i)
            RE = abs(float(i[11]))#abs(float(i[6])-float(i[9]))/float(i[6])
            all_re.append(RE)
            all_CSI.append(float(i[10]))
            area.append(float(i[9]))
            if float(i[10]) >= 0.8:
                CSI_60.append(float(i[10]))
            if RE <= 0.05:
                FAB_NHD_area_0_005[0].append(float((infos[6]))/1)
                FAB_NHD_area_0_005[1].append(float((infos[9]))/1)
            if 0.05 < RE <= 0.1:
                FAB_NHD_area_005_01[0].append(float((infos[6]))/1)
                FAB_NHD_area_005_01[1].append(float((infos[9]))/1)
            if 0.1 < RE <= 0.2:
                FAB_NHD_area_01_02[0].append(float((infos[6]))/1)
                FAB_NHD_area_01_02[1].append(float((infos[9]))/1)
            if 0.2 < RE :
                FAB_NHD_area_02_99[0].append(float((infos[6]))/1)
                FAB_NHD_area_02_99[1].append(float((infos[9]))/1)
        f.close()



    print(np.mean(all_re),np.median(all_re))
    print(len(CSI_60),len(all_CSI),len(CSI_60)/len(all_CSI))
    print(np.mean(all_CSI),np.median(all_CSI))

    print(len(FAB_NHD_area_0_005[0]),len(all_re),len(FAB_NHD_area_0_005[0])/len(all_re))
    print(len(FAB_NHD_area_005_01[0]), len(all_re), len(FAB_NHD_area_005_01[0]) / len(all_re))
    print(len(FAB_NHD_area_01_02[0]), len(all_re), len(FAB_NHD_area_01_02[0]) / len(all_re))
    print((len(FAB_NHD_area_0_005[0])+len(FAB_NHD_area_005_01[0])+len(FAB_NHD_area_01_02[0]))/len(all_re))

    plt.yscale("log")
    plt.xscale("log")
    plt.scatter(FAB_NHD_area_02_99[1], FAB_NHD_area_02_99[0], alpha=0.3, c=color(40,120,181),s=42)
    plt.scatter(FAB_NHD_area_01_02[1],FAB_NHD_area_01_02[0],alpha=0.6,c = color(243,210,102),s=38)
    plt.scatter(FAB_NHD_area_005_01[1], FAB_NHD_area_005_01[0], alpha=0.4, c=color(147,148,231),s=34)
    plt.scatter(FAB_NHD_area_0_005[1], FAB_NHD_area_0_005[0], alpha=0.2, c=color(200,36,35),s=12)
    plt.xticks(fontsize = 12)
    # plt.xlabel("Hydrography30 (Km2)",fontsize = 13)
    # plt.ylabel("GRDC (Km2)",fontsize = 13)
    plt.yticks(fontsize = 12)
    plt.plot([0,5000000],[0,5000000],color = 'black',linewidth = 1)
    plt.savefig(r'F:\青藏高原水体数据集\New_DATA\GRDC\制图\GRDC_Re.svg',transparent=True)
    plt.show()
    plt.close()
    #
    #
    # plt.boxplot(all_re,showfliers=False,meanline=True,medianprops={'color':'black','linewidth':2},showmeans=True,
    #             meanprops={'color':'red','linewidth':2},patch_artist=True,boxprops=dict(facecolor=color(154,201,219),linewidth = 2),
    #             whiskerprops = dict(linewidth = 2),capprops = dict(linewidth = 2))
    # plt.axis('off')
    # plt.savefig(r'F:\青藏高原水体数据集\New_DATA\GRDC\制图\GRDC_Re_box.svg', transparent=True)
    # plt.show()
    #
    # plt.boxplot(all_CSI,showfliers=False,meanline=True,medianprops={'color':'black','linewidth':2},showmeans=True,
    #             meanprops={'color':'red','linewidth':2},patch_artist=True,boxprops=dict(facecolor=color(190,184,220),linewidth = 2),
    #             whiskerprops = dict(linewidth = 2),capprops = dict(linewidth = 2))
    # plt.axis('off')
    # plt.savefig(r'F:\青藏高原水体数据集\New_DATA\GRDC\制图\GRDC_CSI_box.svg', transparent=True)
    # plt.show()



def draw_NHD_Global_CSI_RE_errorbar():
    plt.rcParams["font.family"] = "Times New Roman"  # 设置全局字体

    lon_lat_area = {}
    with open(r'F:\青藏高原水体数据集\New_DATA\GRDC\制图\NHDPlus验证\NHD_AREA.csv','r') as f:
        reader = csv.reader(f)
        for i in reader:
            lon_lat_area.setdefault(i[2]+i[3],float(i[7]))
    # print(lon_lat_area)


    result = []
    NHD_FAB_CSI_result = r'F:\青藏高原水体数据集\New_DATA\GRDC\制图\NHDPlus验证\CSI1.csv'
    with open(NHD_FAB_CSI_result, 'r') as f:
        NHD_CSI_RE = [[],[]]
        FAB_NHD_area_0_005 = [[], []]
        FAB_NHD_area_005_01 = [[], []]
        FAB_NHD_area_01_02 = [[], []]
        FAB_NHD_area_02_99 = [[], []]
        all_re = []
        all_CSI = []
        CSI_60 = []
        area = []
        reader = csv.reader(f)
        for i in reader:


            n = 0


            baseArea = float(lon_lat_area[i[2]+i[3]])
            if baseArea <= 0:
                continue
            calarea = float(i[7])

            RE = abs(2*(calarea-baseArea)/(baseArea+calarea))
            # if RE >= 5:
            #     continue
            all_re.append(RE)
            all_CSI.append(float(i[8]))
            if float(i[8])>=0.6:
                CSI_60.append(float(i[8]))
            area.append(baseArea)
            # print(RE)

            result.append([float(i[2]),float(i[3]),float(i[8])])
            if RE <= 0.05:
                FAB_NHD_area_0_005[0].append(baseArea / 1)
                FAB_NHD_area_0_005[1].append(calarea / 1)

            if 0.05 < RE <= 0.1:
                FAB_NHD_area_005_01[0].append(baseArea / 1)
                FAB_NHD_area_005_01[1].append(calarea / 1)
                # print(RE)
            if 0.1 < RE <= 0.2:
                FAB_NHD_area_01_02[0].append(baseArea / 1)
                FAB_NHD_area_01_02[1].append(calarea / 1)
            if 0.2 < RE:
                FAB_NHD_area_02_99[0].append(baseArea / 1)
                FAB_NHD_area_02_99[1].append(calarea / 1)
        f.close()
    # print(FAB_NHD_area_02_99)
    print(len(all_re))
    print(np.mean(all_re), np.median(all_re))
    print(len(CSI_60), len(all_CSI), len(CSI_60) / len(all_CSI))
    print(np.mean(all_CSI), np.median(all_CSI))

    # plt.boxplot(all_re)
    #
    print(len(FAB_NHD_area_0_005[0]), len(all_re), len(FAB_NHD_area_0_005[0]) / len(all_re))
    print(len(FAB_NHD_area_005_01[0]), len(all_re), len(FAB_NHD_area_005_01[0]) / len(all_re))
    print(len(FAB_NHD_area_01_02[0]), len(all_re), len(FAB_NHD_area_01_02[0]) / len(all_re))
    print((len(FAB_NHD_area_0_005[0]) + len(FAB_NHD_area_005_01[0]) + len(FAB_NHD_area_01_02[0])) / len(all_re))

    # plt.yscale("log")
    # plt.xscale("log")
    # plt.scatter(FAB_NHD_area_02_99[1], FAB_NHD_area_02_99[0], alpha=0.3, c=color(40,120,181),s=42)
    # plt.scatter(FAB_NHD_area_01_02[1],FAB_NHD_area_01_02[0],alpha=0.6,c = color(243,210,102),s=38)
    # plt.scatter(FAB_NHD_area_005_01[1], FAB_NHD_area_005_01[0], alpha=0.4, c=color(147,148,231),s=34)
    # plt.scatter(FAB_NHD_area_0_005[1], FAB_NHD_area_0_005[0], alpha=0.2, c=color(200,36,35),s=12)
    # plt.xticks(fontsize = 12)
    # plt.xlabel("FAB-Hydro (Km2)",fontsize = 13)
    # plt.ylabel("NHDPlus (Km2)",fontsize = 13)
    # plt.yticks(fontsize = 12)
    # plt.plot([0,5000000],[0,5000000],color = 'black',linewidth = 1)
    # plt.savefig(r'F:\青藏高原水体数据集\New_DATA\GRDC\制图\NHD_Re.svg',transparent=True)
    # plt.show()
    # plt.close()

    # with open(r'F:\青藏高原水体数据集\New_DATA\GRDC\制图\NHDPlus验证\NHD_CSI_map.csv','w',newline='') as f:
    #     writer = csv.writer(f)
    #     writer.writerows(result)
    #     f.close()

    plt.boxplot(all_re, showfliers=False, meanline=True, medianprops={'color': 'black', 'linewidth': 2}, showmeans=True,
                meanprops={'color': 'red', 'linewidth': 2}, patch_artist=True,
                boxprops=dict(facecolor=color(154, 201, 219), linewidth=2),
                whiskerprops=dict(linewidth=2), capprops=dict(linewidth=2))
    plt.axis('off')
    plt.savefig(r'F:\青藏高原水体数据集\New_DATA\GRDC\制图\NHD_Re_box.svg', transparent=True)
    plt.show()
    plt.close()

    plt.boxplot(all_CSI, showfliers=False, meanline=True, medianprops={'color': 'black', 'linewidth': 2},
                showmeans=True,
                meanprops={'color': 'red', 'linewidth': 2}, patch_artist=True,
                boxprops=dict(facecolor=color(190, 184, 220), linewidth=2),
                whiskerprops=dict(linewidth=2), capprops=dict(linewidth=2))
    plt.axis('off')
    plt.savefig(r'F:\青藏高原水体数据集\New_DATA\GRDC\制图\NHD_CSI_box.svg', transparent=True)
    plt.show()

def draw_O():


    file = r'C:\Users\张斌\Desktop\阿扎.csv'

    T = []  #温度
    O = []
    t = []   # 溶解氧时间
    tT = [] # 温度时间
    n = 0
    with open(file,'r',encoding='UTF-8') as f:
        reader = csv.reader(f)
        for i in reader:

            if n == 0:
                n += 1
                continue
            if float(i[0]) > 0.3:
                tT.append(i[2].split(' ')[0])
                T.append(float(i[0]))

            if 10 <= float(i[1]) <= 20:
                t.append(i[2].split(' ')[0])
                O.append(float(i[1]))

        f.close()
    print(len(T),len(tT),len(t),len(O))
    plt.subplot(2,1,1)
    plt.plot(range(0,len(T)),T,c = 'black',linewidth = 1)
    # plt.plot(range(0,len(T),int(len(T)/1440)),T[0::int(len(T)/1440)])
    # plt.scatter(range(0,len(T),int(len(T)/1440)),T[0::int(len(T)/1440)],facecolor = 'none',edgecolors='black')
    plt.xticks(range(0, len(T), int(len(T) / 5)), tT[0::int(len(T) / 5)], fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylabel('温度(℃)', fontsize=16)
    plt.title('阿扎站水质监测图',fontsize = 20)

    plt.subplot(2,1,2)
    plt.plot(range(0, len(O)), O, c='black', linewidth=1)
    plt.xticks(range(0, len(O), int(len(O) / 5)), t[0::int(len(O) / 5)], fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylabel('溶解氧浓度(mg/L)', fontsize=16)



    plt.show()




# ---------------------------------------- 河流验证绘图 ----------------------------------------------------
def draw_stream_point(point_csv):
    """
    根据验证的河流点绘制图象
    :param point_csv:
    :return:
    """

    df = pd.read_csv(point_csv)

    dx = df['dx_m']
    dy = df['dy_m']
    ds = df['ds_m']

    Q1 = df["dx_m"].quantile(0.25)
    Q3 = df["dx_m"].quantile(0.75)
    IQR = Q3 - Q1
    lower = Q1 - 1.5 * IQR
    upper = Q3 + 1.5 * IQR
    dx_clean = df[(df["dx_m"] >= lower) & (df["dx_m"] <= upper)]['dx_m']

    Q1 = df["dy_m"].quantile(0.25)
    Q3 = df["dy_m"].quantile(0.75)
    IQR = Q3 - Q1
    lower = Q1 - 1.5 * IQR
    upper = Q3 + 1.5 * IQR
    dy_clean = df[(df["dy_m"] >= lower) & (df["dy_m"] <= upper)]['dy_m']

    Q1 = df["ds_m"].quantile(0.25)
    Q3 = df["ds_m"].quantile(0.75)
    IQR = Q3 - Q1
    lower = Q1 - 1.5 * IQR
    upper = Q3 + 1.5 * IQR
    ds_clean = df[(df["ds_m"] >= lower) & (df["ds_m"] <= upper)]['ds_m']



    # print(dx_clean)
    # sns.violinplot(data=[dx_clean,dy_clean,ds_clean],labels = ['lon','lat','ds'])
    # plt.xticks([0,1,2],['lon','lat','ds'])
    # plt.xlabel('Direction',fontsize = 14)
    # plt.ylabel('Bias (m)',fontsize = 14)


    # # plt.boxplot(dx,showfliers=False)
    # # plt.yscale('log')
    # plt.show()

    return ds_clean

def draw_stream_area_figure():
    """
    绘制五个河流验证区的violin

    :return:
    """
    lon_river = draw_stream_point(r'F:\青藏高原水体数据集\New_DATA\Basin\SWORD_valid_point_result.csv')   # 长江
    congo = draw_stream_point(r'F:\青藏高原水体数据集\New_DATA\Basin\Congo\Congo_valid_point_result.csv')
    danube = draw_stream_point(r'F:\青藏高原水体数据集\New_DATA\Basin\Danube\Danube_valid_point_result.csv')
    mississippi = draw_stream_point(r'F:\青藏高原水体数据集\New_DATA\Basin\Mississippi\Mississippi_valid_point_result.csv')
    amazon = draw_stream_point(r'F:\青藏高原水体数据集\New_DATA\Basin\Amazon\Amazon_valid_point_result.csv')
    # plt.figure(figsize=(5,5))
    # ['#05B9E2','#8983BF','#C76DA2','#F97270','#F97270']
    # sns.violinplot(data=[mississippi,amazon,congo,danube,lon_river],
    #                labels = ['Long River Basin','Congo Basin','Danube Basin','Mississippi Basin','Amazon Basin'],
    #                colors=['white','white','white','white','white'],
    #                linewidth= 0.8)
    # plt.boxplot([mississippi,amazon,congo,danube,lon_river],showfliers=False,showmeans=True)

    plt.boxplot([mississippi, amazon, congo, danube, lon_river], showfliers = False,  showmeans=True,meanline = True, medianprops = {'color': 'black', 'linewidth': 2},
                meanprops={'color':'red','linewidth':2,'linestyle':'--'},patch_artist=True,boxprops=dict(facecolor='white',linewidth = 2),
                whiskerprops = dict(linewidth = 2),capprops = dict(linewidth = 2))



    # plt.xticks([0,1,2,3,4],['Long River','Congo','Danube','Mississippi','Amazon'],fontsize = 18,rotation = 30)
    plt.yticks(range(0,201,50),range(0,201,50),fontsize = 18)
    # plt.xlabel('Direction',fontsize = 14)
    plt.ylabel('Bias (m)',fontsize = 18)

    # paras = ["EU", 'AU', 'AS', 'Siberia', 'NA', 'SA', 'AF']
    # colors = ['#F97270', '#BB9727', '#F97270', '#32B897', '#05B9E2', '#8983BF', '#C76DA2']

    # plt.axis('off')  # 完全隐藏 x/y 轴和刻度
    plt.gca().set_frame_on(False)
    # plt.boxplot(dx,showfliers=False)
    # plt.yscale('log')
    plt.savefig(r'F:\青藏高原水体数据集\New_DATA\Stream\制图\stream_valid.svg',transparent=True)
    plt.show()

def draw_endorheic_bar():
    """
    绘制柱状图.svg (羌塘流域总的内流数)，加入到gis pro里
    :return:
    """

    Zhang = 274  # 张国庆的
    Bin = 473   # 我们的
    MERIT = 294  # MERIT hydro的

    labels = ['','1','2']
    data = [Zhang,MERIT,Bin]
    plt.bar(labels,data,width= 0.5, color = ['#2878B5','#FF8884','red'])

    # 去掉图框
    # 去掉坐标轴
    plt.axis('off')  # 完全隐藏 x/y 轴和刻度
    plt.gca().set_frame_on(False)
    plt.savefig(r'F:\青藏高原水体数据集\New_DATA\endorheic_point\制图\bar.svg',transparent=True)

    plt.show()

def color(r,g,b):

    return r/255,g/255,b/255
def draw_endorheic_pie():



    # labels = ['Past endorheic points', 'New endorheic points']
    # sizes = [182,112]
    # plt.pie(sizes,colors=['#05B9E2','#8983BF'],labels=None, autopct=None)
    # plt.savefig(r'F:\青藏高原水体数据集\New_DATA\endorheic_point\制图\MERIT1.svg',transparent=True)
    # plt.close()
    # print(182/294,112/294)
    # plt.show()

    # labels = ['219', '254']
    # sizes = [219, 254]
    # plt.pie(sizes, labels=None, autopct=None, colors=['#05B9E2','#8983BF'])
    # plt.savefig(r'F:\青藏高原水体数据集\New_DATA\endorheic_point\制图\Bin1.svg',transparent=True)
    # plt.close()

    # 数量占比
    # labels = ['219', '254','219', '254','219', '254','219']
    # sizes = [462, 349,923, 215,337,367,1235]
    sizes = [463, 349, 924, 215, 337, 367, 1235]
    plt.pie(sizes, labels=None, autopct="%1.1f%%", colors=['#C2B1D9', '#FFEBAF','#87BFFF','#32B897','#A8D5BA','#F5B895','#F1D8A7'],wedgeprops={'linewidth': 0.8, 'edgecolor': 'black'},textprops={"fontsize":14})
    plt.savefig(r'F:\青藏高原水体数据集\New_DATA\endorheic_point\制图\global_endorheic.svg', transparent=True)
    plt.show()
    plt.close()

    # 面积占比
    paras = ["EU",'AU','AS','SI','NA','SA','AF']
    venu = r'F:\青藏高原水体数据集\New_DATA\endorheic_point\制图\内流终点csv'
    sizes = []
    for para in paras:
        file = os.path.join(venu,para+'_Endorheic_basin_information.csv')
        df = pd.read_csv(file)
        print(df['Area'].sum())
        sizes.append(df['Area'].sum())
    #
    labels = ['219', '254', '219', '254', '219', '254', '219']
    # sizes = [462, 349, 923, 215, 337, 367, 1235]
    plt.pie(sizes, labels=None, autopct="%1.1f%%",
            colors=['#C2B1D9', '#FFEBAF','#87BFFF','#32B897','#A8D5BA','#F5B895','#F1D8A7'],
            wedgeprops={'linewidth': 0.8, 'edgecolor': 'black'}, textprops={"fontsize": 14})
    plt.savefig(r'F:\青藏高原水体数据集\New_DATA\endorheic_point\制图\global_endorheic_area.svg', transparent=True)
    # plt.show()
    plt.close()

    # 全球大小流域的数量统计
    # df = pd.read_csv(r'F:\青藏高原水体数据集\New_DATA\endorheic_point\制图\Global_Endorheic_basin_information.csv')
    # labels = ['2191', '2541', '219', '254', '2119']
    # sizes = [len(df[df['Area']<=100]), len(df[df['Area']<=1000])-len(df[df['Area']<=100]),len(df[df['Area']<=10000])-len(df[df['Area']<=1000]),len(df[df['Area']<=100000])-len(df[df['Area']<=10000]),len(df[df['Area']>100000])]
    # print(sizes)
    # print(sum(sizes))
    # plt.bar(labels, sizes,width=0.5,color=[color(255,255,0), color(0,112,255), color(255,166,127), color(38,115,0), color(255,0,0)])
    # # plt.gca().invert_yaxis()
    #
    # plt.axis('off')  # 不显示坐标轴，适合导入 ArcGIS Pro
    # plt.savefig(r'F:\青藏高原水体数据集\New_DATA\endorheic_point\制图\global_endorheic_size.svg', transparent=True)
    # plt.show()

def statis_stream_order_length(venu):


    # a = r'F:\青藏高原水体数据集\New_DATA\Stream\制图\各大洲非内流区河网\AF.shp'
    #
    # gdf = gpd.read_file(a)
    #
    # orders = gdf[gdf['order_'] == 1]['length']
    # length = gdf['length'].astype(float)
    # print(length)
    # temp = []
    # for order in range(1, 7):
    #
    #         x = gdf[gdf['order_'] == order]['length']
    #
    #         temp.append(sum(x))
    # print(temp)
    # plt.bar(range(6),temp,color = '#4392BB')


    all = []
    for name in ['AR','NA','SA','EU','AF','SI','AS','AU']:
        file = os.path.join(venu,name+'.shp')
        gdf = gpd.read_file(file)
        temp = []
        for order in range(1,7):

            x = gdf[gdf['order_'] == order]['length']

            temp.append(sum(x))

        all.append(temp)
        print(temp)

    L1 = []
    L2 = []
    L3 = []
    L4 = []
    L5 = []
    L6 = []
    for i in range(len(all)):
        L1.append(all[i][0])
        L2.append(sum(all[i][:1])+all[i][1])
        L3.append(sum(all[i][:2])+all[i][2])
        L4.append(sum(all[i][:3])+all[i][3])
        L5.append(sum(all[i][:4])+all[i][4])
        L6.append(sum(all[i][:5])+all[i][5])


    plt.bar(range(len(all)), L6,color = '#0070FF')
    plt.bar(range(len(all)),L5,color = '#0070FF')
    plt.bar(range(len(all)),L4,color = '#0070FF')
    plt.bar(range(len(all)), L3,color = '#00A9E6')
    plt.bar(range(len(all)), L2,color = '#00C5FF')
    plt.bar(range(len(all)),L1,color = '#73DFFF')
    plt.yticks(range(0,400001,100000),range(0,400001,100000),fontsize = 16)
    plt.gca().set_frame_on(False)
    plt.savefig(r'F:\青藏高原水体数据集\New_DATA\Stream\制图\continent_stream_length.svg',transparent = True)


    # print(orders)

    plt.show()



if __name__ == '__main__':

    # draw_NHD_relative_error()


    # ---------------------------------------- 绘图代码 ----------------------------------------------------------
    # draw_Global_relative_error()
    # draw_NHD_Global_CSI_RE_errorbar()
    # draw_stream_area_figure()


    # draw_endorheic_bar()
    draw_endorheic_pie()

    # statis_stream_order_length(r'F:\青藏高原水体数据集\New_DATA\Stream\制图\各大洲非内流区河网')
    # ---------------------------------------- 绘图代码 ----------------------------------------------------------



    # special_process.Check_self_intersect(r'F:\青藏高原水体数据集\New_DATA\Basin\五大湖\mask')
    # draw_NHD_Global_CSI_RE_errorbar()
    # draw_O()

    pass