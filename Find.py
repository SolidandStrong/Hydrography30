# -*- coding: utf-8 -*-
"""
@Time ： 2024/10/17 12:34
@Auth ：
@File ：Find.py
@IDE ：PyCharm
"""
import os
import shutil
from db import *
from osgeo import gdal
import global_sink
from multiprocessing import Pool

def find(db_name,outvenu):
    """
    查询是否生成该文件
    :param db_name:
    :param outvenu:
    :return:
    """
    infos = query_data_byType(db_name, str(2))

    for i,info in enumerate(infos):
        print(info[0])
        outpath = os.path.join(outvenu,str(info[0])+'.tif')
        if os.path.exists(outpath):
            print("{:d} is save successfully".format(info[0]))
        else:
            print("!!!!!!!!!!!! {:d} !!!!!".format(info[0]))

def copy_dems(originVenu,outVenu):

    if not os.path.exists(outVenu):
        os.mkdir(outVenu)
        os.chmod(outVenu,0o777)

    demsV = os.listdir(originVenu)
    for demsVV in demsV:
        dems = os.listdir(os.path.join(originVenu,demsVV))

        for dem in dems:
            shutil.copyfile(os.path.join(originVenu,demsVV,dem),os.path.join(outVenu,dem))


def clip(input_file,mask_file,output_file):
    """
    矢量裁剪栅格，用于提取后的流域裁剪河流、DEM，制图（需要手动将栅格流域转成矢量）
    :param input_file: .tif
    :param mask_file:  .shp
    :param output_file: .tif
    :return:
    """
    ds = gdal.Warp(output_file, input_file, format='GTiff',
                   cutlineDSName=mask_file,
                   cropToCutline=True,
                   copyMetadata=True,
                   creationOptions=['COMPRESS=LZW', "TILED=True","BIGTIFF=YES"]
                   )

def sbatch_clip(DemPath,baseVenu,maskVenu,flag = "NorthAmerica"):
    if not os.path.exists(baseVenu):
        os.mkdir(baseVenu)
        os.chmod(baseVenu,0o777)
    masks = os.listdir(maskVenu)
    for mask in masks:
        name = mask.split('.')
        if name[-1] != 'shp':
            continue
        maskPath = os.path.join(maskVenu,mask)

        outVenu = os.path.join(baseVenu,flag+"_"+name[0])
        if not os.path.exists(outVenu):
            os.mkdir(outVenu)
            os.chmod(outVenu,0o777)
        outPath = os.path.join(outVenu,flag+"_"+name[0]+'_dem_.tif')
        clip(DemPath,maskPath,outPath)

def sbatch_cal_dir_sink(dirVenu):
    import sink
    params = []
    dirs = os.listdir(dirVenu)
    for dirName in dirs:
        temp = []
        dirPath1 = os.path.join(dirVenu,dirName)

        demName = dirName+"_burnDEM_.tif"
        demPath = os.path.join(dirPath1,demName)
        dirName1 = dirName+"_burndir_.tif"
        dirPath = os.path.join(dirPath1, dirName1)
        sinkName = dirName+"_burnsink_.tif"
        sinkPath = os.path.join(dirPath1, sinkName)
        if not os.path.exists(demPath):
            print(demPath,"not exists!")
            continue

        params.append([demPath,dirPath,sinkPath])
        # global_sink.cal_dir_sink(demPath,dirPath,sinkPath)
        try:
            sink.Cal_dir(demPath, dirPath)
        except Exception as e:
            # 捕获并输出错误信息
            print(dirPath,"计算失败")
            print(f"An error occurred: {e}")
        try:
            sink.Cal_sink(demPath, sinkPath)
        except Exception as e:
            # 捕获并输出错误信息
            print(sinkPath,"计算失败")
            print(f"An error occurred: {e}")
    # po = Pool(3)
    # for info in params:
    #
    #     po.apply_async(global_sink.cal_dir_sink,(info[0],info[1],info[2],))
    # po.close()
    # po.join()

def draw_tif(file_path,saveFigPath):
    import matplotlib.pyplot as plt
    import numpy as np
    from PIL import Image

    # 打开 .tif 文件

    tif_image = Image.open(file_path)

    # 转换为 numpy 数组
    data = np.array(tif_image)

    # 绘图
    # plt.figure(figsize=(8, 8))
    plt.imshow(data, cmap='gray')
    plt.colorbar()
    plt.title('TIF Visualization')
    plt.show()
    plt.savefig(saveFigPath)
    plt.close()

def sbatch_draw(Venu):
    dirsVenu = os.listdir(Venu)
    for dirVenu in dirsVenu:
        dirPath = os.path.join(Venu,dirVenu)
        input_file = os.path.join(dirPath,dirVenu+"_sink_.tif")
        output_Path = os.path.join(dirPath, dirVenu + "_sink_.png")

        draw_tif(input_file,output_Path)


if __name__ == '__main__':
    draw_tif(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\Global_20241114\sink\sink_N52W132_FABDEM_V1-0.tif')
    # find(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\20240923\sinkdb.db')

    # ruoyun1591
    # sink.Cal_dir(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020024600/NorthAmerica_7020024600_dem_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020024600/NorthAmerica_7020024600_dir_.tif")
    # sink.Cal_sink(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020024600/NorthAmerica_7020024600_dem_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020024600/NorthAmerica_7020024600_sink_.tif")

    # 张斌1592
    # sink.Cal_dir(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020038340/NorthAmerica_7020038340_dem_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020038340/NorthAmerica_7020038340_dir_.tif")
    # sink.Cal_sink(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020038340/NorthAmerica_7020038340_dem_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020038340/NorthAmerica_7020038340_sink_.tif")

    # 谢凡耀1595
    # sink.Cal_dir(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020046750/NorthAmerica_7020046750_dem_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020046750/NorthAmerica_7020046750_dir_.tif")
    # sink.Cal_sink(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020046750/NorthAmerica_7020046750_dem_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020046750/NorthAmerica_7020046750_sink_.tif")


    # 大卫1598
    # sink.Cal_dir(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020047840/NorthAmerica_7020047840_dem_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020047840/NorthAmerica_7020047840_dir_.tif")
    # sink.Cal_sink(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020047840/NorthAmerica_7020047840_dem_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020047840/NorthAmerica_7020047840_sink_.tif")


    ################## SouthAmerica ####################
    # cry
    sink.Cal_dir(
        "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4020050470/Asia_4020050470_burnDEM_.tif",
        "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4020050470/Asia_4020050470_burndir_.tif")
    sink.Cal_sink(
        "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4020050470/Asia_4020050470_burnDEM_.tif",
        "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4020050470/Asia_4020050470_burnsink_.tif")
    # xfy
    sink.Cal_dir(
        "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4020050290/Asia_4020050290_burnDEM_.tif",
        "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4020050290/Asia_4020050290_burndir_.tif")
    sink.Cal_sink(
        "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4020050290/Asia_4020050290_burnDEM_.tif",
        "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4020050290/Asia_4020050290_burnsink_.tif")
    "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4020050220/Asia_4020050220_burnDEM_.tif"
    sink.Cal_dir(
        "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4020050220/Asia_4020050220_burnDEM_.tif",
        "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4020050220/Asia_4020050220_burndir_.tif")
    sink.Cal_sink(
        "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4020050220/Asia_4020050220_burnDEM_.tif",
        "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4020050220/Asia_4020050220_burnsink_.tif")

    #
    # sbatch_clip("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/SouthAmerica.tif","/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/","/datanode05/zhangbin/TBP_Stream/DATA/mask/Sa/","SouthAmerica")
    pass