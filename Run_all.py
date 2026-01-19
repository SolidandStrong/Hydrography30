# -*- coding: utf-8 -*-
"""
@Time ： 2025/6/27 10:49
@Auth ：
@File ：Run_all.py
@IDE ：PyCharm
"""
import csv
import shutil
from collections import deque

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import check_counterflow
import Find
import Raster
import application_endorheic
import delineate_basins
import genral_functions
import heap_PF_D8
import main
import process_FAB_SRTM
import process_sink
import search_outlet_manually
import special_process
import Merge
import sink
import endorheic_process
import calculate_sink_area_volume
import postprocess
import os
from osgeo import gdal
from multiprocessing import Pool

import split_cal_acc
import stream_valid
import valid
from compare_NHD import extravt_basin_for_po


def preprocess():
    """
    准备运行数据的
    :return:
    """
    # 复制文件到一个文件夹下
    # main.copy_DEM_flat("/datanode05/zhangbin/FAB_hydrography/Global_FABDEM/","/datanode05/zhangbin/FAB_hydrography/initial_data/FABDEM_flatten")

    # Asia
    # ------------------------------------------------------ Asia ------------------------------------------------------------------------- #
    # Merge.merge_DEM("/datanode05/zhangbin/FAB_hydrography/initial_data/FABDEM_flatten",[56.6,56.9,151.9,0.1],
    #                 "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Asia.tif")
    # main.read_and_save_tif("/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Asia.tif","/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Asia.png")
    # ------------------------------------------------------ Asia ------------------------------------------------------------------------- #

    # # Europe
    # ------------------------------------------------------ Europe ------------------------------------------------------------------------- #
    # Merge.merge_DEM("/datanode05/zhangbin/FAB_hydrography/initial_data/FABDEM_flatten", [-25.54, 82.8, 70.55, 11.59],
    #                 "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Europe_1.tif")
    #
    # # 修正合并的欧洲里海区域
    # dems = ["/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Europe_1.tif","/datanode05/zhangbin/FAB_hydrography/initial_data/SRTM/modified_FABDEM.tif"]
    # # 构造 BuildVRTOptions（正确方式）
    # vrt_options = gdal.BuildVRTOptions(
    #     separate=False,
    #     resolution='highest',
    #     VRTNodata=-9999  # 设定 NoData 值，可根据你的影像情况修改
    # )
    #
    # vrt_path = "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/mergedem.vrt"
    # vrt_ds = gdal.BuildVRT(vrt_path,dems, options=vrt_options)
    # vrt_ds = None  # 保存并关闭
    #
    # # 输出 GeoTIFF + 压缩
    # gdal.Translate(
    #     "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Europe.tif",
    #     vrt_path,
    #     format="GTiff",
    #     creationOptions=[
    #         "COMPRESS=LZW",
    #         "TILED=YES",
    #         "BIGTIFF=IF_SAFER"
    #     ]
    # )
    # os.remove(vrt_path)
    # main.read_and_save_tif("/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Europe.tif",
    #                        "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Europe.png")
    # ------------------------------------------------------ Europe ------------------------------------------------------------------------- #

    # SouthAmerica
    # ------------------------------------------------------ SouthAmerica ------------------------------------------------------------------------- #
    # Merge.merge_DEM("/datanode05/zhangbin/FAB_hydrography/initial_data/FABDEM_flatten", [-93, 15.8, -31.3, -56.9],
    #                 "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/SouthAmerica.tif")
    # main.read_and_save_tif("/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/SouthAmerica.tif",
    #                        "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/SouthAmerica.png")
    # ------------------------------------------------------ SouthAmerica ------------------------------------------------------------------------- #

    # NorthAmerica
    # ------------------------------------------------------ NorthAmerica ------------------------------------------------------------------------- #
    Merge.merge_DEM("/datanode05/zhangbin/FAB_hydrography/initial_data/FABDEM_flatten", [-138.9, 63.7, -51.6, 4.4],
                    "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/NorthAmerica.tif")
    main.read_and_save_tif("/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/NorthAmerica.tif",
                           "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/NorthAmerica.png")
    # ------------------------------------------------------ NorthAmerica ------------------------------------------------------------------------- #

    # Africa
    # ------------------------------------------------------ Africa ------------------------------------------------------------------------- #
    Merge.merge_DEM("/datanode05/zhangbin/FAB_hydrography/initial_data/FABDEM_flatten", [-19.1, 38.5, 55.5, -35.8],
                    "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Africa.tif")
    main.read_and_save_tif("/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Africa.tif",
                           "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Africa.png")
    # ------------------------------------------------------ Africa ------------------------------------------------------------------------- #

    # Greenland
    # ------------------------------------------------------ Greenland ------------------------------------------------------------------------- #
    Merge.merge_DEM("/datanode05/zhangbin/FAB_hydrography/initial_data/FABDEM_flatten", [-74, 84.6, -10.3, 58.7],
                    "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Greenland.tif")
    main.read_and_save_tif("/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Greenland.tif",
                           "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Greenland.png")
    # ------------------------------------------------------ Greenland ------------------------------------------------------------------------- #

    # Siberia
    # ------------------------------------------------------ Siberia ------------------------------------------------------------------------- #
    Merge.merge_DEM("/datanode05/zhangbin/FAB_hydrography/initial_data/FABDEM_flatten",[57.95,82.26,180,44.5],
                    "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Siberia.tif")
    main.read_and_save_tif("/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Siberia.tif",
                           "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Siberia.png")
    # ------------------------------------------------------ Siberia ------------------------------------------------------------------------- #

    # Australia
    # ------------------------------------------------------ Australia ------------------------------------------------------------------------- #
    Merge.merge_DEM("/datanode05/zhangbin/FAB_hydrography/initial_data/FABDEM_flatten", [90, 25.3, 180, -56.1],
                    "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Australia.tif")
    main.read_and_save_tif("/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Australia.tif",
                           "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Australia.png")
    # ------------------------------------------------------ Australia ------------------------------------------------------------------------- #

    # Arctic
    # ------------------------------------------------------ Arctic ------------------------------------------------------------------------- #
    Merge.merge_DEM("/datanode05/zhangbin/FAB_hydrography/initial_data/FABDEM_flatten", [-180, 84.2, -60.09, 50.2],
                    "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Arctic.tif")
    main.read_and_save_tif("/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Arctic.tif",
                           "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Arctic.png")
    # ------------------------------------------------------ Arctic ------------------------------------------------------------------------- #

# ----------------------------------------------- burning -----------------------------------------------#
def clip_code():
    """
    批量按照矢量文件夹来裁剪DEM
    :return:
    """
    # Asia
    # special_process.Check_self_intersect("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Asia")
    # special_process.sbatch_clip_ModifiedDir("/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Asia.tif",
    #                                         "/datanode05/zhangbin/FAB_hydrography/run_data/Asia",
    #                                         "/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Asia")
    # main.Clip_OSM_GSWO("/datanode05/zhangbin/FAB_hydrography/run_data/Asia")


    # Europe
    # special_process.Check_self_intersect("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Europe")
    # special_process.sbatch_clip_ModifiedDir("/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Europe.tif", "/datanode05/zhangbin/FAB_hydrography/run_data/Europe",
    #                                         "/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Europe")
    # main.Clip_OSM_GSWO("/datanode05/zhangbin/FAB_hydrography/run_data/Europe")

    # SouthAmerica
    # special_process.Check_self_intersect("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/SouthAmerica")
    # special_process.sbatch_clip_ModifiedDir(
    #     "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/SouthAmerica.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica",
    #     "/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/SouthAmerica")
    # main.Clip_OSM_GSWO("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica")

    # Siberia
    # special_process.Check_self_intersect(
    #     "/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Siberia")
    # special_process.sbatch_clip_ModifiedDir(
    #     "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Siberia.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/run_data/Siberia",
    #     "/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Siberia")
    # main.Clip_OSM_GSWO("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia")

    # Africa
    # special_process.Check_self_intersect(
    #     "/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Africa")
    # special_process.sbatch_clip_ModifiedDir(
    #     "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Africa.tif",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa",
    #     "/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Africa")
    # main.Clip_OSM_GSWO("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa")
    #
    # # NorthAmerica
    # special_process.Check_self_intersect(
    #     "/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/NorthAmerica")
    # special_process.sbatch_clip_ModifiedDir(
    #     "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/NorthAmerica.tif",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica",
    #     "/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/NorthAmerica")
    # main.Clip_OSM_GSWO("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica")
    #
    # # Australia
    # special_process.Check_self_intersect(
    #     "/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Australia")
    # special_process.sbatch_clip_ModifiedDir(
    #     "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Australia.tif",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia",
    #     "/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Australia")
    # main.Clip_OSM_GSWO("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia")

    # Greenland
    special_process.Check_self_intersect(
        "/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Greenland")
    special_process.sbatch_clip_ModifiedDir(
        "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/GreenLand_AW3D.tif",
        "/datanode05/zhangbin/FAB_hydrography/run_data/Greenland",
        "/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Greenland")
    # main.Clip_OSM_GSWO("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland")

    # Arctic
    # special_process.Check_self_intersect(
    #     "/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Arctic")
    # special_process.sbatch_clip_ModifiedDir(
    #     "/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/Arctic.tif",
    #     "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic",
    #     "/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Arctic")
    # main.Clip_OSM_GSWO("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic")

def prepare_occurance():
    """
    准备Occurance.tif文件
    :return:
    """

    # main.sbatch_occ("/datanode05/zhangbin/FAB_hydrography/run_data/Asia")
    # main.sbatch_occ("/datanode05/zhangbin/FAB_hydrography/run_data/Europe")
    # main.sbatch_occ("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica")
    # main.sbatch_occ("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia")

    # main.sbatch_occ("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica")
    main.sbatch_occ("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia")
    # main.sbatch_occ("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa")
    # main.sbatch_occ("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic")
    # main.sbatch_occ("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland")

def run_burning():
    """
    使用OCC、GDW对DEM进行burning
    :return:
    """

    # main.sbatch_burning('/datanode05/zhangbin/FAB_hydrography/run_data/Asia',
    #                     '/datanode05/zhangbin/FAB_hydrography/initial_data/GlobalDam/Asia')
    # main.sbatch_burning('/datanode05/zhangbin/FAB_hydrography/run_data/Europe',
    #                     '/datanode05/zhangbin/FAB_hydrography/initial_data/GlobalDam/Europe')
    # main.sbatch_burning('/datanode05/zhangbin/FAB_hydrography/run_data/Siberia',
    #                     '/datanode05/zhangbin/FAB_hydrography/initial_data/GlobalDam/Siberia')
    # main.sbatch_burning('/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica',
    #                     '/datanode05/zhangbin/FAB_hydrography/initial_data/GlobalDam/SouthAmerica')

    # main.sbatch_burning("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa",
    #                     '/datanode05/zhangbin/FAB_hydrography/initial_data/GlobalDam/Africa')
    main.sbatch_burning('/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica',
                        '/datanode05/zhangbin/FAB_hydrography/initial_data/GlobalDam/NorthAmerica')
    # main.sbatch_burning('/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia',
    #                     '/datanode05/zhangbin/FAB_hydrography/initial_data/GlobalDam/Australia')
    # main.sbatch_burning('/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland',
    #                     '/datanode05/zhangbin/FAB_hydrography/initial_data/GlobalDam/Greenland')
    # main.sbatch_burning("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic",
    #                     '/datanode05/zhangbin/FAB_hydrography/initial_data/GlobalDam/Arctic')
# ----------------------------------------------- burning -----------------------------------------------#


# ----------------------------- calculate initial fdir and depressions ----------------------------------#
def run_fdir_drpression():
    """
    按照每个流域计算fdir和depression
    :return:
    """
    # Asia
    # Asia_venu = "/datanode05/zhangbin/FAB_hydrography/run_data/Asia/Asia_4020000010/Asia_4020000010"#"/datanode05/zhangbin/FAB_hydrography/run_data/Asia"#
    # for dirName in os.listdir(Asia_venu):
    #     try:
    #         sink.calculate_fdir_depression(os.path.join(Asia_venu,dirName))
    #     except Exception as e:
    #         print(dirName,'is fail:',e)
    #         continue

    # sink.calculate_fdir_depression("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/Asia_4030050410")
    # sink.calculate_fdir_depression("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/Asia_4030050270")
    # sink.calculate_fdir_depression("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/Asia_4030050240")
    # sink.calculate_fdir_depression("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/Asia_4030050230")
    # sink.calculate_fdir_depression("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/Asia_4030050220")
    # sink.calculate_fdir_depression("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/Asia_4030050210")
    # sink.calculate_fdir_depression("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/Asia_4020050470")
    # sink.calculate_fdir_depression("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/Asia_4020050290")
    # sink.calculate_fdir_depression("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/Asia_4020034510")
    # sink.calculate_fdir_depression("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/Asia_4020024190")
    # sink.calculate_fdir_depression("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/Asia_4020015090")
    # sink.calculate_fdir_depression("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/Asia_4020006940")
    # sink.calculate_fdir_depression("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/Asia_4020000010")

    # # Europe
    # Europe_venu = "/datanode05/zhangbin/FAB_hydrography/run_data/Europe"
    # for dirName in os.listdir(Europe_venu):
    #     try:
    #         sink.calculate_fdir_depression(os.path.join(Europe_venu,dirName))
    #     except Exception as e:
    #         print(dirName,'is fail:',e)
    # #
    #
    # # Siberia
    # Siberia_venu = "/datanode05/zhangbin/FAB_hydrography/run_data/Siberia"
    # for dirName in os.listdir(Siberia_venu):
    #     try:
    #         sink.calculate_fdir_depression(os.path.join(Siberia_venu, dirName))
    #     except Exception as e:
    #         print(dirName, 'is fail:', e)
    #
    # # SouthAmerica
    # SouthAmerica_venu = "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica"
    # for dirName in os.listdir(SouthAmerica_venu):
    #     try:
    #         sink.calculate_fdir_depression(os.path.join(SouthAmerica_venu, dirName))
    #     except Exception as e:
    #         print(dirName, 'is fail:', e)

    # # # Africa
    # SouthAmerica_venu = "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa"
    # for dirName in os.listdir(SouthAmerica_venu):
    #     try:
    #         sink.calculate_fdir_depression(os.path.join(SouthAmerica_venu, dirName))
    #     except Exception as e:
    #         print(dirName, 'is fail:', e)
    # #
    # # # # NorthAmerica
    # SouthAmerica_venu = "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica"
    # for dirName in os.listdir(SouthAmerica_venu):
    #     try:
    #         sink.calculate_fdir_depression(os.path.join(SouthAmerica_venu, dirName))
    #     except Exception as e:
    #         print(dirName, 'is fail:', e)
    #
    # # # Australia
    # SouthAmerica_venu = "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia"
    # for dirName in os.listdir(SouthAmerica_venu):
    #     try:
    #         sink.calculate_fdir_depression(os.path.join(SouthAmerica_venu, dirName))
    #     except Exception as e:
    #         print(dirName, 'is fail:', e)

    # # Arctic
    # SouthAmerica_venu = "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic"
    # for dirName in os.listdir(SouthAmerica_venu):
    #     try:
    #         sink.calculate_fdir_depression(os.path.join(SouthAmerica_venu, dirName))
    #     except Exception as e:
    #         print(dirName, 'is fail:', e)
    #
    # # # Greenland
    SouthAmerica_venu = "/datanode05/zhangbin/FAB_hydrography/run_data/Greenland"
    for dirName in os.listdir(SouthAmerica_venu):
        try:
            sink.calculate_fdir_depression(os.path.join(SouthAmerica_venu, dirName))
        except Exception as e:
            print(dirName, 'is fail:', e)

def run_clip_next_level():
    """
    按找shp裁剪下一层级的burnDEM，Occ
    :return:
    """

    # Asia
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Asia/20010",
    #                      "/datanode05/zhangbin/FAB_hydrography/run_data/Asia/Asia_4020000010")

    # Europe
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Europe/4230",
    #                      "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_2020024230")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Europe/6860",
    #                      "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_2030068680")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Europe/0090",
    #                      "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_2030030090")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Europe/6030",
    #                      "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_2030026030")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Europe/1070",
    #                      "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_2040031070")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Europe/3320",
    #                      "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_2050033320")

    # Australia
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Australia/5870",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/Australia_5020055870")

    # Africa
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Africa/0010",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/Africa_1020000010")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Africa/1530",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/Africa_1020011530")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Africa/4170",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/Africa_1020034170")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Africa/8110",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/Africa_1020018110")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Africa/4_1530",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/Africa_1030011530")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Africa/4_8110",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/Africa_1030008110")

    # Siberia
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/8670",
    #                      "/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/Siberia_3020008670")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/3790",
    #                      "/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/Siberia_3020003790")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Siberia/temp",
    #                      "/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/Siberia_3040411570")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Siberia/3880",
    #                      "/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/Siberia_3040483880")
    # NorthAmerica
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/NorthAmerica/7840",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/NorthAmerica_7020047840")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/NorthAmerica/4600",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/NorthAmerica_7020024600")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/NorthAmerica/1430",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/NorthAmerica_7020021430")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/NorthAmerica/0010",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/NorthAmerica_7020000010")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/NorthAmerica/4_0010",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/NorthAmerica_7030000010")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/NorthAmerica/3520",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/NorthAmerica_7030034520")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/NorthAmerica/4_4600",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/NorthAmerica_7030024600")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/NorthAmerica/9280",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/NorthAmerica_7030049280")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/NorthAmerica/2240",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/NorthAmerica_7030022240")

    # SouthAmerica
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/SouthAmerica/1870",
    #                      "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/SouthAmerica_6020021870")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/SouthAmerica/14330",
    #                      "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/SouthAmerica_6020014330")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/SouthAmerica/6540",
    #                      "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/SouthAmerica_6020006540")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/SouthAmerica/9280",
    #                      "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/SouthAmerica_6020029280")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/SouthAmerica/1871",
    #                      "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/SouthAmerica_6030021871")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/SouthAmerica/5990",
    #                      "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/SouthAmerica_6040285990")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/SouthAmerica/9040",
    #                      "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/SouthAmerica_6050029040")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/SouthAmerica/5560",
    #                      "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/SouthAmerica_6040015560")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/SouthAmerica/6_14330",
    #                      "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/SouthAmerica_6040014330")


    # Arctic
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Arctic/9560",
    #                      "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/Arctic_8030009560")
    # main.clip_next_level("/datanode05/zhangbin/FAB_hydrography/initial_data/clip_05_buffer/Arctic/6860",
    #                      "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/Arctic_8030016860")

    pass

# ----------------------------- calculate initial fdir and depressions ----------------------------------#


# ----------------------------- Modified_fdir ----------------------------------#

def run_calculate_sink_infos():
    """
    计算sink的属性，构建数据库
    :return:
    """

    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa")
    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/FAB_hydrography/run_data/Asia")
    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland")
    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/FAB_hydrography/run_data/Europe")
    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia")
    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica")
    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica")
    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic")
    calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia")

def run_special_endorheic_fdir():

    endorheic_process.process_special_endorheic_fdir("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/Asia_4030050230")
    endorheic_process.process_special_endorheic_fdir("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/Asia_4030050240")


def run_special_exorheic_fdir():
    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/SouthAmerica_6050491720")
    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/SouthAmerica_6060015530")
    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/SouthAmerica_6060015120")
    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/SouthAmerica_6060014910")
    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/SouthAmerica_6050016970")
    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/SouthAmerica_6050015950")
    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/SouthAmerica_6050015570")

    # process_sink.priority_D8_cal_dir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/Africa_1040837890")

    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_2050033320")

    # process_sink.priority_D8_cal_dir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/NorthAmerica_7030021430")
    # process_sink.priority_D8_cal_dir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/NorthAmerica_7040022240")
    # process_sink.priority_D8_cal_dir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/NorthAmerica_7040054660")
    # process_sink.priority_D8_cal_dir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/NorthAmerica_7040000010")
    # process_sink.priority_D8_cal_dir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/NorthAmerica_7040192560")
    # process_sink.priority_D8_cal_dir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/NorthAmerica_7040192420")
    # process_sink.priority_D8_cal_dir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/NorthAmerica_7040402430")
    # process_sink.priority_D8_cal_dir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/NorthAmerica_7040392520")

    # Siberia
    process_sink.priority_D8_cal_dir("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/Siberia_3040411570")
    process_sink.priority_D8_cal_dir("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/Siberia_3040441670")
    process_sink.priority_D8_cal_dir("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/Siberia_3040483880")
    process_sink.priority_D8_cal_dir("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/Siberia_3040483881")

    pass

def run_process_endorheic():
    """
    复制内流区文件到endorheic文件夹，作为目视解译的参考
    :return:
    """

    # process_sink.sbatch_main1("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland")
    process_sink.sbatch_main1("/datanode05/zhangbin/FAB_hydrography/run_data/Asia")
    # process_sink.sbatch_main1("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia")
    # process_sink.sbatch_main1("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic")
    # process_sink.sbatch_main1("/datanode05/zhangbin/FAB_hydrography/run_data/Europe")
    # process_sink.sbatch_main1("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica")
    # process_sink.sbatch_main1("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa")
    # process_sink.sbatch_main1("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica")
    # process_sink.sbatch_main1("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia")

# ----------------------------- Modified_fdir ----------------------------------#

def check_sink(venu):
    """
    检查文件夹是否存在sink，输出不存在的文件名
    :return:
    """

    names = os.listdir(venu)
    for name in names:
        if len(name.split('.')) != 1:
            continue
        if not os.path.exists(os.path.join(venu,name,name+'_burnsink.tif')):
            print(name)

def po_cal_acc():

    venu = "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica"
    names = os.listdir(venu)
    for name in names:
        if len(name.split('.')) != 1:
            continue
        if not os.path.exists(os.path.join(venu, name, name + '_ModifiedDir.tif')):
            continue
        sink.Cal_acc(os.path.join(venu, name, name + '_ModifiedDir.tif'),os.path.join(venu, name, name + '_acc.tif'))
        sink.Cal_Stream(os.path.join(venu, name, name + '_acc.tif'),os.path.join(venu, name, name + '_stream_2000.tif'),2000)

def run_process_exorheic():
    """
    处理外流流域
    :return:
    """

    # postprocess.sbatch_VisualCheck("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica")
    postprocess.sbatch_VisualCheck("/datanode05/zhangbin/FAB_hydrography/run_data/Asia")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/FAB_hydrography/run_data/Europe")
    # postprocess.sbatch_VisualCheck("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland")
    # postprocess.sbatch_VisualCheck("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic")
    # postprocess.sbatch_VisualCheck("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia")
    # postprocess.sbatch_VisualCheck("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica")


# ----------------------------- Modified_fdir ----------------------------------#

# ----------------------------- Merge_fdir ----------------------------------#
def run_clip_outDir():
    """
    用0.01°缓冲区修建呕吐DIR，得到clip_outDIR.tif
    :return:
    """
    # NA_clip_venu = "/datanode05/xiefy/zhangb/FAB_Hydro/clip_mask/NA"
    # clip_venu = "/datanode05/xiefy/zhangb/FAB_Hydro/clip_mask/temp"
    # venu = "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica"

    # clip_venu = "/datanode05/zhangbin/FAB_hydrography/initial_data/clip_mask_001/Asia/"
    # venu = "/datanode05/zhangbin/FAB_hydrography/run_data/Asia/"

    # clip_venu = "/datanode05/zhangbin/FAB_hydrography/initial_data/clip_mask_001/Europe/"
    # venu = "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/"

    # clip_venu = "/datanode05/tangjj/zhangb/FAB_Hydro/initial_data/clip_mask_001/Arctic/"
    # venu = "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/"

    clip_venu = "/datanode05/tangjj/zhangb/FAB_Hydro/initial_data/clip_mask_001/Greenland/"
    venu = "/datanode05/zhangbin/FAB_hydrography/run_data/Greenland/"

    # clip_venu = "/datanode05/xiefy/zhangb/FAB_Hydro/clip_mask/Africa/"
    # venu = "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/"

    # clip_venu = "/datanode05/zhangbin/FAB_hydrography/initial_data/clip_mask_001/Siberia/"
    # venu = "/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/"

    # clip_venu = "/datanode05/xiefy/zhangb/FAB_Hydro/clip_mask/Australia/"
    # venu = "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/"

    # clip_venu = "/datanode05/zhangbin/FAB_hydrography/initial_data/clip_mask_001/SouthAmerica/"
    # venu = "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/"


    for filename in os.listdir(venu):

        if len(filename.split('.')) != 1:
            continue

        baseName = filename.split('.')[0]

        infile = os.path.join(venu,baseName,baseName+'_ModifiedDir.tif')  #
        if not os.path.exists(infile):
            print(infile)
            infile = os.path.join(venu,baseName,baseName+'_outDIR_.tif')
        clip_file = os.path.join(clip_venu,baseName+'.shp')
        if not os.path.exists(clip_file):
            print(clip_file,'缺失')
            continue
        outfile = os.path.join(venu,baseName,baseName+'_ClipdDir.tif')

        # if os.path.exists(outfile):
        #     continue

        # print(infile,clip_file,outfile)
        # Find.clip(infile, clip_file, outfile)
        try:
            Find.clip(infile,clip_file,outfile)
        except Exception as e:
            print(baseName,e)
            continue

def merge(venu):

    lists = []
    for fileName in os.listdir(venu):
        if len(fileName.split('.')) > 1:
            continue
        clipfile = os.path.join(venu,fileName,fileName+'_ClipdDir.tif')
        #
        # A = Raster.get_raster(clipfile)
        # proj,geo,nodata = Raster.get_proj_geo_nodata(clipfile)
        # A[A == nodata] = 255
        # Raster.save_raster(clipfile,A,proj,geo,gdal.GDT_Byte,255)

        if os.path.exists(clipfile):
            lists.append(clipfile)

    # 构造 BuildVRTOptions（正确方式）
    vrt_options = gdal.BuildVRTOptions(
        separate=False,
        resolution='highest',
        VRTNodata=255  # 设定 NoData 值，可根据你的影像情况修改
    )
    # 构建 VRT
    vrt_path = os.path.join(venu, "temp_merged.vrt")
    vrt_ds = gdal.BuildVRT(vrt_path, lists, options=vrt_options)
    vrt_ds = None  # 保存并关闭

    # 输出 GeoTIFF + 压缩
    gdal.Translate(
        os.path.join(venu, "initial_dir.tif"),
        vrt_path,
        format="GTiff",
        creationOptions=[
            "COMPRESS=LZW",
            "TILED=YES",
            "BIGTIFF=IF_SAFER"
        ]
    )

def merge_files(lists,out_path):

    baseDir = os.path.dirname(out_path)

    # 构造 BuildVRTOptions（正确方式）
    vrt_options = gdal.BuildVRTOptions(
        separate=False,
        resolution='highest',
        VRTNodata=255  # 设定 NoData 值，可根据你的影像情况修改
    )
    # 构建 VRT
    vrt_path = os.path.join(baseDir,"temp_merged.vrt")
    vrt_ds = gdal.BuildVRT(vrt_path, lists, options=vrt_options)
    vrt_ds = None  # 保存并关闭

    # 输出 GeoTIFF + 压缩
    gdal.Translate(
        out_path,
        vrt_path,
        format="GTiff",
        creationOptions=[
            "COMPRESS=LZW",
            "TILED=YES",
            "BIGTIFF=IF_SAFER"
        ]
    )
    os.remove(vrt_path)


def run_merge():
    """
    合并clip_outDIR.tif
    :return:
    """
    # merge("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica")
    # Merge.merge_burnDEM("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica")

    # merge("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/")
    # Merge.merge_burnDEM("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/")

    merge("/datanode05/zhangbin/FAB_hydrography/run_data/Greenland/")
    # Merge.merge_burnDEM("/datanode05/zhangbin/FAB_hydrography/run_data/Greenland/")

    # merge("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/")
    # Merge.merge_burnDEM("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/")

    # merge("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/")
    # Merge.merge_burnDEM("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/")
    #
    # merge("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/")
    # Merge.merge_burnDEM("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/")

    # merge("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/")
    # Merge.merge_burnDEM("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/")

    # merge("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/")
    # Merge.merge_burnDEM("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/")
    #
    # merge("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/")
    # Merge.merge_burnDEM("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/")
# ----------------------------- Merge_fdir ----------------------------------#

# ----------------------------- process_counterflow_fdir ----------------------------------#
def check_huan(fdir_file,out_file):

    fdir = Raster.get_raster(fdir_file)
    proj,geo,nodata = Raster.get_proj_geo_nodata(fdir_file)
    result = check_counterflow.detect_cycles(fdir, nodata)
    Raster.save_raster(out_file,result,proj,geo,gdal.GDT_Byte,0)



def run_produce_counterflower():

    # flow = Raster.get_raster("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/initial_dir.tif")
    # proj,geo,nodata = Raster.get_proj_geo_nodata("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/initial_dir.tif")
    # result = check_counterflow.detect_cycles(flow,nodata)
    # df = pd.DataFrame(result)
    # df.to_csv("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/error_xy.csv")
    # A("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/split_fdir/region/","/datanode05/zhangbin/FAB_hydrography/run_data/Asia/split_fdir/counterflow/")

    # postprocess.counterflow_process("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/initial_dir.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/Asia/counterflow.tif")
    # postprocess.divide_conuterflow_region("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/counterflow.tif",
    #                                       "/datanode05/zhangbin/FAB_hydrography/run_data/Asia/BURN_DEM.tif",
    #                                       "/datanode05/zhangbin/FAB_hydrography/run_data/Asia/initial_dir.tif",
    #                                       "/datanode05/zhangbin/FAB_hydrography/run_data/Asia/counter_mask")
    # postprocess.sbatch_PFD8_exorheic("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/counter_mask/")
    # postprocess.check_all_0("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/counter_mask/")
    # postprocess.merge_counterflow_outdir("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/counter_mask/",
    #                                      "/datanode05/zhangbin/FAB_hydrography/run_data/Asia/initial_dir.tif")
    # check_huan("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/Final_FDIR.tif","/datanode05/zhangbin/FAB_hydrography/run_data/Asia/huan.tif")
    # check_counterflow.split_check("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/ACC/region_table/",
    #                               "/datanode05/zhangbin/FAB_hydrography/run_data/Asia/ACC/region_tree/",
    #                               "/datanode05/tangjj/zhangb/FAB_Hydro/check_csv.csv")
    # check_counterflow.C("/datanode05/tangjj/zhangb/FAB_Hydro/check_csv.csv")
    # check_counterflow.sbatch_D("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/ACC/region/")


    # postprocess.counterflow_process("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/initial_dir.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counterflow.tif")

    # postprocess.divide_conuterflow_region("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counterflow.tif",
    #                                       "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/BURN_DEM.tif",
    #                                       "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/initial_dir.tif",
    #                                       "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask")
    # postprocess.sbatch_PFD8_exorheic("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask")
    # postprocess.merge_counterflow_outdir("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask",
    #                                      "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/initial_dir.tif")
    # check_counterflow.split_check("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/ACC/region_table/",
    #                               "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/ACC/region_tree/",
    #                               "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/ACC/check_csv.csv")
    # # postprocess.check_all_0("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask")
    # check_counterflow.C("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/ACC/check_csv.csv")
    # check_counterflow.sbatch_D("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/ACC/region/")
    # check_counterflow.D("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/SouthAmerica_fdr.tif")
    # postprocess.sbatch_process_huan("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/SA_huan_global.xlsx","/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/_BURN_DEM.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/Final_FDIR.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/huan_global")
    # postprocess.sbatch_PFD8_exorheic("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/huan_global")
    # postprocess.check_all_0("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/huan_global")
    # special_process.PFD8_exorheic("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/huan_global/1/1_temp_dem.tif",
    #                               "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/huan_global/1/1_temp_sink.tif",
    #                               "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/huan_global/1/1_temp_fdr.tif",
    #                               "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/huan_global/1/outdir_1_.tif")
    # postprocess.counterflow_get_upstream("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/SA_huan_global.xlsx","/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_fdr.tif")
    # postprocess.divide_conuterflow_region("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/huan_counterflow.tif",
    #                                       "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/_BURN_DEM.tif",
    #                                       "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_fdr.tif",
    #                                       "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask1")
    # postprocess.sbatch_PFD8_exorheic("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask1")

    # special_process.PFD8_exorheic("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask1/1/dem_temp1.tif",
    #                             "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask1/1/sink_temp1.tif",
    #                             "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask1/1/fdr_temp1.tif",
    #                             "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask1/1/outfdr1.tif")
    # Find.clip("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask1/3/dir_3_.tif",
    #           "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask1/3/2.shp",
    #           "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask1/3/2_fdr.tif")
    #
    # Find.clip("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask1/3/dem_3_.tif",
    #           "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask1/3/2.shp",
    #           "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask1/3/2_dem.tif")
    #
    # Find.clip("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask1/3/sink_3_.tif",
    #           "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask1/3/2.shp",
    #           "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask1/3/2_sink.tif")
    # special_process.PFD8_exorheic(
    #     "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask1/3/2_dem.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask1/3/2_sink.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask1/3/2_fdr.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask1/3/outfdr3.tif")
    # special_process.PFD8_exorheic(
    #     r"F:\青藏高原水体数据集\3\SA\3_1\2_dem.tif",
    #     r"F:\青藏高原水体数据集\3\SA\3_1\2_sink.tif",
    #     r"F:\青藏高原水体数据集\3\SA\3_1\2_fdr.tif",
    #     r"F:\青藏高原水体数据集\3\SA\3_1\outfdr3.tif")
    # Find.clip("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final1/region/SouthAmerica_fdr.tif",
    #           "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final1/special/SA_补.shp",
    #           "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final1/special/SA_counter_fdr.tif")
    # check_counterflow.D("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final1/special/SA_counter_fdr.tif")
    # special_process.read_and_save_tif("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final1/region/SouthAmerica_fdr.tif",
    #                                   "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final1/special/SA_fdr.png")
    # Find.clip("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/SouthAmerica_6020006540.1/SouthAmerica_6020006540_burnDEM.tif",
    #           "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final1/special/SA_补.shp",
    #           "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final1/special/SA_补_dem.tif")
    # special_process.Fill_fdr_nodata(r'F:\青藏高原水体数据集\3\SA\特殊处理\SA_counter_fdr.tif',
    #                                 r'F:\青藏高原水体数据集\3\SA\SA_补_dem_fdr.tif')

    # special_process.PFD8_exorheic(
    #         "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final1/special/PF/sink_dem.tif",
    #         "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final1/special/PF/sink.tif",
    #         "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final1/special/PF/sink_fdr.tif",
    #         "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final1/special/PF/outdir.tif")



    # Find.clip("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/Final_FDIR.tif",
    #           "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/huan_global/1/1.shp",
    #           "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/huan_global/1/1_temp_fdr.tif")

    # A = Raster.get_raster(r'F:\青藏高原水体数据集\3\43\44_fdr1.tif')
    # proj,geo,nodata = Raster.get_proj_geo_nodata(r'F:\青藏高原水体数据集\3\43\44_fdr1.tif')
    # row,col = A.shape
    # for i in range(row):
    #     A[i,1] = 4
    # Raster.save_raster(r'F:\青藏高原水体数据集\3\43\44_fdr_补.tif',A,proj,geo,gdal.GDT_Byte,255)
    # postprocess.check_all_0("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/huan_global")
    # postprocess.merge_counterflow_outdir("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/counter_mask1",
    #                                      "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_fdr.tif",
    #                                      "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_fdr_2.tif")

    # A = Raster.get_raster("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_fdr_2.tif")
    # proj,geo,nodata = Raster.get_proj_geo_nodata("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_fdr_2.tif")
    # paras = []
    # with open("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/nodata_position.csv",'r') as f:
    #     reader = csv.reader(f)
    #     n = 0
    #     for i in reader:
    #         if n == 0:
    #             n += 1
    #             continue
    #         paras.append([float(i[2]),float(i[3])])
    #     f.close()
    # for para in paras:
    #     row,col = valid.lonlat_to_pixel(geo,para[0]+geo[1]/3,para[1]+geo[5]/3)
    #     A[row,col] = nodata
    #
    # Raster.save_raster("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_fdr_3.tif",A,proj,geo,gdal.GDT_Byte,nodata)
    # split_cal_acc.merge_update_acc("","/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_final_fdr.tif")




    # check_counterflow.D("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_fdr.tif")
    postprocess.sbatch_process_huan("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/SA_huan_global1.xlsx",
                                    "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/_BURN_DEM.tif",
                                    "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/SouthAmerica_fdr.tif",
                                    "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/huan_global1")


    print('----------------------------------------------------------------------------------------------------')
    #
    # postprocess.counterflow_process("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/initial_dir.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/counterflow.tif")
    # postprocess.divide_conuterflow_region("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/counterflow.tif",
    #                                       "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/BURN_DEM.tif",
    #                                       "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/initial_dir.tif",
    #                                       "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/counter_mask")
    # postprocess.sbatch_PFD8_exorheic("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/counter_mask")
    # postprocess.check_all_0("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/counter_mask")
    # postprocess.merge_counterflow_outdir("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/counter_mask",
    #                                      "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/initial_dir.tif")
    # check_counterflow.split_check("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/ACC/region_table/",
    #                               "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/ACC/region_tree/",
    #                               "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/ACC/check_csv.csv")
    # check_counterflow.C("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/ACC/check_csv.csv")
    # check_counterflow.sbatch_D("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/ACC/region/")
    check_counterflow.D("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Final_FDIR.tif")
    # postprocess.counterflow_process("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/counter_mask/15/dir_15_.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/counter_mask/15/counterflow_15_.tif")
    # postprocess.counterflow_process("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/counter_mask/2/dir_2_.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/counter_mask/2/counterflow_2_.tif")
    # postprocess.sbatch_process_huan("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/EU_huan_global.xlsx","/datanode05/zhangbin/FAB_hydrography/run_data/Europe/_BURN_DEM.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Final_FDIR.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/huan_global")
    # postprocess.sbatch_PFD8_exorheic("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/huan_global")
    # postprocess.check_all_0("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/huan_global")
    # Find.clip("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Final_FDIR.tif",
    #           "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/huan_global/233/233.shp",
    #           "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/huan_global/233/233_fdr1.tif")
    # Find.clip("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/huan_global/3/dir_3_.tif",
    #           "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/huan_global/3/3.shp",
    #           "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/huan_global/3/temp_dir_3_.tif")
    # Find.clip("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/huan_global/3/dem_3_.tif",
    #           "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/huan_global/3/3.shp",
    #           "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/huan_global/3/temp_dem_3_.tif")
    # Find.clip("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/huan_global/3/sink_3_.tif",
    #           "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/huan_global/3/3.shp",
    #           "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/huan_global/3/temp_sink_3_.tif")
    # special_process.PFD8_exorheic("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/huan_global/3/temp_dem_3_.tif",
    #                               "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/huan_global/3/temp_sink_3_.tif",
    #                               "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/huan_global/3/temp_dir_3_.tif",
    #                               "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/huan_global/3/outdir_3_.tif")
    # special_process.PFD8_exorheic(r"F:\青藏高原水体数据集\3\EU\3\temp_dem_3_.tif",
    #                               r"F:\青藏高原水体数据集\3\EU\3\temp_sink_3_.tif",
    #                               r"F:\青藏高原水体数据集\3\EU\3\temp_dir_3_.tif",
    #                               r"F:\青藏高原水体数据集\3\EU\3\outdir_3_.tif")
    # A = Raster.get_raster(r'F:\青藏高原水体数据集\3\EU\233_fdr1.tif')
    # proj,geo,nodata = Raster.get_proj_geo_nodata(r'F:\青藏高原水体数据集\3\EU\233_fdr1.tif')
    # row,col = A.shape
    # # for i in range(row):
    # #     A[i,1] = 4
    # A[0,0] = 64
    # A[1, 0] = 64
    # Raster.save_raster(r'F:\青藏高原水体数据集\3\EU\233_fdr1_补.tif',A,proj,geo,gdal.GDT_Byte,255)
    # check_counterflow.D("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_ACC_final/Europe_fdr.tif")
    # postprocess.sbatch_process_huan("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/EU_huan_global1.xlsx",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/_BURN_DEM.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_ACC_final/Europe_fdr.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/huan_global1")








    #
    # postprocess.counterflow_process("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/initial_dir.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/counterflow.tif")
    # postprocess.divide_conuterflow_region("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/counterflow.tif",
    #                                       "/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/BURN_DEM.tif",
    #                                       "/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/initial_dir.tif",
    #                                       "/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/counter_mask")
    # postprocess.sbatch_PFD8_exorheic("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/counter_mask")
    # postprocess.check_all_0("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/counter_mask")
    # postprocess.merge_counterflow_outdir("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/counter_mask",
    #                                      "/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/initial_dir.tif")
    # print('------------------------------------------------- Siberia -------------------------------------------------------------------')
    # check_counterflow.split_check("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/ACC/region_table/",
    #                               "/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/ACC/region_tree/",
    #                               "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Siberia_check_csv.csv")
    # check_counterflow.C("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Siberia_check_csv.csv")
    # check_counterflow.sbatch_D("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/ACC/region/")
    # check_counterflow.D("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/Final_FDIR.tif")
    # postprocess.sbatch_process_huan("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/SI_huan_global.xlsx",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/_BURN_DEM.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/Final_FDIR.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/huan_global")
    # postprocess.sbatch_PFD8_exorheic("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/huan_global")
    # postprocess.merge_counterflow_outdir("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/huan_global",
    #                                      "/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/Final_FDIR.tif",
    #                                      "/datanode05/zhangbin/FAB_hydrography/Result/Siberia_fdr.tif")


    # postprocess.counterflow_process("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/initial_dir.tif",
    #                                 "//datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/counterflow.tif")
    # postprocess.divide_conuterflow_region("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/counterflow.tif",
    #                                       "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/BURN_DEM.tif",
    #                                       "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/initial_dir.tif",
    #                                       "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/counter_mask")
    # postprocess.sbatch_PFD8_exorheic("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/counter_mask")
    # postprocess.check_all_0("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/counter_mask")
    # postprocess.merge_counterflow_outdir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/counter_mask",
    #                                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/initial_dir.tif")
    # postprocess.counterflow_process_check("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/Final_FDIR.tif")
    # check_huan("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/Final_FDIR.tif","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/huan.tif")
    # check_counterflow.split_check("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/ACC/region_table/",
    #                               "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/ACC/region_tree/",
    #                               "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/ACC/check_csv.csv")
    # check_counterflow.C("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/ACC/check_csv.csv")
    # check_counterflow.sbatch_D("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/ACC/region/")
    # postprocess.sbatch_process_huan("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/AF_huan.xlsx","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/huan_fdr")
    # check_counterflow.D("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/Final_FDIR.tif")

    # postprocess.sbatch_process_huan("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/AF_huan_global.xlsx","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/_BURN_DEM.tif",
    #                                 "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/Final_FDIR.tif","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/huan_global")
    # postprocess.sbatch_PFD8_exorheic("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/huan_global")
    # postprocess.check_all_0("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/huan_global")
    postprocess.merge_counterflow_outdir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/huan_global",
                                         "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/Final_FDIR.tif",
                                         "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/FDIR_process1.tif")
    # fdr = Raster.get_raster("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/FDIR_process1.tif")
    # proj,geo,nodata = Raster.get_proj_geo_nodata("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/FDIR_process1.tif")
    # fdr[29171, 162665] = 0    # 目视解译到此点为内流终点，手动修正
    # Raster.save_raster("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_fdr.tif",fdr,proj,geo,gdal.GDT_Byte,nodata)




    # postprocess.counterflow_process("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/initial_dir.tif",
    #                                 "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/counterflow.tif")
    # postprocess.divide_conuterflow_region("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/counterflow.tif",
    #                                       "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/BURN_DEM.tif",
    #                                       "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/initial_dir.tif",
    #                                       "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/counter_mask")
    # postprocess.check_all_0("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/counter_mask")
    # postprocess.merge_counterflow_outdir("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/counter_mask/",
    #                                      "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/initial_dir.tif")
    # check_counterflow.D("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/Final_FDIR.tif")
    # postprocess.counterflow_process("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/counter_mask/19/dir_19_.tif",
    #                                 "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/counter_mask/19/counterflow_19_.tif")
    # postprocess.counterflow_process("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/counter_mask/36/dir_36_.tif",
    #                                 "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/counter_mask/36/counterflow_36_.tif")

    # postprocess.counterflow_process("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/counter_mask/36/dir_36_.tif",
    #                                 "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/counter_mask/36/counterflow36.tif")
    # postprocess.counterflow_process("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/counter_mask/36/dir_36_.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/Asia/counterflow36.tif")
    # postprocess.sbatch_PFD8_exorheic("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/counter_mask")
    # postprocess.merge_counterflow_outdir("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/counter_mask",
    #                                      "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/initial_dir.tif")
    # postprocess.counterflow_process_check("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/Final_FDIR.tif")
    # check_huan("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/Final_FDIR.tif","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/huan.tif")
    # check_counterflow.split_check("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/ACC/region_table/",
    #                               "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/ACC/region_tree/",
    #                               "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/ACC/check_csv.csv")

    # check_counterflow.C("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/ACC/check_csv.csv")
    # check_counterflow.sbatch_D("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/ACC/region/")
    # check_counterflow.D("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/Final_FDIR.tif")

    # postprocess.sbatch_process_huan("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/AR_huan_global.xlsx","/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/_BURN_DEM.tif",
    #                                 "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/Final_FDIR.tif",
    #                                 "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/huan_global")
    # postprocess.sbatch_PFD8_exorheic("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/huan_global")
    # postprocess.check_all_0("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/huan_global")
    # Find.clip("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/Final_FDIR.tif",
    #           "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/huan_global/54/54.shp",
    #           "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/huan_global/54/54_fdr.tif")
    # special_process.PFD8_exorheic("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/huan_global/1/dem_1_.tif",
    #                               "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/huan_global/1/sink_1_.tif",
    #                               "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/huan_global/1/dir_1_.tif",
    #                               "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/huan_global/1/outdir_1_.tif",-135.6474676,68.7559798,64)
    # special_process.PFD8_exorheic("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/huan_global/54/dem_54_.tif",
    #                               "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/huan_global/54/sink_54_.tif",
    #                               "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/huan_global/54/dir_54_.tif",
    #                               "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/huan_global/54/outdir_54_.tif",-115.0928235,61.3897567,16)
    # A = Raster.get_raster(r'F:\青藏高原水体数据集\3\AR\3_fdr1.tif')
    # proj,geo,nodata = Raster.get_proj_geo_nodata(r'F:\青藏高原水体数据集\3\AR\3_fdr1.tif')
    # row,col = A.shape
    #
    # for i in range(col):
    #     A[0,i] = 16
    # Raster.save_raster(r'F:\青藏高原水体数据集\3\AR\3_fdr1_补.tif',A,proj,geo,gdal.GDT_Byte,255)
    # postprocess.merge_counterflow_outdir("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/huan_global",
    #                                      "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/Final_FDIR.tif",
    #                                      "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Arctic_fdr.tif")




    # postprocess.counterflow_process("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/initial_dir.tif",
    #                                 "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/counterflow.tif")
    # postprocess.divide_conuterflow_region("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/counterflow.tif",
    #                                       "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/BURN_DEM.tif",
    #                                       "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/initial_dir.tif",
    #                                       "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/counter_mask")
    # postprocess.sbatch_PFD8_exorheic("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/counter_mask")
    # postprocess.check_all_0("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/counter_mask")
    # postprocess.merge_counterflow_outdir("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/counter_mask",
    #                                      "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/initial_dir.tif")
    # check_counterflow.split_check("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/ACC/region_table/",
    #                               "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/ACC/region_tree/",
    #                               "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/ACC/check_csv.csv")
    # check_counterflow.C("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/ACC/check_csv.csv")
    # check_counterflow.sbatch_D("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/ACC/region/")
    # postprocess.sbatch_process_huan("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/GR_huan_global.xlsx","/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/_BURN_DEM.tif",
    #                                 "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/Final_FDIR.tif",
    #                                 "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/huan_global")
    # postprocess.sbatch_PFD8_exorheic("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/huan_global")
    # postprocess.merge_counterflow_outdir("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/huan_global",
    #                                      "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Greenland/Final_FDIR.tif",
    #                                      "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Greenland_fdr.tif")







    # postprocess.counterflow_process("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/initial_dir.tif",
    #                                 "//datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/counterflow.tif")
    # postprocess.divide_conuterflow_region("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/counterflow.tif",
    #                                       "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/BURN_DEM.tif",
    #                                       "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/initial_dir.tif",
    #                                       "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/counter_mask")
    # postprocess.sbatch_PFD8_exorheic("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/counter_mask")
    # postprocess.check_all_0("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/counter_mask")
    # postprocess.merge_counterflow_outdir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/counter_mask",
    #                                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/initial_dir.tif")
    # print('------------------------------------------------- NorthAmerica -------------------------------------------------------------------')
    # check_counterflow.split_check("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/ACC/region_table/",
    #                               "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/ACC/region_tree/",
    #                               "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/ACC/check_csv.csv")
    # check_counterflow.C("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/ACC/check_csv.csv")
    # check_counterflow.sbatch_D("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/ACC/region/")
    # check_counterflow.D("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/Final_FDIR.tif")
    # postprocess.sbatch_process_huan("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/NA_huan_global.xlsx",
    #                                 "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/_BURN_DEM.tif",
    #                                 "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/Final_FDIR.tif",
    #                                 "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/huan_global")
    # postprocess.sbatch_PFD8_exorheic("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/huan_global")
    # postprocess.merge_counterflow_outdir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/huan_global",
    #                                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/Final_FDIR.tif",
    #                                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_fdr.tif")



    # #
    # postprocess.counterflow_process("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/initial_dir.tif",
    #                                 "//datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/counterflow.tif")
    # postprocess.divide_conuterflow_region("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/counterflow.tif",
    #                                       "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/BURN_DEM.tif",
    #                                       "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/initial_dir.tif",
    #                                       "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/counter_mask")
    # postprocess.sbatch_PFD8_exorheic("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/counter_mask")
    # postprocess.check_all_0("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/counter_mask")
    # postprocess.merge_counterflow_outdir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/counter_mask",
    #                                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/initial_dir.tif")
    # check_counterflow.split_check("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/ACC/region_table/",
    #                               "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/ACC/region_tree/",
    #                               "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/ACC/check_csv.csv")
    # check_counterflow.C("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/ACC/check_csv.csv")
    # check_counterflow.sbatch_D("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/ACC/region/")
    # postprocess.sbatch_process_huan("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/AU_huan_global.xlsx","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/_BURN_DEM.tif",
    #                                 "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/Final_FDIR.tif",
    #                                 "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/huan_global")
    # postprocess.sbatch_PFD8_exorheic("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/huan_global")
    # postprocess.check_all_0("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/huan_global")
    # postprocess.merge_counterflow_outdir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/huan_global",
    #                                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/Final_FDIR.tif",
    #                                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Australia_fdr.tif")

    # postprocess.merge_all_processed_fdr("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/huan_fdr/",
    #                                     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/ACC/region/",
    #                                     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Australia_fdr.tif")
    # check_counterflow.D("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/Final_FDIR.tif")



    def check_out_dir(venu):
        counter_venu = os.path.join(venu,'counter_mask')
        for i in os.listdir(counter_venu):
            if os.path.exists(os.path.join(counter_venu,i,'outdir_'+str(i)+'_.tif')):
                continue
            print(os.path.join(counter_venu,i,'outdir_'+str(i)+'_.tif'))

    # check_out_dir("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/")
    # check_out_dir("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/")
    # check_out_dir("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/")
    # check_out_dir("/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/")
    # check_out_dir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/")
    # check_out_dir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/")
    # check_out_dir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/")
    # check_out_dir("/datanode05/zhangbin/FAB_hydrography/run_data/Asia/")

    # postprocess.F2("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/counter_mask/1/sink_1_.tif",
    #                "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/counter_mask/1/dir_1_.tif",
    #                 "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/counter_mask/1/dem_1_.tif",
    #                 "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/counter_mask/1/outdir_1_.tif")
    # print('Over')
    # postprocess.F2("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/counter_mask/90/sink_90_.tif",
    #                "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/counter_mask/90/dir_90_.tif",
    #                "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/counter_mask/90/dem_90_.tif",
    #                "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/counter_mask/90/outdir_90_.tif")


# ----------------------------- process_counterflow_fdir ----------------------------------#

# ----------------------------- process_calculate_acc_stream ----------------------------------#

def run_acc_calculation():
    """
    计算汇流累积量
    :return:
    """



    # split_cal_acc.cal_normal_acc_work_flow("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_fdr2.tif",
    #                                        "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final")
    # #
    # split_cal_acc.cal_normal_acc_work_flow("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_final_fdr.tif",
    #                                        "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final")

    # split_cal_acc.cal_normal_acc_work_flow("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_fdr2.tif",
    #                                        "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final")

    # split_cal_acc.cal_normal_acc_work_flow("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_fdr4.tif",
    #                                        "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final")


    # split_cal_acc.cal_normal_acc_work_flow("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_fdr.tif",
    #                                        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final")
    #
    # split_cal_acc.cal_normal_acc_work_flow("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_fdr2.tif",
    #                                        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final")
    #
    # split_cal_acc.cal_normal_acc_work_flow("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_fdr.tif",
    #                                        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final")


    # split_cal_acc.cal_normal_acc_work_flow("/datanode05/zhangbin/FAB_hydrography/run_data/Greenland/initial_dir.tif",
    #                                        "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Greenland/Greenland_ACC_final")

    # split_cal_acc.cal_normal_acc_work_flow("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_fdr3.tif",
    #                                        "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_ACC_final")

    # split_cal_acc.cal_normal_acc_work_flow(
    #     "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic_north_island/Arctic_north_island_fdr.tif",
    #     "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic_north_island/Arctic_north_island_ACC_final")


    # 亚欧大陆综合计算汇流累积量
    # split_cal_acc.cal_normal_acc_work_flow(
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia_Europe_Siberia/Asia_Europe_Siberia_fdr.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia_Europe_Siberia/ACC_final")

    split_cal_acc.cal_normal_acc_work_flow(
        "/datanode05/zhangbin/Nandita_watershed_delineation/Canadian_river/rundata/Canadian_fdr.tif",
        "/datanode05/zhangbin/Nandita_watershed_delineation/Canadian_river/rundata/Canadian_ACC_final")

def run_extract_stream():
    """
    计算河网
    :return:
    """

    # Merge.merge_update_fdr("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/region/",
    #                        "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/Siberia_fdr.tif")
    # sink.Cal_Stream("/datanode05/zhangbin/FAB_hydrography/Result/Siberia_ACC_final/FINAL_ACC.tif",
    #                 "/datanode05/zhangbin/FAB_hydrography/Result/Siberia_ACC_final/Siberia_stream_10000w.tif", 11111111)
    # sink.streamFeature(
    #     "/datanode05/zhangbin/FAB_hydrography/Result/Siberia_ACC_final/Siberia_fdr.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/Result/Siberia_ACC_final/Siberia_stream_10000w.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/Result/Siberia_ACC_final/Siberia_stream_10000w.shp")
    #

    # Merge.merge_update_fdr("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/region/",
    #                        "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/Asia_fdr.tif")
    # sink.Cal_Stream("/datanode05/zhangbin/FAB_hydrography/Result/Asia_ACC_final1/FINAL_ACC.tif",
    #                 "/datanode05/zhangbin/FAB_hydrography/Result/Asia_ACC_final1/Asia_stream_10000km2.tif",11111111)
    # sink.streamFeature("/datanode05/zhangbin/FAB_hydrography/Result/Asia_ACC_final1/Asia_fdr.tif",
    #                    "/datanode05/zhangbin/FAB_hydrography/Result/Asia_ACC_final1/Asia_stream_10000km2.tif",
    #                    "/datanode05/zhangbin/FAB_hydrography/Result/Asia_ACC_final1/Asia_stream_10000km2.shp")
    #
    # sink.Cal_Stream("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/FINAL_ACC.tif",
    #                 "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/SouthAmerica_stream_10000km2.tif",11111111)
    # sink.streamFeature("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/SouthAmerica_fdr.tif",
    #                    "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/SouthAmerica_stream_10000km2.tif",
    #                    "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/SouthAmerica_stream_10000km2.shp")

    # Merge.merge_update_fdr("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_ACC_final/region",
    #                        "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_ACC_final/Arctic_fdr.tif")
    # sink.Cal_Stream("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Arctic_ACC_final/FINAL_ACC.tif",
    #                 "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Arctic_ACC_final/Arctic_stream_10000km2.tif",11111111)
    # sink.streamFeature("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Arctic_ACC_final/Arctic_fdr.tif",
    #                    "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Arctic_ACC_final/Arctic_stream_10000km2.tif",
    #                    "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Arctic_ACC_final/Arctic_stream_10000km2.shp")

    # Merge.merge_update_fdr("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic_north_island/Arctic_north_island_ACC_final/region/",
    #                        "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic_north_island/Arctic_north_island_ACC_final/Arctic_north_island_fdr.tif")

    # Merge.merge_update_fdr("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/region/",
    #                        "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/Europe_fdr.tif")
    # sink.Cal_Stream("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_ACC_final/FINAL_ACC.tif",
    #                 "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_ACC_final/Europe_stream_10000km2.tif",
    #                 11111111)
    # A = Raster.get_raster("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_ACC_final/Europe_fdr.tif")
    # print(A.shape)
    # B = Raster.get_raster("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_ACC_final/Europe_stream_10000km2.tif")
    # print(B.shape)

    # sink.streamFeature("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_ACC_final/Europe_fdr.tif",
    #                    "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_ACC_final/Europe_stream_10000km2.tif",
    #                    "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_ACC_final/Europe_stream_10000km2.shp")


    # Merge.merge_update_fdr("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/region",
    #                        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/Africa_fdr.tif")

    # sink.Cal_Stream("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/FINAL_ACC.tif",
    #                 "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/Africa_stream_10000km2.tif",11111111)
    # sink.streamFeature("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/Africa_fdr.tif",
    #                    "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/Africa_stream_10000km2.tif",
    #                    "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/Africa_stream_10000km2.shp")
    #
    # Merge.merge_update_fdr("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/region/",
    #                        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/NorthAmerica_fdr.tif")
    # sink.Cal_Stream("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final2/FINAL_ACC.tif",
    #                 "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final2/NorthAmerica_stream_10000km2.tif", 1111111)
    # sink.streamFeature("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final2/NorthAmerica_fdr.tif",
    #                    "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final2/NorthAmerica_stream_10000km2.tif",
    #                    "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final2/NorthAmerica_stream_10000km2.shp")

    # Merge.merge_update_fdr("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/region",
    #                        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/Australia_fdr.tif")

    # sink.Cal_Stream("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Australia_ACC_final1/FINAL_ACC.tif",
    #                 "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Australia_ACC_final1/Australia_stream_10000km2.tif", 11111111)
    # sink.streamFeature("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr1/Australia_fdr.tif",
    #                    "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr1/Australia_ACC_final/Australia_stream_10000km2.tif",
    #                    "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr1/Australia_ACC_final/Australia_stream_10000km2.shp")

    # Merge.merge_update_fdr("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Greenland/Greenland_ACC_final/region/",
    #                        "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Greenland/Greenland_ACC_final/Greenland_fdr.tif")
    # sink.Cal_Stream("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Greenland_ACC_final/FINAL_ACC.tif",
    #                 "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Greenland_ACC_final/Greenland_stream_10km2.tif", 11111)
    # sink.streamFeature("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Greenland_ACC_final/Greenland_fdr.tif","/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Greenland_ACC_final/Greenland_stream_10km2.tif",
    #                    "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Greenland_ACC_final/Greenland_stream_10km2.shp")

    # sink.Cal_Stream("/datanode05/zhangbin/Analysis_3North/DATA/FABDEM/acc.tif",
    #                 "/datanode05/zhangbin/Analysis_3North/DATA/FABDEM/stream.tif",8888000)
    # sink.streamFeature("/datanode05/zhangbin/Analysis_3North/DATA/FABDEM/fdr.tif","/datanode05/zhangbin/Analysis_3North/DATA/FABDEM/stream.tif",
    #                    "/datanode05/zhangbin/Analysis_3North/DATA/FABDEM/stream_link.tif")

    # Merge.merge_update_fdr("/datanode05/zhangbin/Nandita_watershed_delineation/Canadian_river/rundata/Canadian_ACC_final/region/",
    #                        "/datanode05/zhangbin/Nandita_watershed_delineation/Canadian_river/rundata/Canadian_ACC_final/Canadian_fdr.tif")

    # Merge.merge_update_fdr(
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia_Europe_Siberia/ACC_final/region/",
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia_Europe_Siberia/ACC_final/Asia_Europe_Siberia_fdr.tif")

    pass

def run_extract_stream_valid_basin():


    def get_basin(fdr_file,lon,lat,out_file):
        fdr = Raster.get_raster(fdr_file)
        proj,geo,f_nodata = Raster.get_proj_geo_nodata(fdr_file)
        row, col = fdr.shape

        out_row,out_col = valid.lonlat_to_pixel(geo,lon,lat)
        extravt_basin_for_po(fdr, out_row, out_col, row, col, geo, proj, out_file)


    # 北美的密西西比河流域
    # get_basin("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final1.tif/NorthAmerica_fdr.tif",
    #           -89.3512464, 29.2693055, "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Mississippi_basin1.tif")
    # valid.raster2Feature("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Mississippi_basin1.tif",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Mississippi_basin1.shp")
    # Find.clip("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final1.tif/FINAL_ACC.tif",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Mississippi_basin1.shp",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Mississippi_basin_acc.tif")
    # Find.clip("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final1.tif/NorthAmerica_fdr.tif",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Mississippi_basin1.shp",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Mississippi_basin_fdr.tif")
    # sink.Cal_Stream("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Mississippi_basin_acc.tif",
    #                 "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Mississippi_basin_streams_1000km2.tif",
    #                 1111111)
    # sink.streamFeature("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Mississippi_basin_fdr.tif",
    #                    "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Mississippi_basin_streams_1000km2.tif",
    #                    "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Mississippi_basin_streams_1000km2.shp")



    # 欧洲的多瑙河流域
    # get_basin("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_ACC_final/Europe_fdr.tif",
    #           29.6710296, 45.3456181, "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Danube/Danube_basin.tif")
    # valid.raster2Feature("/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Danube/Danube_basin.tif",
    #                      "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Danube/Danube_basin.shp")

    # Find.clip("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_ACC_final/FINAL_ACC.tif",
    #           "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Danube/Danube_basin.shp",
    #           "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Danube/Danube_basin_acc.tif")
    # Find.clip("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_ACC_final/Europe_fdr.tif",
    #           "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Danube/Danube_basin.shp",
    #           "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Danube/Danube_basin_fdr.tif")
    # sink.Cal_Stream("/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Danube/Danube_basin_acc.tif",
    #                 "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Danube/Danube_basin_streams_1000km2.tif",
    #                 1111111)
    # sink.streamFeature("/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Danube/Danube_basin_fdr.tif",
    #                    "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Danube/Danube_basin_streams_1000km2.tif",
    #                    "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Danube/Danube_basin_streams_1000km2.shp")
    sink.stream_class("/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Danube/Danube_basin_streams_1000km2.tif",
                      "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Danube/Danube_basin_fdr.tif")



    # 非洲的刚果流域
    # get_basin("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/Africa_fdr.tif",
    #           13.1047297, -5.8897298, "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_basin1.tif")
    # valid.raster2Feature("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_basin1.tif",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_basin.shp")
    # Find.clip("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/FINAL_ACC.tif",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_basin.shp",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_basin_acc.tif")
    # Find.clip("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/Africa_fdr.tif",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_basin.shp",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_basin_fdr.tif")
    # sink.Cal_Stream("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_basin_acc.tif",
    #                 "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_basin_streams_1000km2.tif",
    #                 1111111)
    # sink.streamFeature("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_basin_fdr.tif",
    #                    "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_basin_streams_1000km2.tif",
    #                    "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_basin_streams_1000km2.shp")
    # A = Raster.get_raster("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_basin_fdr.tif")
    # print(A.shape)
    # #
    # B = Raster.get_raster("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_basin_streams_1000km2.tif")
    # # proj,geo,nodata = Raster.get_proj_geo_nodata("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_basin_streams_1000km2.tif")
    # print(B.shape)
    # B = B[:,:-1]
    # Raster.save_raster("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_basin_streams_1000km2.tif",B,proj,geo,gdal.GDT_Byte,nodata)
    # sink.streamFeature("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_basin_fdr.tif",
    #                    "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_basin_streams_1000km2.tif",
    #                    "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_basin_streams_1000km2.shp")





    # 南美的Amazon流域
    # get_basin("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/SouthAmerica_fdr.tif",
    #           -50.9914011, -0.0389006,"/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Amazon_basin.tif")
    # valid.raster2Feature("/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Amazon_basin.tif",
    #                      "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Amazon/Amazon_basin.shp")
    # Find.clip("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/FINAL_ACC.tif",
    #           "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Amazon/Amazon_basin.shp",
    #           "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Amazon/Amazon_basin_acc.tif")
    # Find.clip("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/SouthAmerica_fdr.tif",
    #           "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Amazon/Amazon_basin.shp",
    #           "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Amazon/Amazon_basin_fdr.tif")
    # sink.Cal_Stream("/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Amazon/Amazon_basin_acc.tif",
    #                 "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Amazon/Amazon_basin_streams_1000km2.tif",
    #                 1111111)
    # sink.streamFeature("/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Amazon/Amazon_basin_fdr.tif",
    #                    "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Amazon/Amazon_basin_streams_1000km2.tif",
    #                    "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Amazon/Amazon_basin_streams_1000km2.shp")
    sink.stream_class("/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Amazon/Amazon_basin_streams_1000km2.tif",
                      "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Amazon/Amazon_basin_fdr.tif")

    # 亚洲的长江流域
    # get_basin("/datanode05/zhangbin/FAB_hydrography/Result/Asia_ACC_final/Asia_fdr.tif",
    #           121.0907447, 31.7616748,"/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/LongRiver_basin.tif")
    # valid.raster2Feature("/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/LongRiver_basin.tif",
    #                      "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/LongRiver_basin.shp")
    # Find.clip("/datanode05/zhangbin/FAB_hydrography/Result/Asia_ACC_final/FINAL_ACC.tif",
    #           "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/LongRiver_basin.shp",
    #           "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/LongRiver_basin_acc.tif")
    # Find.clip("/datanode05/zhangbin/FAB_hydrography/Result/Asia_ACC_final/Asia_fdr.tif",
    #           "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/LongRiver_basin.shp",
    #           "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/LongRiver_basin_fdr.tif")
    # sink.Cal_Stream("/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/LongRiver_basin_acc.tif",
    #                 "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/LongRiver_basin_streams_1000km2.tif",
    #                 1111111)
    # sink.streamFeature( "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/LongRiver_basin_fdr.tif",
    #                    "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/LongRiver_basin_streams_1000km2.tif",
    #                    "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/LongRiver_basin_streams_1000km2.shp")
    sink.stream_class("/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/LongRiver_basin_streams_1000km2.tif",
                      "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/LongRiver_basin_fdr.tif")




    # ---------------------------------北美流向修正流域-----------------------------------------
    # get_basin("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final.tif/NorthAmerica_fdr.tif",
    #           -93.9381891, 47.5138299, "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/midify_basin_south.tif")
    # get_basin("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final.tif/NorthAmerica_fdr.tif",
    #           -92.7423585, 56.9328957,"/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/midify_basin_north.tif")
    # merge_files(["/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/midify_basin_north.tif","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/midify_basin_south.tif"],
    #             "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/mask.tif")
    # valid.raster2Feature("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/mask.tif",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/mask.shp")
    # Find.clip("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final.tif/NorthAmerica_fdr.tif",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/mask.shp","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/mask_fdr.tif")
    # Find.clip("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/_BURN_DEM.tif",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/mask.shp",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/mask_dem.tif")
    # special_process.PFD8_exorheic(
    #         "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/mask_dem.tif",
    #         "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/mask.tif",
    #         "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/mask_fdr.tif",
    #         "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/outdir_NA.tif",-92.7424151,56.9329588,128)

    # get_basin("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final.tif/NorthAmerica_fdr.tif",
    #           -109.9945862, 50.9043072,"/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/midify_basin_tri.tif")
    # valid.raster2Feature("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/midify_basin_tri.tif",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/midify_basin_tri.shp")
    # Find.clip("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final.tif/NorthAmerica_fdr.tif",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/midify_basin_tri.shp","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/mask_fdr1.tif")
    # Find.clip("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/_BURN_DEM.tif",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/midify_basin_tri.shp",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/mask_dem1.tif")
    # special_process.PFD8_exorheic(
    #         "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/mask_dem1.tif",
    #         "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/midify_basin_tri.tif",
    #         "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/mask_fdr1.tif",
    #         "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/outdir_NA_tri.tif",-109.9945862, 50.9043072,1)

    # get_basin("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final.tif/NorthAmerica_fdr.tif",
    #           -80.4677589, 51.3279026, "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/basin1.tif")
    # get_basin("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final.tif/NorthAmerica_fdr.tif",
    #           -82.0548704, 52.1312459,"/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/basin2.tif")
    # get_basin("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final.tif/NorthAmerica_fdr.tif",
    #           -73.1551590, 46.0351532,"/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/basin3.tif")
    # get_basin("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final.tif/NorthAmerica_fdr.tif",
    #           -85.6029169,49.5395870,"/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/basin4.tif")
    # get_basin("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final.tif/NorthAmerica_fdr.tif",
    #           -78.1535106,48.3276834,"/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/basin5.tif")
    # valid.raster2Feature("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/basin1.tif",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/basin1.shp")
    # valid.raster2Feature("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/basin2.tif",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/basin2.shp")
    # valid.raster2Feature("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/basin3.tif",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/basin3.shp")
    # valid.raster2Feature("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/basin4.tif",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/basin4.shp")
    # valid.raster2Feature("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/basin5.tif",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/basin5.shp")

    # Find.clip("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final1.tif/FINAL_ACC.tif",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/yousahng.shp",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/yousahng_upa.tif")

    # Find.clip("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final.tif/NorthAmerica_fdr.tif",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/mask_big.shp","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/fdr_big.tif")
    # Find.clip("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/_BURN_DEM.tif",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/mask_big.shp",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/dem_big.tif")
    #
    # Find.clip("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final.tif/NorthAmerica_fdr.tif",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/small/basin5.shp",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/small/fdr_5.tif")
    # Find.clip("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/_BURN_DEM.tif",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/small/basin5.shp",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/small/dem_5.tif")
    # merge_files(["/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_modified_FABDEM.tif","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/dem_big.tif"],
    #             "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/merge_greatlake_modified_FABDEM.tif")
    #
    # Find.clip("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/merge_greatlake_modified_FABDEM.tif",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/mask_big.shp",
    #           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/mask_dem_big.tif")
    # A = Raster.get_raster("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/mask_dem_big.tif")
    # print(A.shape)
    #
    # B = Raster.get_raster("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/sink_big.tif")
    # print(B.shape)

    # special_process.PFD8_exorheic(
    #         "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/mask_dem_big.tif",
    #         "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/sink_big.tif",
    #         "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/fdr_big.tif",
    #         "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/outdir_big_.tif",-73.1551590, 46.0351532,32)

    # special_process.read_and_save_tif("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/small/dem_5.tif",
    #                                   "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/small/dem_5.png")
    # special_process.read_and_save_tif(
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/dem_big.tif",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/dem_big.png")

    # special_process.PFD8_exorheic(
    #         "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/small/dem_5.tif",
    #         "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/small/sink_5.tif",
    #         "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/small/fdr_5.tif",
    #         "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/small/outdir_5_.tif",-78.0215916,48.6285362,2)

    # process_FAB_SRTM.modify_SRTM_notda("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatLakeSRTM/")
    # process_FAB_SRTM.merge_SRTM_FABDEM("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatLakeSRTM/", "/data/zhangbin/TBP_Stream_Data/FABDEM",
    #                   "/data/zhangbin/TBP_Stream_Data/FAB.tif", "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_srtm.tif")
    #
    # A = Raster.get_raster("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_srtm.tif")
    # proj,geo,nodata = Raster.get_proj_geo_nodata("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_srtm.tif")
    #
    #
    # row1,col1 = valid.lonlat_to_pixel(geo,-87.1,48.1)
    # row2, col2 = valid.lonlat_to_pixel(geo,-85.9,46.9 )
    #
    # A[row1:row2,col1:col2] = 179
    # Raster.save_raster("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_srtm1.tif",A,proj,geo,gdal.GDT_Float32,nodata)
    # process_FAB_SRTM.clip_FABDEM("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_srtm1.tif",
    #             "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/dem_big.tif",
    #             "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_FAB.tif")
    # process_FAB_SRTM.cancha_model_smooth_connect("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_FAB.tif",
    #                             "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_srtm1.tif",
    #                             "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_modified_FABDEM.tif")

    # special_process.read_and_save_tif("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_modified_FABDEM.tif","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_modified_FABDEM.png")
    # special_process.read_and_save_tif("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_FAB.tif","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_FAB.png")
    # special_process.read_and_save_tif("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_srtm1.tif","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_srtm1.png")
    # special_process.read_and_save_tif("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_srtm.tif","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_srtm.png")
    # Find.clip("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_modified_FABDEM.tif","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/test.shp","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/test_dem.tif")


    def AA(file,out_file):
        # file = r'F:\青藏高原水体数据集\3\NA\NA支流\all.tif'
        A = Raster.get_raster(file)
        proj,geo,nodata = Raster.get_proj_geo_nodata(file)

        row,col = A.shape

        B = np.zeros((row,col),dtype = np.int8)
        B[A != nodata] = 1
        # B[0,:] = nodata
        # B[row-1,:] = nodata
        # B[:,col-1] = nodata
        Raster.save_raster(out_file,B,proj,geo,gdal.GDT_Byte,0)

    # AA("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/small/fdr_5.tif","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/small/sink_5.tif")
    # AA("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/mask_dem_big.tif","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/sink_big.tif")

    # special_process.PFD8_exorheic(r'F:\青藏高原水体数据集\3\NA\NA支流\tri_dem.tif',
    #                               r'F:\青藏高原水体数据集\3\NA\NA支流\sink.tif',
    #                               r'F:\青藏高原水体数据集\3\NA\NA支流\tri_fdr.tif',
    #                               r'F:\青藏高原水体数据集\3\NA\NA支流\outdir.tif',-109.9865708,50.9051636,64)
    # merge_files(["/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final.tif/NorthAmerica_fdr.tif","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/outdir_NA.tif",
    #              "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/outdir_NA_tri.tif","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/NA支流/outdir.tif"],
    #             "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_fdr1.tif")

    # A = Raster.get_raster("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/outdir_big_.tif")
    # proj,geo,nodata = Raster.get_proj_geo_nodata("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/outdir_big_.tif")
    # A[A == 0] = nodata
    # Raster.save_raster("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/outdir_big_1.tif",A,proj,geo,gdal.GDT_Byte,nodata)
    # merge_files(["/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final1.tif/NorthAmerica_fdr.tif",
    #              "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/outdir_big_1.tif",
    #              "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/small/outdir_5_.tif"],"/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_fdr2.tif")

    # paras = ["/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Australia_ACC_final/region/Australia_fdr.tif"]
    # for file in os.listdir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/Australia_5020049720/endorheic"):
    #     out_dfie_ = os.path.join("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/Australia_5020049720/endorheic",file,file+'_outdir_1.tif')
    #     if os.path.exists(out_dfie_):
    #         paras.append(out_dfie_)
    #
    # merge_files(paras,"/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Australia_fdr1.tif")


def run_stream_valid():

    # stream_valid.valid_stream("/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/LongRiver/SWORD_LongRiver_nodes_greater_1000km2.shp",
    #                           "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/LongRiver/LongRiver_basin_streams_1000km2.shp",
    #                           "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/LongRiver/SWORD_valid_point_result.shp",
    #                           "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/LongRiver/SWORD_valid_point_result.csv")

    # stream_valid.valid_stream(
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/SWORD_Congo_nodes_greater_1000km2.shp",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_basin_streams_1000km2.shp",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_SWORD_valid_point_result.shp",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Congo_valid_point_result.csv")

    # stream_valid.valid_stream(
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Danube/Danube_basin_SWORD_streams_nodes.shp",
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Danube/Danube_basin_streams_1000km2.shp",
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Danube/Danube_SWORD_valid_point_result.shp",
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Danube/Danube_valid_point_result.csv")

    # stream_valid.valid_stream(
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Mississippi_basin_streams_nodes_1000km21.shp",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Mississippi_basin_streams_1000km2.shp",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Mississippi_SWORD_valid_point_result.shp",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/Mississippi_valid_point_result.csv")

    stream_valid.valid_stream(
        "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Amazon/Amazon_basin_streams_nodes_1000km21.shp",
        "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Amazon/Amazon_basin_streams_1000km2.shp",
        "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Amazon/Amazon_SWORD_valid_point_result.shp",
        "/datanode05/zhangbin/FAB_hydrography/Final_results/stream_valid_basin/Amazon/Amazon_valid_point_result.csv")



# ----------------------------- process_calculate_acc_stream ----------------------------------#


# -------------------------------- check endorheic point---------------------------------------------

def run_check_endorheic_point():


    # application_endorheic.po_get_endorheic_final_point(
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/Africa_fdr.tif",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/Endoheic_application")
    # application_endorheic.extract_endorheic_watershed(
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/Africa_fdr.tif",
    #       "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/Endoheic_application/endorheic_point.csv",
    #       "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/Endoheic_application/Endorheic_basins")
    # application_endorheic.sbatch_raster_endorheic_venu(
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/Endoheic_application/Endorheic_basins",
    #                                                    "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/Endoheic_application/Endorheic_basins_shp")
    # application_endorheic.statis_endorheic_information("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/Endoheic_application/endorheic_point.csv",
    #                                                    "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/Endoheic_application/Endorheic_basins_shp")


    # application_endorheic.po_get_endorheic_final_point(
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/Asia_fdr.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/Endoheic_application")
    # application_endorheic.extract_endorheic_watershed("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/Asia_fdr.tif",
    #                                                   "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/Endoheic_application/endorheic_point.csv",
    #                                                   "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/Endoheic_application/Endorheic_basins")
    # application_endorheic.sbatch_raster_endorheic_venu("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/Endoheic_application/Endorheic_basins",
    #                                                    "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/Endoheic_application/Endorheic_basins_shp")
    # application_endorheic.statis_endorheic_information("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/Endoheic_application/endorheic_point.csv",
    #                                                    "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/Endoheic_application/Endorheic_basins_shp")





    # application_endorheic.po_get_endorheic_final_point(
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/Europe_fdr.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/Endoheic_application")
    # application_endorheic.extract_endorheic_watershed("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/Europe_fdr.tif",
    #                                                   "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/Endoheic_application/endorheic_point.csv",
    #                                                   "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/Endoheic_application/Endorheic_basins")
    # application_endorheic.sbatch_raster_endorheic_venu("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/Endoheic_application/Endorheic_basins",
    #                                                    "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/Endoheic_application/Endorheic_basins_shp")
    # application_endorheic.statis_endorheic_information("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/Endoheic_application/endorheic_point.csv",
    #                                                    "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/Endoheic_application/Endorheic_basins_shp")



    # application_endorheic.po_get_endorheic_final_point(
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/Siberia_fdr.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/Endoheic_application")
    # application_endorheic.extract_endorheic_watershed("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/Siberia_fdr.tif",
    #                                                   "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/Endoheic_application/endorheic_point.csv",
    #                                                   "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/Endoheic_application/Endorheic_basins")
    # application_endorheic.sbatch_raster_endorheic_venu("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/Endoheic_application/Endorheic_basins",
    #                                                    "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/Endoheic_application/Endorheic_basins_shp")
    # application_endorheic.statis_endorheic_information("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/Endoheic_application/endorheic_point.csv",
    #                                                    "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/Endoheic_application/Endorheic_basins_shp")



    # application_endorheic.po_get_endorheic_final_point(
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/Australia_fdr.tif",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/Endoheic_application")
    # application_endorheic.extract_endorheic_watershed(
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/Australia_fdr.tif",
    #                                                   "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/Endoheic_application/endorheic_point.csv",
    #                                                   "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/Endoheic_application/Endorheic_basins")
    # application_endorheic.sbatch_raster_endorheic_venu("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/Endoheic_application/Endorheic_basins",
    #                                                    "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/Endoheic_application/Endorheic_basins_shp")
    # application_endorheic.statis_endorheic_information("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/Endoheic_application/endorheic_point.csv",
    #                                                    "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/Endoheic_application/Endorheic_basins_shp")
    #
    #
    # # application_endorheic.po_get_endorheic_final_point(
    # #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/NorthAmerica_fdr.tif",
    # #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/Endoheic_application")
    # application_endorheic.extract_endorheic_watershed(
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/NorthAmerica_fdr.tif",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/Endoheic_application/endorheic_point.csv",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/Endoheic_application/Endorheic_basins")
    # application_endorheic.sbatch_raster_endorheic_venu(
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/Endoheic_application/Endorheic_basins",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/Endoheic_application/Endorheic_basins_shp")
    # application_endorheic.statis_endorheic_information(
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/Endoheic_application/endorheic_point.csv",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/Endoheic_application/Endorheic_basins_shp")


    # application_endorheic.po_get_endorheic_final_point(
    #     "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/SouthAmerica_fdr.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/Endoheic_application")
    # application_endorheic.extract_endorheic_watershed("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/SouthAmerica_fdr.tif",
    #                                                   "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/Endoheic_application/endorheic_point.csv",
    #                                                   "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/Endoheic_application/Endorheic_basins")
    # application_endorheic.sbatch_raster_endorheic_venu("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/Endoheic_application/Endorheic_basins",
    #                                                    "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/Endoheic_application/Endorheic_basins_shp")
    # application_endorheic.statis_endorheic_information("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/Endoheic_application/endorheic_point.csv",
    #                                                    "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/Endoheic_application/Endorheic_basins_shp")
    pass


def run_endorheic_basin_postprocess():



    paras = []
    shp_venu = "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/Endoheic_application/Endorheic_basins_shp/"
    with open("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/Endoheic_application/Endorheic_basin_already_delete.csv",'r') as f:

        reader = csv.reader(f)
        n = 0
        for i in reader:
            if n==0:
                n+=1
                continue
            paras.append(i[0])
        f.close()

    for file in os.listdir(shp_venu):

        ids = file.split('.')[0].split('_')[-1]
        if ids not in paras:
            os.remove(os.path.join(shp_venu,file))

    application_endorheic.merge_endorheic_basin_to_gpkg(
        "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/Endoheic_application/Endorheic_basins_shp/")
    # application_endorheic.merge_endorheic_basin_to_gpkg("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/Endoheic_application/Endorheic_basins_shp")
    # application_endorheic.merge_endorheic_basin_to_gpkg("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/Endoheic_application/Endorheic_basins_shp/")
    # application_endorheic.merge_endorheic_basin_to_gpkg("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/Endoheic_application/Endorheic_basins_shp/")
    # application_endorheic.merge_endorheic_basin_to_gpkg("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/Endoheic_application/Endorheic_basins_shp/")
    # application_endorheic.merge_endorheic_basin_to_gpkg("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/Endoheic_application/Endorheic_basins_shp/")
    # application_endorheic.merge_endorheic_basin_to_gpkg("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/Endoheic_application/Endorheic_basins_shp/")

    pass
def get_stream_order():

    # sink.Cal_Stream("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/FINAL_ACC.tif",
    #                 "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/stream_5000.tif",2222222)
    # sink.stream_class("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/stream_5000.tif",
    #                   "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/Asia_fdr.tif")
    # #
    # sink.Cal_Stream("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/FINAL_ACC.tif",
    #                 "/datanode05/zhangbin/FAB_hydrography/Final_results/stream/SouthAmerica/stream_5000.tif", 2222222)
    # sink.stream_class("/datanode05/zhangbin/FAB_hydrography/Final_results/stream/SouthAmerica/stream_5000.tif",
    #                   "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/SouthAmerica_fdr.tif")
    # #
    # sink.Cal_Stream("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/FINAL_ACC.tif",
    #                 "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/stream_2000.tif", 2222222)
    # sink.stream_class("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/stream_2000.tif",
    #                   "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/Siberia_fdr.tif")
    # #
    sink.Cal_Stream("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/FINAL_ACC.tif",
                    "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/stream_5000.tif", 2222222)
    sink.stream_class("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/stream_5000.tif",
                       "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/Europe_fdr.tif")


    # sink.Cal_Stream("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/FINAL_ACC.tif",
    #                 "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/stream_2000.tif", 2222222)
    # sink.stream_class("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/stream_2000.tif",
    #                   "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/NorthAmerica_fdr.tif")
    #
    # sink.Cal_Stream("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/FINAL_ACC.tif",
    #                 "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Stream/Africa/stream_5000.tif",
    #                 2222222)
    # sink.stream_class("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Stream/Africa/stream_5000.tif",
    #                   "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/Africa_fdr.tif")
    # sink.Cal_Stream("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Australia_ACC_final1/FINAL_ACC.tif",
    #                 "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Stream/Australia/stream_5000.tif", 2222222)
    # sink.stream_class("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Stream/Australia/stream_5000.tif",
    #                   "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Australia_ACC_final1/Australia_fdr.tif")


    # sink.Cal_Stream("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_ACC_final/FINAL_ACC.tif",
    #                 "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_ACC_final/stream_5000.tif", 2222222)
    # sink.stream_class("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_ACC_final/stream_5000.tif",
    #                   "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_ACC_final/Arctic_fdr.tif")


    # sink.Cal_Stream("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic_north_island/Arctic_north_island_ACC_final/FINAL_ACC.tif",
    #                 "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic_north_island/Arctic_north_island_ACC_final/stream_2000.tif", 246913)
    # sink.stream_class("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic_north_island/Arctic_north_island_ACC_final/stream_2000.tif",
    #                   "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic_north_island/Arctic_north_island_ACC_final/Arctic_north_island_fdr.tif")



    # sink.Cal_Stream("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Greenland/Greenland_ACC_final/FINAL_ACC.tif",
    #                 "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Greenland/Greenland_ACC_final/stream_2000.tif",
    #                 246913)
    # sink.stream_class("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Greenland/Greenland_ACC_final/stream_2000.tif",
    #                   "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Greenland/Greenland_ACC_final/Greenland_fdr.tif")



    pass

# -------------------------------- check endorheic point---------------------------------------------





# --------------------------------------- valid of NHDPLUS -----------------------------------------

def run_valid_NHDPLUS():
    import compare_NHD

    gauge = "/datanode05/zhangbin/FAB_hydrography/initial_data/NHDPLUS/WBD_shp/NHDGauge.shp"
    gdb_path = "/datanode05/zhangbin/FAB_hydrography/initial_data/NHDPLUS/WBD_shp/Export_Output.shp"
    # huc12_layer = "WBDHU12"  # 流域边界图层
    # compare_NHD.construct_topoTree(gdb_path,huc12_layer,"/datanode05/xiefy/zhangb/FAB_Hydro/NHD_valid/tree.json","/datanode05/xiefy/zhangb/FAB_Hydro/NHD_valid/nodownids.json")
    # compare_NHD.construct_HU_area(gdb_path,huc12_layer,"/datanode05/xiefy/zhangb/FAB_Hydro/NHD_valid/HU_area.json")
    # compare_NHD.cal_csi(gauge, gdb_path, "/datanode05/xiefy/zhangb/FAB_Hydro/NHD_valid/tree.json",
    #         "/datanode05/xiefy/zhangb/FAB_Hydro/NHD_valid/HU_area.json",
    #         "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_fdir.tif",
    #         "/datanode05/xiefy/zhangb/FAB_Hydro/NHD_valid/FABUPS",
    #         "/datanode05/xiefy/zhangb/FAB_Hydro/NHD_valid/NHDPlusUPS",
    #         "/datanode05/xiefy/zhangb/FAB_Hydro/NHD_valid/NHD_area_trial.csv")

    # compare_NHD.modify_outlet("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/NorthAmerica_fdr.tif",
    #                           "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/FINAL_ACC.tif",
    #                           gauge,"/datanode05/xiefy/zhangb/FAB_Hydro/NHD_valid/NHD_area_trial.csv")
    #
    valid.one_time_extract_basin(
        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/NorthAmerica_fdr.tif",
        "/datanode05/xiefy/zhangb/FAB_Hydro/NHD_valid/NHD_gauge_area_modified_outlet.csv")

    # compare_NHD.sbatch_extract_basin()


def GRDC_valid():

    import valid

    # valid.get_basin_main('', "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/Europe_fdr.tif",
    #                      "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/FINAL_ACC.tif",
    #                      "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/GRDC_valid",'EU')
    valid.one_time_extract_basin(
        "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/Europe_fdr.tif",
                        "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/GRDC_valid/GRDC_station.csv")
    # valid.extract_basin("/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_ACC_final/Europe_fdr.tif",
    #                     "/datanode05/zhangbin/FAB_hydrography/GRDC_valid/Europe/Europe_GRDC_station.csv")


    # valid.get_basin_main('', "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_ACC_final/Arctic_fdr.tif",
    #                      "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_ACC_final/FINAL_ACC.tif",
    #                      "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_ACC_final/GRDC_valid",'AR')

    # valid.extract_basin("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_ACC_final-结果版/Arctic_fdr.tif",
    #                     "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_ACC_final-结果版/GRDC_valid/GRDC_station.csv")
    # valid.one_time_extract_basin("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_ACC_final/Arctic_fdr.tif",
    #                              "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_ACC_final/GRDC_valid/GRDC_station.csv")
    # valid.sbatch_cal_area("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_ACC_final-结果版/GRDC_valid")


    # valid.get_basin_main('',"/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/Asia_fdr.tif",
    #                      "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/FINAL_ACC.tif",
    #                      "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/GRDC_valid/","AS")
    # valid.extract_basin("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/Asia_fdr.tif",
    #                     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/GRDC_valid/GRDC_station.csv")
    # valid.one_time_extract_basin("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/Asia_fdr.tif",
    #                     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/GRDC_valid/GRDC_station.csv")
    # valid.sbatch_cal_area("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/GRDC_valid")
    # valid.get_CSI_lower_05("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/GRDC_valid/CSI.csv",
    #                        "/datanode05/zhangbin/FAB_hydrography/Result/Asia_ACC_final/FINAL_ACC.tif",
    #                        "/datanode05/zhangbin/FAB_hydrography/GRDC_valid/Asia/FAB_basins_vector",
    #                        "/datanode05/zhangbin/FAB_hydrography/GRDC_valid/Asia/FAB_basins_vector_lower05")

    # valid.get_basin_main('',"/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/NorthAmerica_fdr.tif",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/FINAL_ACC.tif",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/GRDC_valid",'NA')
    # valid.extract_basin("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/NorthAmerica_fdr.tif",
    #                     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/GRDC_valid/NorthAmerica/NorthAmerica_GRDC_station.csv")
    # valid.one_time_extract_basin("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/NorthAmerica_fdr.tif",
    #                              "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/GRDC_valid/GRDC_station.csv")
    # valid.sbatch_cal_area("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/GRDC_valid/NorthAmerica")
    # valid.get_CSI_lower_05("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/GRDC_valid/NorthAmerica/CSI.csv",
    #                        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final.tif/FINAL_ACC.tif",
    #                        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/GRDC_valid/NorthAmerica/FAB_basins_vector",
    # #                        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/GRDC_valid/NorthAmerica/FAB_basins_vector_lower05")

    # Merge.merge_update_fdr("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/region/",
    #                         "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/SouthAmerica_fdr.tif")
    # valid.get_basin_main('',"/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final1/region/SouthAmerica_fdr.tif",
    #                      "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final1/FINAL_ACC.tif",
    #                      "/datanode05/zhangbin/FAB_hydrography/GRDC_valid/SouthAmerica")
    # valid.extract_basin("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final1/region/SouthAmerica_fdr.tif",
    #                     "/datanode05/zhangbin/FAB_hydrography/GRDC_valid/SouthAmerica/SouthAmerica_GRDC_station.csv")
    # valid.sbatch_cal_area("/datanode05/zhangbin/FAB_hydrography/GRDC_valid/SouthAmerica")



    # Merge.merge_update_fdr("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Greenland_ACC_final/region",
    #                        "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Greenland_ACC_final/Greenland_fdr.tif")
    # Merge.mask_fdr_acc("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Greenland_ACC_final/Greenland_fdr.tif",
    #                    "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Greenland_ACC_final/FINAL_ACC.tif",
    #                    "")
    # valid.get_basin_main('', "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Greenland_ACC_final/Greenland_fdr.tif",
    #                      "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Greenland_ACC_final/FINAL_ACC.tif",
    #                      "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/GRDC_valid/Greenland")
    # valid.sbatch_cal_area("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/GRDC_valid/Greenland")

    # valid.get_basin_main('', "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/Siberia_fdr.tif",
    #                      "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/FINAL_ACC.tif",
    #                      "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/Siberia","SI")
    # valid.extract_basin("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/Siberia_fdr.tif",
    #                     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/GRDC_valid/GRDC_station.csv")
    # valid.one_time_extract_basin("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/Siberia_fdr.tif",
    #                     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/GRDC_valid/GRDC_station.csv")
    # valid.sbatch_cal_area("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_ACC_final/GRDC_valid")
    # valid.get_CSI_lower_05("/datanode05/zhangbin/FAB_hydrography/GRDC_valid/Siberia/CSI.csv",
    #                        "/datanode05/zhangbin/FAB_hydrography/Result/Siberia_ACC_final/FINAL_ACC.tif",
    #                        "/datanode05/zhangbin/FAB_hydrography/GRDC_valid/Siberia/FAB_basins_vector",
    #                        "/datanode05/zhangbin/FAB_hydrography/GRDC_valid/Siberia/FAB_basins_vector_lower05")
    # search_outlet_manually.sbatch_midify("/datanode05/zhangbin/FAB_hydrography/GRDC_valid/Siberia/FAB_basins_vector_lower05",
    #                                      "/datanode05/zhangbin/FAB_hydrography/Result/Siberia_ACC_final/Siberia_fdr.tif")
    # valid.sbatch_cal_area("/datanode05/zhangbin/FAB_hydrography/GRDC_valid/Siberia/FAB_Basin2_lower05")


    # valid.get_basin_main("","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/Africa_fdr.tif",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/FINAL_ACC.tif",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/GRDC_valid",'AF')
    # valid.extract_basin("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/Africa_fdr.tif",
    #                     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/GRDC_valid/GRDC_station.csv")
    # valid.one_time_extract_basin("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/Africa_fdr.tif",
    #                     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/GRDC_valid/GRDC_station.csv")
    # valid.sbatch_cal_area("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/GRDC_valid")
    # valid.get_CSI_lower_05("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/GRDC_valid/Africa/CSI.csv",
    #                        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/FINAL_ACC.tif",
    #                        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/GRDC_valid/Africa/FAB_basins_vector",
    #                        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/GRDC_valid/Africa/FAB_basins_vector_lower05")

    #
    # valid.get_basin_main("", "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/Australia_fdr.tif",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/FINAL_ACC.tif",
    #                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/GRDC_valid",'AU')
    # valid.extract_basin("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/Australia_fdr.tif",
    #                     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/GRDC_valid/GRDC_station.csv")
    # valid.one_time_extract_basin("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/Australia_fdr.tif",
    #                     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/GRDC_valid/GRDC_station.csv")
    # valid.sbatch_cal_area("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/GRDC_valid")
    # valid.get_CSI_lower_05("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/GRDC_valid/Australia/CSI.csv",
    #                        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Australia_ACC_final/FINAL_ACC.tif",
    #                        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/GRDC_valid/Australia/FAB_basins_vector",
    #                        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/GRDC_valid/Australia/FAB_basins_vector_lower05")

# --------------------------------------- valid of NHDPLUS -----------------------------------------


def copy_file():


    shutil.copy("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Australia_ACC_final/region/Australia_fdr.tif",
                "/datanode05/zhangbin/FAB_hydrography/Final_results/fdr/Australia_fdr.tif")
    shutil.copy("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Australia_ACC_final/FINAL_ACC.tif",
                "/datanode05/zhangbin/FAB_hydrography/Final_results/upa/Australia_upa.tif")

    shutil.copy("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final.tif/NorthAmerica_fdr.tif",
                "/datanode05/zhangbin/FAB_hydrography/Final_results/fdr/NorthAmerica_fdr.tif")
    shutil.copy("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final.tif/FINAL_ACC.tif",
                "/datanode05/zhangbin/FAB_hydrography/Final_results/upa/NorthAmerica_upa.tif")

    shutil.copy("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/Africa_fdr.tif",
                "/datanode05/zhangbin/FAB_hydrography/Final_results/fdr/Africa_fdr.tif")
    shutil.copy("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/FINAL_ACC.tif",
                "/datanode05/zhangbin/FAB_hydrography/Final_results/upa/Africa_upa.tif")

    shutil.copy("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/SouthAmerica_fdr.tif",
                "/datanode05/zhangbin/FAB_hydrography/Final_results/fdr/SouthAmerica_fdr.tif")
    shutil.copy("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/FINAL_ACC.tif",
                "/datanode05/zhangbin/FAB_hydrography/Final_results/upa/SouthAmerica_upa.tif")

    shutil.copy("/datanode05/zhangbin/FAB_hydrography/Result/Siberia_ACC_final/Siberia_fdr.tif",
                "/datanode05/zhangbin/FAB_hydrography/Final_results/fdr/Siberia_fdr.tif")
    shutil.copy("/datanode05/zhangbin/FAB_hydrography/Result/Siberia_ACC_final/FINAL_ACC.tif",
                "/datanode05/zhangbin/FAB_hydrography/Final_results/upa/Siberia_upa.tif")

    shutil.copy("/datanode05/zhangbin/FAB_hydrography/Result/Asia_ACC_final/Asia_fdr.tif",
                "/datanode05/zhangbin/FAB_hydrography/Final_results/fdr/Asia_fdr.tif")
    shutil.copy("/datanode05/zhangbin/FAB_hydrography/Result/Asia_ACC_final/FINAL_ACC.tif",
                "/datanode05/zhangbin/FAB_hydrography/Final_results/upa/Asia_upa.tif")

    shutil.copy("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Arctic_ACC_final/Arctic_fdr.tif",
                "/datanode05/zhangbin/FAB_hydrography/Final_results/fdr/Arctic_fdr.tif")
    shutil.copy("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Arctic_ACC_final/FINAL_ACC.tif",
                "/datanode05/zhangbin/FAB_hydrography/Final_results/upa/Arctic_upa.tif")

    shutil.copy("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Greenland_ACC_final/Greenland_fdr.tif",
                "/datanode05/zhangbin/FAB_hydrography/Final_results/fdr/Greenland_fdr.tif")
    shutil.copy("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Greenland_ACC_final/FINAL_ACC.tif",
                "/datanode05/zhangbin/FAB_hydrography/Final_results/upa/Greenland_upa.tif")




def check_huan_all():
    """
    检查流向中有没有环
    :return:
    """

    # check_counterflow.D("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_fdr2.tif")
    # check_counterflow.D("/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/SouthAmerica_fdr.tif")
    # check_counterflow.D("/datanode05/zhangbin/FAB_hydrography/Result/Asia_ACC_final1/Asia_fdr.tif")
    # check_counterflow.D("/datanode05/zhangbin/FAB_hydrography/Result/Siberia_ACC_final/Siberia_fdr.tif")
    # check_counterflow.D("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final2/NorthAmerica_fdr.tif")
    # check_counterflow.D("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/Africa_fdr.tif")
    # check_counterflow.D("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Australia_ACC_final1/Australia_fdr.tif")
    # check_counterflow.D("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Arctic_ACC_final/Arctic_fdr.tif")
    # check_counterflow.D("/datanode05/zhangbin/FAB_hydrography/run_data/Greenland/initial_dir.tif")




    # -------------------------------------------------------------------------------------------
    # postprocess.sbatch_process_huan("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/AR_huan.xlsx",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/_BURN_DEM.tif",
    #                                 "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_fdr1.tif",
    #                                 "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/huan")
    # postprocess.sbatch_process_huan("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/As_huan.xlsx",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/_BURN_DEM.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_fdr1.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/huan")
    # postprocess.sbatch_process_huan("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Si_huan.xlsx",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/_BURN_DEM.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_fdr1.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/huan")
    # postprocess.sbatch_process_huan("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Eu_huan.xlsx",
    #                                 "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/_BURN_DEM.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_fdr2.tif",
    #                                 "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/huan")
    # postprocess.sbatch_process_huan(
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/Af_huan.xlsx",
    #     "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/_BURN_DEM.tif",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/Africa_fdr.tif",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/huan")
    # postprocess.sbatch_process_huan(
    #     "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Greenland_ACC_final/Gr_huan.xlsx",
    #     "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/_BURN_DEM.tif",
    #     "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Greenland_ACC_final/Greenland_fdr.tif",
    #     "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Greenland_ACC_final/huan")
    # postprocess.sbatch_process_huan(
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NA_huan.xlsx",
    #     "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/_BURN_DEM.tif",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_fdr1.tif",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/huan")
    # postprocess.sbatch_process_huan(
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Australia_ACC_final1/AU_huan.xlsx",
    #     "/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/_BURN_DEM.tif",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Australia_ACC_final1/Australia_fdr.tif",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Australia_ACC_final1/huan")



    # ---------------------------------把sink转成shp，合并到gpkg---------------------------------------
    def convert_raster_to_gpkg(venu,out):
        import geopandas as gpd
        shp_files = []
        for file in os.listdir(venu):
            if file.split('.')[-1] == 'csv':
                continue
            sink_file = os.path.join(venu,file,'sink_'+str(file)+'_.tif')
            out_shp = os.path.join(venu,file,'sink_'+str(file)+'_.shp')
            if os.path.exists(os.path.join(venu,file,'modified_dir_'+str(file)+'_.tif')):
                continue

            valid.raster2Feature(sink_file,out_shp)
            shp_files.append(out_shp)


        # 读取并合并
        gdf_list = []
        for shp in shp_files:
            gdf = gpd.read_file(shp)
            gdf_list.append(gdf)

        merged_gdf = gpd.GeoDataFrame(
            pd.concat(gdf_list, ignore_index=True),
            crs="EPSG:4326"
        )

        # 输出到 GPKG
        out_gpkg = out
        layer_name = "merged_sink"

        merged_gdf.to_file(out_gpkg, layer=layer_name, driver="GPKG")

    # convert_raster_to_gpkg("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/huan/",
    #                        "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/huan/sink.gpkg")

    # convert_raster_to_gpkg("/datanode05/zhangbin/FAB_hydrography/Result/Siberia_ACC_final/huan",
    #                        "/datanode05/zhangbin/FAB_hydrography/Result/Siberia_ACC_final/huan/sink.gpkg")

    # convert_raster_to_gpkg("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/huan",
    #                        "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/huan/sink.gpkg")

    # convert_raster_to_gpkg("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/huan/",
    #                        "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/huan/sink.gpkg")
    # convert_raster_to_gpkg("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Australia_ACC_final1/huan/",
    #                        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Australia_ACC_final1/sink.gpkg")
    # convert_raster_to_gpkg("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/huan/",
    #                        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/huan/sink.gpkg")

    # convert_raster_to_gpkg("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/huan/",
    #                        "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/huan/sink.gpkg")

    # convert_raster_to_gpkg("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/huan/",
    #                        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/huan/sink.gpkg")


    def clip_dem(venu,DEM_file):
        DEM = Raster.get_raster(DEM_file)
        proj1,geo1,d_nodata = Raster.get_proj_geo_nodata(DEM_file)

        for file in os.listdir(venu):
            if len(file.split('.')) > 1:
                continue


            fdr_file = os.path.join(venu, file, 'dir_' + str(file) + '_.tif')
            AA = Raster.get_raster(fdr_file)
            proj,geo,f_nodata = Raster.get_proj_geo_nodata(fdr_file)
            AA[AA!= f_nodata] = 1
            AA[AA != 1] = 0
            sink_file = os.path.join(venu, file, 'sink_' + str(file) + '_.tif')
            Raster.save_raster(sink_file,AA,proj,geo,gdal.GDT_Byte,0)


            out_dem = os.path.join(venu, file, 'dem_' + str(file) + '_.tif')
            if os.path.exists(os.path.join(venu, file, 'modified_dir_' + str(file) + '_.tif')):
                continue

            proj,geo,nodata = Raster.get_proj_geo_nodata(sink_file)
            sink = Raster.get_raster(sink_file)
            row,col = sink.shape
            s_row,s_col = valid.lonlat_to_pixel(geo1,geo[0],geo[3])
            dem = DEM[s_row:s_row+row,s_col:s_col+col].copy()
            Raster.save_raster(out_dem,dem,proj,geo,gdal.GDT_Float32,d_nodata)

    # clip_dem("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正/",
    #          "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/_BURN_DEM.tif")

    # clip_dem("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/修正/","/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/_BURN_DEM.tif")

    # clip_dem("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/huan/",
    #          "/datanode05/zhangbin/FAB_hydrography/run_data/Asia/_BURN_DEM.tif")
    #
    # clip_dem("/datanode05/zhangbin/FAB_hydrography/Result/Siberia_ACC_final/huan",
    #          "/datanode05/zhangbin/FAB_hydrography/run_data/Siberia/_BURN_DEM.tif")

    # clip_dem("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Arctic_ACC_final/huan",
    #          "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Arctic/_BURN_DEM.tif")
    #
    # clip_dem("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/修正/",
    #          "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/_BURN_DEM.tif")
    # clip_dem("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/修正/","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Africa/_BURN_DEM.tif")

    # clip_dem("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正缺失的海","/datanode05/zhangbin/FAB_hydrography/run_data/Europe/_BURN_DEM.tif")

    # clip_dem("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正缺失的海","/datanode05/zhangbin/FAB_hydrography/run_data/Europe/_BURN_DEM.tif")
    #
    #
    #
    #
    # postprocess.sbatch_huan_based_on_dem("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/huan/")
    # postprocess.sbatch_huan_based_on_dem("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/huan/")
    # postprocess.sbatch_huan_based_on_dem("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/修正")
    # postprocess.sbatch_huan_based_on_dem("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Arctic_ACC_final/huan")
    # postprocess.sbatch_huan_based_on_dem("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/NorthAmerica_ACC_final2/huan")


    # postprocess.merge_counterflow_outdir("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/huan/",
    #                                      "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_fdr2.tif",
    #                                      "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_fdr3.tif")
    #
    # postprocess.merge_counterflow_outdir("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正/",
    #                                      "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_fdr3.tif",
    #                                      "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_fdr4.tif")


    # check_counterflow.D("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_fdr1.tif")

    # postprocess.merge_counterflow_outdir("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/huan/",
    #                                      "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_fdr1.tif",
    #                                      "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_fdr2.tif")
    # check_counterflow.D("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_fdr1.tif")

    # postprocess.merge_counterflow_outdir("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/修正/",
    #                                      "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_fdr1.tif",
    #                                      "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_fdr2.tif")
    # check_counterflow.D("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Siberia/Siberia_fdr1.tif")

    postprocess.merge_counterflow_outdir("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/修正1/",
                                         "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_ACC_final-结果版/Arctic_fdr.tif",
                                         "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_fdr4.tif")

    # check_counterflow.D("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_fdr1.tif")
    #
    # postprocess.merge_counterflow_outdir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/修正/",
    #                                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_fdr1.tif",
    #                                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_fdr2.tif")
    # check_counterflow.D("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_fdr1.tif")

    # postprocess.merge_counterflow_outdir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Australia_ACC_final1/huan/",
    #                                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Australia_ACC_final1/Australia_fdr.tif",
    #                                      "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_fdr.tif")

    postprocess.merge_counterflow_outdir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/修正/",
                                         "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_fdr.tif",
                                         "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_fdr1.tif")


    pass

def run_extract_region_stream():


    def get_fdr_acc(acc_file,fdr_file,mask,occ_file):

        basename = os.path.basename(mask)
        name = basename.split('.')[0]

        basevenu = os.path.dirname(mask)
        out_fdr = os.path.join(basevenu,name+'_fdr.tif')
        out_acc = os.path.join(basevenu,name+'_acc.tif')
        out_occ = os.path.join(basevenu, name + '_Occ.tif')
        out_occshp = os.path.join(basevenu, name + '_Occ.shp')
        out_stream = os.path.join(basevenu, name + '_stream.tif')
        out_outlet = os.path.join(basevenu, name + '_rivermouth.shp')
        Find.clip(acc_file,mask,out_acc)
        Find.clip(fdr_file,mask,out_fdr)
        Find.clip(occ_file,mask,out_occ)
        A = Raster.get_raster(out_occ)
        proj,geo,nodata = Raster.get_proj_geo_nodata(out_occ)
        A[A!=nodata] = 1
        Raster.save_raster(out_occ,A,proj,geo,gdal.GDT_Byte,nodata)

        # valid.raster2Feature(out_occ,out_occshp)


        sink.Cal_Stream(out_acc,
                        out_stream, 11111)
        sink.stream_class(out_stream,
                          out_fdr)


        application_endorheic.get_river_mouth_point(out_fdr,out_stream,out_outlet,'')
    #
    get_fdr_acc("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/FINAL_ACC.tif",
                "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/Asia_fdr.tif",
                "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/珠江/ZJ.shp",
                "/datanode05/zhangbin/FAB_hydrography/run_data/Asia/Asia_4020006940/Asia_4020006940_Occ.tif")
    # get_fdr_acc("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/FINAL_ACC.tif",
    #             "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/Asia_fdr.tif",
    #             "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia/Asia_ACC_final/青藏内流湖/TBP.shp",
    #             "/datanode05/zhangbin/FAB_hydrography/run_data/Asia/Asia_4020050470/Asia_4020050470_Occ.tif")
    #
    # get_fdr_acc("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/FINAL_ACC.tif",
    #             "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/Europe_fdr.tif",
    #             "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/北德维纳河/NDvina.shp",
    #             "/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_2020041390/Europe_2020041390_Occ.tif")

    # get_fdr_acc("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/FINAL_ACC.tif",
    #             "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/Australia_fdr.tif",
    #             "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/澳大利亚北海岸/Adf.shp",
    #             "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/Australia/Australia_5020049720/Australia_5020049720_Occ.tif")


def last():
    # #
    # def A(file,out):
    #
    #     arr = Raster.get_raster(file)
    #     proj,geo, nodata = Raster.get_proj_geo_nodata(file)
    #
    #     arr[arr != nodata] = 1
    #     arr[arr != 1] = 0
    #     Raster.save_raster(out,arr,proj,geo,gdal.GDT_Byte,0)
    # A("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正/4/dir_4_.tif","/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正/4/sink_4_.tif")
    # A("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正/3/dir_3_.tif","/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正/3/sink_3_.tif")
    # A("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正/2/dir_2_.tif",
    #   "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正/2/sink_2_.tif")
    # A("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正/1/dir_1_.tif",
    #   "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正/1/sink_1_.tif")


    # dem_file = "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正/4/dem_4_.tif"
    # sink_file = "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正/4/sink_4_.tif"
    # dem = Raster.get_raster(dem_file)
    # proj,geo,d_nodata = Raster.get_proj_geo_nodata(dem_file)
    # sink = Raster.get_raster(sink_file)
    # proj,geo,s_nodata = Raster.get_proj_geo_nodata(sink_file)
    # sink[sink!=1] = 0
    # maskDEM = dem.copy()
    # maskDEM[sink != 1] = d_nodata
    # lihai = (49.107176,39.53331 )
    # col = int((float(lihai[0]) - geo[0]) / geo[1])
    # row = int((geo[3] - float(lihai[1])) / abs(geo[5]))
    # min_cell = (row,col)
    # A = heap_PF_D8.optimized_flow_repair(maskDEM, sink, 0, (min_cell[0], min_cell[1]))  # 堆栈优化后的版本
    # A[sink != 1] = 255
    # Raster.save_raster("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正/4/modified_dir_4_.tif", A, proj, geo, gdal.GDT_Byte, 255)

    # from genral_functions import *
    # A = Raster.get_raster("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正缺失的海/3.tif")
    # proj,geo,nodata = Raster.get_proj_geo_nodata("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正缺失的海/3.tif")
    # row_num,col_num = A.shape
    # vis = np.zeros((row_num,col_num),dtype = np.int8)
    #
    # lihai = (38.707981,58.098372 )
    # col = int((float(lihai[0]) - geo[0]) / geo[1])
    # row = int((geo[3] - float(lihai[1])) / abs(geo[5]))
    # # popcells = [(row,col)]
    # popcells = deque()
    # popcells.append((row,col))
    #
    # while popcells:
    #
    #     popcell = popcells.popleft()
    #     # print(popcell)
    #
    #     for k in range(8):
    #         next_cell = (popcell[0]+dmove[k][0],popcell[1]+dmove[k][1])
    #         if vis[next_cell[0], next_cell[1]] == 1:
    #             continue
    #         if A[next_cell[0],next_cell[1]] == nodata:
    #
    #             if vis[next_cell[0],next_cell[1]] == 1:
    #                 continue
    #             A[next_cell[0], next_cell[1]] = 2 ** ((k + 4) % 8)
    #             vis[next_cell[0],next_cell[1]] = 1
    #             popcells.append(next_cell)
    #
    # Raster.save_raster("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正缺失的海/3/modified_dir_3.tif",A,proj,geo,gdal.GDT_Byte,nodata)


    # data = Raster.get_raster("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Arctic_ACC_final/FINAL_ACC.tif")
    # data = data[data>0]
    # plt.boxplot(data, showfliers=False)
    # plt.savefig("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/final_fdr/Arctic_ACC_final/FINAL_ACC.svg")


    # out = {}
    # with open("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/huan/lookup_80_.csv",'r') as f:
    #     reader = csv.reader(f)
    #     for i in reader:
    #         out.setdefault(float(i[0]),(int(float(i[6])),int(float(i[7]))))
    #     f.close()
    #
    # for file in os.listdir("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/huan/"):
    #     if len(file.split('.')) > 1:
    #         continue
    #     if os.path.exists(os.path.join("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/huan/",file,"modified_dir_"+file+'_.tif')):
    #         continue
    #     fdr_file = os.path.join("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/huan/",file,"dir_"+file+'_.tif')
    #
    #     A = Raster.get_raster(fdr_file)
    #     proj,geo,nodata = Raster.get_proj_geo_nodata(fdr_file)
    #     row,col = A.shape
    #     now_cell = out[float(file)]
    #     A[now_cell[0],now_cell[1]] = 16
    #     for i in range(10000):
    #         now_cell = [now_cell[0],now_cell[1]-1]
    #         if not genral_functions.check_boundary(now_cell[0],now_cell[1],row,col):
    #             break
    #         if A[now_cell[0],now_cell[1]] == nodata:
    #             break
    #         A[now_cell[0], now_cell[1]] = 16
    #
    #     Raster.save_raster(os.path.join("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Africa_ACC_final/huan/",file,"modified_dir_"+file+'_.tif'),A,proj,geo,gdal.GDT_Byte,nodata)
    pass

def delineated_all_basin():
    """
    划分所有的流域，除了河流流域，还有内流终点，还有入海口
    :return:
    """

    # delineate_basins.delineate_all_basins("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/NorthAmerica_fdr.tif",
    #                                       "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final/FINAL_ACC.tif",
    #                                       "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/Basins/NA/stream_10w.tif",
    #                                       "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/Basins/NA/stream_10w_link.tif",
    #                                       "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/Basins/NA/basin_10w.tif",
    #                                       "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/Basins/NA/basin_10w.shp",
    #                                       90000000)

    # delineate_basins.delineate_all_basins(
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/Africa_fdr.tif",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final/FINAL_ACC.tif",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/Basins/AF/stream_10w.tif",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/Basins/AF/stream_10w_link.tif",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/Basins/AF/basin_10w.tif",
    #     "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/Basins/AF/basin_10w.shp",
    #     90000000)


    # delineate_basins.delineate_all_basins(
    #     "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/SouthAmerica_fdr.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final/FINAL_ACC.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Basins/SA/stream_10w.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Basins/SA/stream_10w_link.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Basins/SA/basin_10w.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Basins/SA/basin_10w.shp",
    #     90000000)

    # delineate_basins.delineate_all_basins(
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia_Europe_Siberia/ACC_final/Asia_Europe_Siberia_fdr.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia_Europe_Siberia/ACC_final/FINAL_ACC.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Basins/AS_EU_SI/stream_10w.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Basins/AS_EU_SI/stream_10w_link.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Basins/AS_EU_SI/basin_10w.tif",
    #     "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Basins/AS_EU_SI/basin_10w.shp",
    #     90000000)

    delineate_basins.delineate_all_basins(
        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/Australia_fdr.tif",
        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final/FINAL_ACC.tif",
        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/Basins/AU/stream_10w.tif",
        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/Basins/AU/stream_10w_link.tif",
        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/Basins/AU/basin_10w.tif",
        "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/Basins/AU/basin_10w.shp",
        90000000)

    # delineate_basins.delineate_all_basins(
    #     "/datanode05/zhangbin/Analysis_3North/DATA/FABDEM/final/fdr.tif",
    #     "/datanode05/zhangbin/Analysis_3North/DATA/FABDEM/final/acc.tif",
    #     "/datanode05/zhangbin/Analysis_3North/DATA/FABDEM/final/stream_10w.tif",
    #     "/datanode05/zhangbin/Analysis_3North/DATA/FABDEM/final/stream_10w_link.tif",
    #     "/datanode05/zhangbin/Analysis_3North/DATA/FABDEM/final/basin_10w.tif",
    #     "/datanode05/zhangbin/Analysis_3North/DATA/FABDEM/final/basin_10w.shp",
    #     8888888)



def run_zip():
    """
    压缩fdr和acc，用于发布
    :return:
    """

    # sink.zip("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/NA/NorthAmerica_ACC_final")
    # sink.zip("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AU/Australia_ACC_final")
    # sink.zip("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Fdr/AF/Africa_ACC_final")
    sink.zip("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_ACC_final/")
    sink.zip("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic_north_island/Arctic_north_island_ACC_final/")
    # sink.zip("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/SouthAmerica/")
    # sink.zip("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Asia_Europe_Siberia/ACC_final/")





if __name__ == "__main__":

    # GPU代码
    # merge_GSWO("/data/zhangbin/TBP_Stream_Data/GSWO/","/data/zhangbin/TBP_Stream_Data/origin_raster/GSWO_global.tif")
    # merge_OSM("/data/zhangbin/TBP_Stream_Data/OSM/OSM/OSM_WaterLayer_tif","/data/zhangbin/TBP_Stream_Data/origin_raster/OSM_global.tif")
    # read_and_save_tif("/data/zhangbin/TBP_Stream_Data/origin_raster/GSWO_global.tif","/data/zhangbin/TBP_Stream_Data/origin_raster/GSWO_global.png")
    # read_and_save_tif("/data/zhangbin/TBP_Stream_Data/origin_raster/OSM_global.tif","/data/zhangbin/TBP_Stream_Data/origin_raster/OSM_global.png")
    # copy_DEM_flat("/data/zhangbin/TBP_Stream_Data/Global_FABDEM/","/data/zhangbin/TBP_Stream_Data/Global_FABDEM_Flatten")
    # mergeDEM("/data/zhangbin/TBP_Stream_Data/Global_FABDEM_Flatten","/data/zhangbin/TBP_Stream_Data/origin_raster/FABDEM_global.tif")
    # read_and_save_tif("/data/zhangbin/TBP_Stream_Data/origin_raster/FABDEM_global.tif","/data/zhangbin/TBP_Stream_Data/origin_raster/FABDEM_global.png")
    # preprocess()
    # clip_code()
    # prepare_occurance()
    # run_burning()
    # run_clip_next_level()
    # run_fdir_drpression()
    # run_calculate_sink_infos()
    # run_process_endorheic()
    # run_process_exorheic()
    # run_clip_outDir()
    # run_merge()
    # run_check_endorheic_point()



    # run_produce_counterflower()
    # check_huan_all()



    # 特殊处理
    # run_special_endorheic_fdir()
    # run_special_exorheic_fdir()

    # check_sink("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica")
    # check_sink("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica")
    # po_cal_acc()



    # --------------------------------- valid ----------------------------------------
    # run_valid_NHDPLUS()
    # GRDC_valid()
    # copy_file()

    # run_extract_stream_valid_basin()

    # merge_files(["/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/补/result_fdr.tif","/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/Europe_fdr.tif"],
    #             "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_fdr4.tif")
    # run_acc_calculation()
    # # #
    # run_extract_stream()
    # get_stream_order()

    # run_check_endorheic_point()

    # application_endorheic.extract_endorheic_watershed("/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_fdr1.tif","","/datanode05/zhangbin/FAB_hydrography/run_data/Europe/Europe_ACC_final/Endoheic_application")

    # run_endorheic_basin_postprocess()
    # run_stream_valid()
    # run_extract_region_stream()

    # delineated_all_basin()
    # run_zip()

    pass
