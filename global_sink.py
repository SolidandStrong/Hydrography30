# -*- coding: utf-8 -*-
"""
@Time ： 2024/10/25 10:19
@Auth ：
@File ：global_sink.py
@IDE ：PyCharm
"""
import math
import os.path
import shutil
from genral_functions import *
from osgeo import gdal,ogr
import Raster
import numpy as np
from collections import defaultdict
import time
from db import *
import whitebox

wbt = whitebox.WhiteboxTools()
wbt.verbose = False

def process_0(demfile,outfile):
    """
    优先搜索边界上的0值加入队列，使用图搜索将8邻域内的0值栅格赋值为nodata
    :param demfile:
    :return:
    """
    dem = Raster.get_raster(demfile)
    proj,geo,nodata = Raster.get_proj_geo_nodata(demfile)

    row,col = dem.shape
    cells = []
    vis = np.zeros((row, col))
    dmove = [(0,1),(1,0),(0,-1),(-1,0)]
    # search original iterate cells
    for i in [0,row-1]:
        temp = []
        for j in range(col-1):

            if dem[i,j] == 0:
                temp.append((i,j))
            else:
                # 中断
                if len(temp) == 0:
                    continue

                if temp[-1][1]-temp[0][1]+1 >= 30:
                    cells += temp.copy()
                temp = []
    if len(temp) != 0:
        if temp[-1][1] - temp[0][1] + 1 >= 30:
            cells += temp.copy()

    for j in [0,col-1]:
        temp = []
        for i in range(row-1):
            if dem[i, j] == 0:
                temp.append((i, j))
            else:
                # 中断
                if len(temp) == 0:
                    continue
                if temp[-1][0] - temp[0][0] + 1 >= 30:
                    cells += temp.copy()
                temp = []
    if len(temp) != 0:
        if temp[-1][0] - temp[0][0] + 1 >= 30:
            cells += temp.copy()
    for cell in cells:
        vis[cell[0],cell[1]] = 1
    # start iterate for searching 4-neibor equal cells
    while cells:
        popCell = cells.pop()
        dem[popCell[0],popCell[1]] = nodata
        for k in dmove:
            nextCell  = (popCell[0]+k[0],popCell[1]+k[1])
            if 0<=nextCell[0]< row and 0<=nextCell[1]<col:
                if vis[nextCell[0],nextCell[1]] == 1:
                    continue
                if dem[nextCell[0],nextCell[1]] != 0:
                    continue
                cells.append(nextCell)
                vis[nextCell[0],nextCell[1]] = 1

    Raster.save_raster(outfile,dem,proj,geo,gdal.GDT_Float32,nodata)


def cal_dir_sink(demfile,dirVenu,sinkVenu):
    dembasename = os.path.basename(demfile)

    outdirfile = os.path.join(dirVenu,'dir_'+dembasename)
    outsinkfile = os.path.join(sinkVenu,'sink_'+dembasename)

    wbt.d8_pointer(demfile,outdirfile,esri_pntr=True)
    wbt.sink(demfile,outsinkfile)
    return outdirfile,outsinkfile

def get_sink(sinkFile,db_name):

    minLon,minLat,maxLon,maxLat = get_sink_name(sinkFile)
    print(minLon,minLat,maxLon,maxLat)
    sinkInfos = {}
    # print(sinkFile)
    sink = Raster.get_raster(sinkFile)
    proj,geo,sinkNodata = Raster.get_proj_geo_nodata(sinkFile)
    rownum,colnum = sink.shape
    # print(geo)

    t1 = time.time()
    # 初始化默认字典，以列表存储每个键的索引



    # 构建字典
    # # 获取唯一值及其对应的索引
    # unique_values, indices = np.unique(sink, return_inverse=True)
    # print(unique_values)
    # rows, cols = np.indices(sink.shape)
    # extentSink = defaultdict(list)
    # for value, row, col in zip(indices, rows.ravel(), cols.ravel()):
    #     extentSink[unique_values[value]].append((row, col))
    # extentSink.pop(sinkNodata)
    # print(len(extentSink))
    t2 = time.time()

    extentSink = {}
    for i in range(rownum):
        for j in range(colnum):
            if sink[i,j] == sinkNodata:
                continue
            extentSink.setdefault(sink[i,j],[]).append((i,j))
    neighborSink_table = []
    Sink_table = []
    for sinkId in extentSink:
        tempSinks = extentSink[sinkId]
        left,down = math.inf,math.inf
        right, up= 0,0

        for cell in tempSinks:
            up = max(up,cell[0]) # 最小纬度，最大行
            left = min(left,cell[1]) # 最小经度
            down = min(down,cell[0]) # 最大纬度
            right = max(right,cell[1])


        tempGeo = [geo[0] + geo[1] * left, geo[1], geo[2], geo[3] + geo[5] * down, geo[4], geo[5]]
        extends = [1 if left==0 else 0, 1 if right==colnum-1 else 0, 1 if up==rownum-1 else 0, 1 if down==0 else 0]  # 左 右 下 上
        lrud = [geo[0] + geo[1] * left, geo[0] + geo[1] * right, geo[3] + geo[5] * up, geo[3] + geo[5] * down]   # 0:最小经度 1：最大经度 2：最小纬度 3：最大纬度
        sinkInfos.setdefault((lrud[2],lrud[0],lrud[3],lrud[1]),tempGeo+extends)
        # print(extends)

        if sum(extends) == 0:
            # 完整的sink
            Sink_table.append(lrud+tempGeo+[int(down),int(up),int(left),int(right)]+[sinkFile])
        else:
            # print(left, right, up, down)
            # print([minLon+ geo[1] * left,minLon+ geo[1] * right,maxLat+ geo[5] * up,maxLat+ geo[5] * down])
            neighborSink_table.append([round(minLon+ geo[1] * left,8),round(minLon+ geo[1] * right,8),round(maxLat+ geo[5] * up,8),round(maxLat+ geo[5] * down,8)]+extends)

        # --------------------------------------------------------------------------------
        # -                             NeighborSink Table                               -
        # --------------------------------------------------------------------------------
        # - 最小经度 | 最大经度 | 最小纬度 | 最大纬度 | Geo 6参数 | 左邻接 | 右邻接 | 下邻接 |上邻接 |
        # --------------------------------------------------------------------------------
        # -        |   -1    |        |    -1   | Geo 6参数 | 0/1 | 0/1    | 0/1   |0/1   |
        # --------------------------------------------------------------------------------

        # --------------------------------------------------------------------------------
        # -                             Sink Table                                       -
        # --------------------------------------------------------------------------------
        # - 最小经度 | 最大经度 | 最小纬度 | 最大纬度 | Geo 6参数 | 四至行列号  | 所在的栅格文件     |
        # --------------------------------------------------------------------------------
        # -        |         |        |         | Geo 6参数 |                             -
        # --------------------------------------------------------------------------------

    if not os.path.exists(db_name):
        create_global_sq(db_name)
    # 如果该块内有邻接邻接关系则不记录完整的sink信息，减少后续的更新量
    # if len(neighborSink_table) == 0:
    #     if len(Sink_table) > 0:
    #         insert_GlobalSink_table(db_name,Sink_table,'globalSink')
    # else:
    #     os.remove(sinkFile)
    if len(neighborSink_table) > 0:
        insert_GlobalSink_table(db_name,neighborSink_table,'neighborSink')

def cal_neighbor(db_name,demVenu):
    waitMerge = set()
    infos = query_neighborSink(db_name)
    for info in infos:
        # print(info)
        nowLon = math.floor(info[0])
        nowLat = math.floor(info[2])

        left = (nowLon - info[-4],nowLat)
        right = (nowLon + info[-3],nowLat)
        down = (nowLon,nowLat-info[-2])
        up = (nowLon,nowLat+info[-1])


        temp_outsource_rec = (min(nowLon,left[0]),min(nowLat,down[1]),
                              max(nowLon,right[0]),max(nowLat,up[1]))

        # print(temp_outsource_rec)
        waitMerge.add(temp_outsource_rec)
    print(waitMerge)
    print('____________')

    # 寻找范围最大的合并外包矩形
    waitMerge = list(waitMerge)
    waitMerge.sort(key = lambda x:(abs(x[0]-x[2])+1)*(abs(x[1]-x[3])+1))
    deleteInfos = []
    A = search_max_rec(waitMerge)
    print(A)
    A.sort(key=lambda x: (abs(x[0] - x[2]) + 1) * (abs(x[1] - x[3]) + 1), reverse=True)

    Vis = query_Vis(db_name)
    extend = []
    flag = False
    for _ in A.copy():
        if _ not in Vis:
            # print(_,Vis)
            extend = _
            regionNum = (abs(extend[0]-extend[2])+1)*(abs(extend[1]-extend[3])+1)
            lons = [i for i in range(min(extend[0], extend[2]), max(extend[0], extend[2]) + 1)]
            lats = [i for i in range(min(extend[1], extend[3]), max(extend[1], extend[3]) + 1)]
            input_rasters = []
            for lon in lons:
                for lat in lats:
                    tempFileName = get_dem_name(lon, lat)
                    print(tempFileName)
                    tempFilePath = os.path.join(demVenu, tempFileName)
                    if not os.path.exists(tempFilePath):
                        continue
                    input_rasters.append(tempFilePath)
            # print(input_rasters)
            if len(input_rasters) != regionNum:
                continue
            A.remove(_)
            insert_Vis_table(db_name,[_],'Vis')
            flag = True
            break


    # print(A)
    for info in infos:
        # print(info)
        nowLon = math.floor(info[0])
        nowLat = math.floor(info[2])

        left = (nowLon - info[-4], nowLat)
        right = (nowLon + info[-3], nowLat)
        up = (nowLon, nowLat + info[-2])
        down = (nowLon, nowLat - info[-1])

        temp_outsource_rec = (
        min(nowLon, left[0], right[0], up[0], down[0]), min(nowLat, left[1], right[1], up[1], down[1]),
        max(nowLon, left[0], right[0], up[0], down[0]), max(nowLat, left[1], right[1], up[1], down[1]))
        if temp_outsource_rec not in A:
            deleteInfos.append([str(info[0]),str(info[1]),str(info[2]),str(info[3])])

    delete_GlobalSink_table(db_name, deleteInfos, 'neighborSink')

    ############
    if not flag:
        return (), ()

    return input_rasters,extend


def mosaic1(raster_files,output_path,nodata_value=-9999):
    """
    批量镶嵌栅格文件，将多个栅格合并成一个。

    Parameters:
    input_dir (str): 输入栅格文件所在的文件夹路径。
    output_path (str): 输出镶嵌后栅格文件的保存路径。
    nodata_value (float, optional): 设置无数据值。如果未指定，保持默认值。
    """
    # 原来的代码
    # venu = os.path.dirname(output_path)
    # # 使用gdal.Warp镶嵌栅格
    # vrt_options = gdal.BuildVRTOptions(srcNodata=nodata_value, VRTNodata=nodata_value)
    # vrt_path = os.path.join(venu, 'mosaic.vrt')
    # gdal.BuildVRT(vrt_path, raster_files, options=vrt_options)
    #
    # # 将VRT文件保存为实际的GeoTIFF格式
    # # gdal.Translate(output_path, vrt_path)
    # gdal.Translate(
    #     output_path,
    #     vrt_path,
    #     format="GTiff",
    #     creationOptions=["COMPRESS=LZW"]  # 使用 LZW 压缩
    # )
    # print(f"镶嵌完成，结果保存在: {output_path}")

    venu = os.path.dirname(output_path)

    # print(result)
    # 构造 BuildVRTOptions（正确方式）
    vrt_options = gdal.BuildVRTOptions(
        separate=False,
        resolution='highest',
        VRTNodata=nodata_value  # 设定 NoData 值，可根据你的影像情况修改
    )

    # 构建 VRT
    vrt_path = os.path.join(venu, 'dem.vrt')
    vrt_ds = gdal.BuildVRT(vrt_path, raster_files, options=vrt_options)
    vrt_ds = None  # 保存并关闭

    # 输出 GeoTIFF + 压缩
    gdal.Translate(
        output_path,
        vrt_path,
        format="GTiff",
        creationOptions=[
            "COMPRESS=LZW",
            "TILED=YES",
            "BIGTIFF=IF_SAFER"
        ]
    )
    os.remove(vrt_path)


def Merge_cal_sir_sink(extend,input_rasters,mergeDEMPath):


    # 定义输出镶嵌文件路径
    output_raster = os.path.join(mergeDEMPath,get_merge_name(extend))
    if len(input_rasters) == 0:
        return False,output_raster
    # 调用 mosaic 工具
    mosaic1(input_rasters,output_raster)
    return True,output_raster

def record_sinkinfos(sinkVenu,db_name):

    sinkFiles = os.listdir(sinkVenu)
    for sinkName in sinkFiles:
        if sinkName.split('.')[-1] != 'tif':
            continue
        print(sinkName)
        sinkFile = os.path.join(sinkVenu,sinkName)
        minLon, minLat, maxLon, maxLat = get_sink_name(sinkFile)

        sinkInfos = {}
        sink = Raster.get_raster(sinkFile)
        proj, geo, sinkNodata = Raster.get_proj_geo_nodata(sinkFile)
        rownum, colnum = sink.shape
        # print(geo)

        t1 = time.time()
        # 初始化默认字典，以列表存储每个键的索引
        extentSink = defaultdict(list)

        # 获取唯一值及其对应的索引
        unique_values, indices = np.unique(sink, return_inverse=True)
        rows, cols = np.indices(sink.shape)

        # 构建字典
        for value, row, col in zip(indices, rows.ravel(), cols.ravel()):
            extentSink[unique_values[value]].append((row, col))
        t2 = time.time()
        extentSink.pop(sinkNodata)
        neighborSink_table = []
        Sink_table = []
        for sinkId in extentSink:
            tempSinks = extentSink[sinkId]
            left, down = math.inf, math.inf
            right, up = 0, 0
            for cell in tempSinks:
                up = max(up, cell[0])  # 最小纬度，最大行
                left = min(left, cell[1])  # 最小经度
                down = min(down, cell[0])  # 最大纬度
                right = max(right, cell[1])
            # print(left,right,up,down)
            tempGeo = [geo[0] + geo[1] * left, geo[1], geo[2], geo[3] + geo[5] * down, geo[4], geo[5]]
            extends = [1 if left == 0 else 0, 1 if right == colnum - 1 else 0, 1 if up == rownum - 1 else 0,
                       1 if down == 0 else 0]  # 左 右 下 上
            lrud = [geo[0] + geo[1] * left, geo[0] + geo[1] * right, geo[3] + geo[5] * up,
                    geo[3] + geo[5] * down]  # 0:最小经度 1：最大经度 2：最小纬度 3：最大纬度
            sinkInfos.setdefault((lrud[2], lrud[0], lrud[3], lrud[1]), tempGeo + extends)
            # print(extends)
            if sum(extends) == 0:
                # 完整的sink
                Sink_table.append(lrud + tempGeo + [int(down), int(up), int(left), int(right)] + [sinkFile])
            else:
                neighborSink_table.append(
                    [int(minLon) + geo[1] * left, int(minLon) + geo[1] * right, int(maxLat) + geo[5] * up,
                     int(maxLat) + geo[5] * down] + tempGeo + extends)

            # --------------------------------------------------------------------------------
            # -                             NeighborSink Table                               -
            # --------------------------------------------------------------------------------
            # - 最小经度 | 最大经度 | 最小纬度 | 最大纬度 | Geo 6参数 | 左邻接 | 右邻接 | 下邻接 |上邻接 |
            # --------------------------------------------------------------------------------
            # -        |   -1    |        |    -1   | Geo 6参数 | 0/1 | 0/1    | 0/1   |0/1   |
            # --------------------------------------------------------------------------------

            # --------------------------------------------------------------------------------
            # -                             Sink Table                                       -
            # --------------------------------------------------------------------------------
            # - 最小经度 | 最大经度 | 最小纬度 | 最大纬度 | Geo 6参数 | 四至行列号  | 所在的栅格文件     |
            # --------------------------------------------------------------------------------
            # -        |         |        |         | Geo 6参数 |                             -
            # --------------------------------------------------------------------------------

        if not os.path.exists(db_name):
            create_global_sq(db_name)

        if len(Sink_table) > 0:
            insert_GlobalSink_table(db_name, Sink_table, 'globalSink')



# 开花算法
if __name__ == '__main__':
    # get_sink(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\globalRegion\studyarea\sink_N25E068_FABDEM_V1-0.tif',r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\globalRegion\sink1.db')
    # inputRasters, extend = cal_neighbor(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\globalRegion\sink1.db', r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\globalRegion\studydem')
    # print(extend)
    # process_0(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\globalRegion\mregedem\N2526E065066_FABDEM_V1-0.tif',r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\globalRegion\4\1.tif')
    # record_sinkinfos(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\globalRegion\3',r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\globalRegion\sink.db')

    process_0(r'F:\FABDEM\Global_FABDEM\N20W090-N30W080_FABDEM_V1-0\N25W082_FABDEM_V1-0.tif',r'F:\FABDEM\temp\3.tif')
    pass