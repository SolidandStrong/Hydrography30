# -*- coding: utf-8 -*-
"""
@Time ： 2024/9/20 21:12
@Auth ：
@File ：calculate_sink_area_volume.py
@IDE ：PyCharm
"""
from osgeo import gdal
import numpy as np
import os
import time
import math
import Raster
from genral_functions import *
from db import *
def calculate_sink_volume(dem_file,dir_file,sink_file):

    dem = Raster.get_raster(dem_file)
    dir = Raster.get_raster(dir_file)
    sink = Raster.get_raster(sink_file)

    proj,geo,sinkNodata = Raster.get_proj_geo_nodata(sink_file)
    proj, geo, demNodata = Raster.get_proj_geo_nodata(dem_file)
    proj, geo, dirNodata = Raster.get_proj_geo_nodata(dir_file)
    row,col = dem.shape
    origin_sink = sink.copy()
    # (i,j,h)
    sinkInfo = {}
    sinkExtent = {}
    delte_ids = set() # 记录与其他块相交的sink_id，后续会根据这个id进行删除

    for i in range(row):
        for j in range(col):
            sinkId = sink[i,j]
            if sinkId == sinkNodata:
                continue
            sinkCellType = check_sink_boundary((i,j),sink,sinkNodata)  # Record whether a sink cell is at its boundary : 0 is No while 1 is Yes
            sinkInfo.setdefault(sinkId,[]).append((i,j,dem[i,j],sinkCellType))

            # 计算四至
            sinkExtent.setdefault(sinkId,[math.inf,-1,math.inf,-1])  # minX maxX minY maxY
            if sinkExtent[sinkId][0] > i:
                sinkExtent[sinkId][0] = i
            if sinkExtent[sinkId][1] < i:
                sinkExtent[sinkId][1] = i
            if sinkExtent[sinkId][2] > j:
                sinkExtent[sinkId][2] = j
            if sinkExtent[sinkId][3] < j:
                sinkExtent[sinkId][3] = j

            if i==0 or i==row-1 or j==0 or j==col-1:
                delte_ids.add(sinkId)
    # 删除边界污染的sink
    for delte_id in delte_ids:
        sinkExtent.pop(delte_id)
        sinkInfo.pop(delte_id)




    print('Sinks are searched all successfully, sum {:d} sinks.'.format(len(sinkInfo)))


    # Start calculate area and volume
    # 1. Search the lowest cell
    # 2. Search the outlet : the lowest cell at the boundaries
    # 3. Calculate volume
    # 4. Record the information

    areaVolumeInfo = []
    EXTENT = np.zeros((row, col))
    EXTENT[:, :] = -9999

    outls = np.zeros((row, col))
    outls[:, :] = -9999

    for sinkId in sinkInfo:
        typeClass = 0
        sinkCells = sinkInfo[sinkId]

        # 1. / 2.
        lowestCell = [-1,-1,math.inf]
        outletCell = [-1,-1,math.inf]
        for cell in sinkCells:
            tempH = dem[cell[0],cell[1]]
            if tempH < lowestCell[2]:
                lowestCell = [cell[0],cell[1],tempH]

        EXTENT, temp_extent_cells = Get_extent_out(dem, sink, demNodata, sinkCells, row, col,
                                                   EXTENT)  # 搜索sink的外围边界，最后输出会有重叠，所以对每个sink单独处理
        outletCell = getOutlet(temp_extent_cells,sink,dem,dir,demNodata,sinkId)
        outls[outletCell[0],outletCell[1]] = 2
        volume = 0
        for volumeCell in sinkCells:
            if dem[volumeCell[0],volumeCell[1]] <= dem[outletCell[0],outletCell[1]]:
                volume += cellSize*cellSize*(dem[outletCell[0],outletCell[1]]-dem[volumeCell[0],volumeCell[1]])/1000000000

        # # 3.
        # if outletCell[0] == -1:
        #     print('***********')
        #     # # F1 寻找外围边界最低点流出
        #     EXTENT, temp_extent_cells = Get_extent_out(dem, sink, demNodata, sinkCells, row, col,
        #                                                EXTENT)  # 搜索sink的外围边界，最后输出会有重叠，所以对每个sink单独处理
        #     outletCell, volume = Get_OUT(dir, dem, sinkId, origin_sink, sinkNodata, temp_extent_cells,
        #                                                sinkCells, demNodata,row, col)
        #     # print(temp_outlet)
        #
        #     # F2 寻找边界的最高点
        #     # outletCell,outletCellDir = find_sink_outlet(sink,dem,sinkCells,sinkId)
        #     # # print(outletCell,dir[outletCell[0],outletCell[1]],outletCellDir)
        #     # dir[outletCell[0],outletCell[1]] = outletCellDir

        area = len(sinkCells)*cellSize*cellSize/1000/1000
        # volume = 0
        # for cell in sinkCells:
        #     volume += (outletCell[2]-cell[2])*cellSize*cellSize/1000/1000/1000
        if volume == 0:
            typeClass = 1  # 平地
        elif area > 50 and volume > 0.5:
            typeClass = 2  # 内流
        else:
            typeClass = 3  # 外流
        # 4.
        areaVolumeInfo.append([int(sinkId),area,volume,outletCell[0],outletCell[1],lowestCell[0],lowestCell[1],typeClass,min(outletCell[0],sinkExtent[sinkId][0]),max(outletCell[0],sinkExtent[sinkId][1]),min(outletCell[1],sinkExtent[sinkId][2]),max(outletCell[1],sinkExtent[sinkId][3]),dem_file])
    print('Area and volume are calculated successfully.')


    # Raster.save_raster(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\OUT1.tif',outls,proj,geo,gdal.GDT_Float32,-9999)
    # Write it to db
    path = os.path.dirname(sink_file)
    dbPath = os.path.join(path,'sinkdb.db')
    # if os.path.exists(dbPath):
    #     os.remove(dbPath)
    sinkdelete = []
    for info in delte_ids:
        sinkdelete.append([info,dem_file])
    if len(areaVolumeInfo) != 0:
        insert_sink_table(dbPath,areaVolumeInfo)
    if len(sinkdelete) != 0:
        insert_sinkdelete_table(dbPath,sinkdelete)
    print('Write sink information to {:s} successfully.'.format(dbPath))

def calculate_sink_volume_byDir(Venu):
    basename = os.path.basename(Venu)
    dem_file = os.path.join(Venu,basename+"_burnDEM.tif")
    dir_file = os.path.join(Venu,basename+"_burndir.tif")
    sink_file = os.path.join(Venu,basename+"_burnsink.tif")
    # Occ_file = os.path.join(Venu,basename+"_Occ_.tif")


    dbPath = os.path.join(Venu, basename + "_burndb.db")
    if os.path.exists(dbPath):
        return

    if not os.path.exists(sink_file):
        return

    dem = Raster.get_raster(dem_file)
    dir = Raster.get_raster(dir_file)
    sink = Raster.get_raster(sink_file)
    # Occ = Raster.get_raster(Occ_file)

    proj,geo,sinkNodata = Raster.get_proj_geo_nodata(sink_file)
    proj, geo, demNodata = Raster.get_proj_geo_nodata(dem_file)
    proj, geo, dirNodata = Raster.get_proj_geo_nodata(dir_file)
    row,col = dem.shape

    outdirmak_file = os.path.join(Venu, basename + "_mergemask.tif")
    outdirmask = np.zeros((row,col))  # 合并流向时用到的掩膜
    outdirmask[dem!=demNodata] = 1

    # origin_sink = sink.copy()
    # (i,j,h)
    sinkInfo = {}
    sinkExtent = {}
    delte_ids = set() # 记录与其他块相交的sink_id，后续会根据这个id进行删除

    for i in range(row):
        for j in range(col):
            sinkId = sink[i,j]
            if sinkId == sinkNodata:
                continue
            sinkCellType = check_sink_boundary((i,j),sink,sinkNodata)  # Record whether a sink cell is at its boundary : 0 is No while 1 is Yes
            sinkInfo.setdefault(sinkId,[]).append((i,j,dem[i,j],sinkCellType))

            # 计算四至
            sinkExtent.setdefault(sinkId,[math.inf,-1,math.inf,-1])  # minX maxX minY maxY
            if sinkExtent[sinkId][0] > i:
                sinkExtent[sinkId][0] = i
            if sinkExtent[sinkId][1] < i:
                sinkExtent[sinkId][1] = i
            if sinkExtent[sinkId][2] > j:
                sinkExtent[sinkId][2] = j
            if sinkExtent[sinkId][3] < j:
                sinkExtent[sinkId][3] = j

            if i==0 or i==row-1 or j==0 or j==col-1:
                delte_ids.add(sinkId)
    # 将边界上的sink掩膜去掉，制作后续合并流向的掩膜
    for delte_id in delte_ids:
        sinkCells = sinkInfo[delte_id]
        for sinkCell in sinkCells:
            outdirmask[sinkCell[0],sinkCell[1]] = 0

    Raster.save_raster(outdirmak_file,outdirmask,proj,geo,gdal.GDT_Byte,0)
    # 删除边界污染的sink
    for delte_id in delte_ids:
        sinkExtent.pop(delte_id)
        sinkInfo.pop(delte_id)



    print('Sinks are searched all successfully, sum {:d} sinks.'.format(len(sinkInfo)))


    # Start calculate area and volume
    # 1. Search the lowest cell
    # 2. Search the outlet : the lowest cell at the boundaries
    # 3. Calculate volume
    # 4. Record the information

    areaVolumeInfo = []
    EXTENT = np.zeros((row, col),dtype=np.float32)
    EXTENT[:, :] = -9999

    outls = np.zeros((row, col),dtype=np.float32)
    outls[:, :] = -9999

    for sinkId in sinkInfo:
        typeClass = 0
        sinkCells = sinkInfo[sinkId]

        # 1. / 2.
        lowestCell = [-1,-1,math.inf]
        outletCell = [-1,-1,math.inf]
        for cell in sinkCells:
            tempH = dem[cell[0],cell[1]]
            if tempH < lowestCell[2]:
                lowestCell = [cell[0],cell[1],tempH]

        EXTENT, temp_extent_cells = Get_extent_out(dem, sink, demNodata, sinkCells, row, col,
                                                   EXTENT)  # 搜索sink的外围边界，最后输出会有重叠，所以对每个sink单独处理
        outletCell = getOutlet(temp_extent_cells,sink,dem,dir,demNodata,sinkId)
        outls[outletCell[0],outletCell[1]] = 2
        volume = 0
        for volumeCell in sinkCells:
            if dem[volumeCell[0],volumeCell[1]] <= dem[outletCell[0],outletCell[1]]:
                volume += cellSize*cellSize*(dem[outletCell[0],outletCell[1]]-dem[volumeCell[0],volumeCell[1]])/1000000000

        # # 3.
        # if outletCell[0] == -1:
        #     print('***********')
        #     # # F1 寻找外围边界最低点流出
        #     EXTENT, temp_extent_cells = Get_extent_out(dem, sink, demNodata, sinkCells, row, col,
        #                                                EXTENT)  # 搜索sink的外围边界，最后输出会有重叠，所以对每个sink单独处理
        #     outletCell, volume = Get_OUT(dir, dem, sinkId, origin_sink, sinkNodata, temp_extent_cells,
        #                                                sinkCells, demNodata,row, col)
        #     # print(temp_outlet)
        #
        #     # F2 寻找边界的最高点
        #     # outletCell,outletCellDir = find_sink_outlet(sink,dem,sinkCells,sinkId)
        #     # # print(outletCell,dir[outletCell[0],outletCell[1]],outletCellDir)
        #     # dir[outletCell[0],outletCell[1]] = outletCellDir
        # 加入遥感判断外流湖泊：构建出口点2km缓冲区，检查有无水体穿过：有则外流
        area = len(sinkCells)*cellSize*cellSize/1000/1000
        # volume = 0
        # for cell in sinkCells:
        #     volume += (outletCell[2]-cell[2])*cellSize*cellSize/1000/1000/1000
        dH = dem[outletCell[0],outletCell[1]] - dem[lowestCell[0],lowestCell[1]]
        if volume == 0:
            typeClass = 1  # 平地
        elif area > 50 and volume > 0.5 and dH > 10:
            typeClass = 2  # 内流
        else:
            typeClass = 3  # 外流
        # 4.
        areaVolumeInfo.append([int(sinkId),area,volume,outletCell[0],outletCell[1],lowestCell[0],lowestCell[1],typeClass,min(outletCell[0],sinkExtent[sinkId][0]),max(outletCell[0],sinkExtent[sinkId][1]),min(outletCell[1],sinkExtent[sinkId][2]),max(outletCell[1],sinkExtent[sinkId][3]),dem_file])
    print('Area and volume are calculated successfully.')


    # Raster.save_raster(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\OUT1.tif',outls,proj,geo,gdal.GDT_Float32,-9999)
    # Write it to db
    path = os.path.dirname(sink_file)
    # dbPath = os.path.join(path,'sinkdb.db')
    dbPath = os.path.join(Venu, basename + "_burndb.db")
    # if os.path.exists(dbPath):
    #     os.remove(dbPath)
    sinkdelete = []
    for info in delte_ids:
        sinkdelete.append([info,dem_file])
    if len(areaVolumeInfo) != 0:
        insert_sink_table(dbPath,areaVolumeInfo)
    if len(sinkdelete) != 0:
        insert_sinkdelete_table(dbPath,sinkdelete)
    os.chmod(dbPath,0o777)
    print('Write sink information to {:s} successfully.'.format(dbPath))


def sbatch_calculate_sink_volume_byDir(venu):

    print(venu)
    dirnames = os.listdir(venu)
    for dirname in dirnames:
        dirpath = os.path.join(venu,dirname)

        print(dirpath)
        if len(dirname.split('.')) > 1:
            continue
        if os.path.exists(os.path.join(dirpath,dirname+'_ModifiedDir.tif')):
            continue
        try:
            calculate_sink_volume_byDir(dirpath)
        except Exception as e:
            print(f"An error occurred: {e}")


if __name__ == '__main__':

    # dem_file = r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\DEM.tif'
    # dir_file = r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\Dir1.tif'
    # sink_file = r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\sink1.tif'

    dem_file = r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\globalRegion\C_trail\N27E085_FABDEM_V1-0.tif'
    dir_file = r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\globalRegion\C_trail\dir_N27E085_FABDEM_V1-0.tif'
    sink_file = r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\globalRegion\C_trail\sink_N27E085_FABDEM_V1-0.tif'
    #
    # t1 = time.time()
    # calculate_sink_volume(dem_file,dir_file,sink_file)
    # t2 = time.time()
    # print('Time consuming {:2} h'.format((t2-t1)/60/60))

    # dem_file = "/datanode05/zhangbin/TBP_Stream/TBP/TBP_FABDEM.tif"
    # dir_file = "/datanode05/zhangbin/TBP_Stream/TBP/TBP_FABDEM_Dir.tif"
    # sink_file = "/datanode05/zhangbin/TBP_Stream/DATA/20240818/TBP_sink.tif"

    t1 = time.time()
    basename = os.path.basename("/datanode05/zhangbin/TBP_Stream/TBP/TBP_FABDEM_Dir")
    dem_file = os.path.join("/datanode05/zhangbin/TBP_Stream/TBP/TBP_FABDEM_Dir", basename + "_burnDEM_.tif")
    print(dem_file)
    # calculate_sink_volume(dem_file, dir_file, sink_file)
    # split1(dem_file, dir_file, sink_file)
    # split2(dem_file, dir_file, sink_file)
    t2 = time.time()
    print('Time consuming {:2} h'.format((t2 - t1) / 60 / 60))

    pass