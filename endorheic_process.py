# -*- coding: utf-8 -*-
"""
@Time ： 2024/12/19 11:11
@Auth ：
@File ：endorheic_process.py
@IDE ：PyCharm
"""
import math
import os
import shutil

import Raster
import numpy as np
from scipy.ndimage import label
from osgeo import gdal
from multiprocessing import Pool
from genral_functions import *
from process_sink import *
dmove=[(0,1),(1,1),(1,0),(1,-1),(0,-1),(-1,-1),(-1,0),(-1,1)]
dmove_dic = {1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1), 16: (0, -1), 32: (-1, -1), 64: (-1, 0), 128: (-1, 1)}
cellSize = 30 # Size of cell



def delete(A):

    AA = Raster.get_raster(A)
    proj,geo,nodata = Raster.get_proj_geo_nodata(A)
    num, la = check_disjoint_regions(AA)
    la[la!=1] = 0


    Raster.save_raster(r'F:\青藏高原水体数据集\DATA\TBP_endorheic\109524\109524_mask1_.tif',la,proj,geo,gdal.GDT_Byte,0)

    print(num)


def check_disjoint_regions(mask,waterLevel=0):
    """
    输入文件，判断有多少连通区
    :param data:
    :return:
    """
    # 打开栅格文件
    # data = Raster.get_raster(raster_file)

    # 创建掩膜：筛选感兴趣的区域（例如值为1的像素）
    # mask = data <= waterLevel  # 替换为你的目标条件，例如：data == 1
    # print(mask)
    # 使用连通性分析
    structure = np.ones((3, 3))  # 8邻域连通性
    labeled_array, num_features = label(mask, structure=structure)
    if num_features >= 1:
        # print(f"存在 {num_features} 个不相邻的图形集合。")
        return num_features,labeled_array
    else:
        # print("图形集合是相邻的。")
        return 0,labeled_array

def lower_water(rasterfile):
    raster = Raster.get_raster(rasterfile)
    proj,geo,nodata = Raster.get_proj_geo_nodata(rasterfile)
    # raster = np.round(raster, decimals=1)  # 0.1m
    raster = np.array(raster,dtype = np.int32)
    uniques = np.unique(raster) # 仅保留唯一值
    values = np.sort(uniques[uniques!=nodata], axis=None)  # 快速排序
    print(values)
    row,col = raster.shape
    output = np.zeros((row,col))

        # print(num_features)
    popMasks = [raster]
    confluceTree = []
    while popMasks:
        popMask = popMasks.pop()
        # print(popMask)
        uniques = np.unique(popMask)
        values = np.sort(uniques[uniques != nodata], axis=None)  # 快速排序
        for i in range(len(values) - 1, -1, -1):
            # print(values[i])
            waterLevel = values[i]
            num_features, labeled_array = check_disjoint_regions(raster, waterLevel)
            print(num_features)
            for j in range(1, num_features + 1):
                temp_mask = labeled_array == j
                output[temp_mask] = j
                # print(temp_mask,output.max())
                area = np.sum(labeled_array == j)
                areaNum = int(50 * 1000000 / 30 / 30)   # 50km2
                # print(areaNum,area)
                if area < areaNum:
                    # 水满
                    raster[temp_mask] = waterLevel
                    print(waterLevel)
                else:
                    # 继续划分
                    popMasks.append(raster)
                    pass
            break

# F2:遥感识别法
# def initial_sink_upstream_extract(dirfile,outsinkfile,outbasinfile):
#
#     fdir = Raster.get_raster(dirfile)
#     proj,geo,f_nodata = Raster.get_proj_geo_nodata(dirfile)
#     mask = fdir==0
#
#     numbers, labeled = check_disjoint_regions(mask)
#
#     row,col = labeled.shape
#     labeled[labeled==0] = -9999
#     basins = labeled.copy()
#
#     for i in range(row):
#         for j in range(col):
#             if fdir[i,j] == f_nodata:
#                 continue
#             if basins[i,j] != -9999:
#                 continue
#
#             # 开始记录流路
#             popCells = [(i,j)]
#             flowPath = [(i,j)]
#             flag = -1
#             while popCells:
#                 popCell = popCells.pop()
#                 popdir = fdir[popCell[0],popCell[1]]
#                 if popdir == 0:
#                     # 遇到洼地，停止,回溯
#                     flag = basins[popCell[0],popCell[1]]
#                     break
#                 nextCell = (popCell[0]+dmove_dic[popdir][0],popCell[1]+dmove_dic[popdir][1])
#                 if basins[nextCell[0],nextCell[1]] != -9999:
#                     flag = basins[nextCell[0],nextCell[1]]
#                     break
#                 popCells.append(nextCell)
#                 flowPath.append(nextCell)
#             if flag == -1:
#                 continue
#
#             else:
#                 for cell in flowPath:
#                     basins[cell[0], cell[1]] = flag
#
#     Raster.save_raster(outsinkfile,labeled,proj,geo,gdal.GDT_Float32,-9999)
#     Raster.save_raster(outbasinfile,basins,proj,geo,gdal.GDT_Float32,-9999)

def initial_sink_upstream_extract(fdir,f_nodata):

    # fdir = Raster.get_raster(dirfile)
    # proj,geo,f_nodata = Raster.get_proj_geo_nodata(dirfile)
    mask = fdir==0

    numbers, labeled = check_disjoint_regions(mask)

    row,col = labeled.shape
    labeled[labeled==0] = -9999
    basins = labeled.copy()

    for i in range(row):
        for j in range(col):
            if fdir[i,j] == f_nodata:
                continue
            if basins[i,j] != -9999:
                continue

            # 开始记录流路
            popCells = [(i,j)]
            flowPath = [(i,j)]
            flag = -1
            while popCells:
                popCell = popCells.pop()
                popdir = fdir[popCell[0],popCell[1]]
                if popdir not in dmove_dic:
                    # 遇到洼地，停止,回溯
                    flag = basins[popCell[0],popCell[1]]
                    break
                nextCell = (popCell[0]+dmove_dic[popdir][0],popCell[1]+dmove_dic[popdir][1])
                # print(nextCell)
                if not check_boundary(nextCell[0],nextCell[1],row,col):
                    break
                if basins[nextCell[0],nextCell[1]] != -9999:
                    flag = basins[nextCell[0],nextCell[1]]
                    break
                popCells.append(nextCell)
                flowPath.append(nextCell)
            if flag == -1:
                continue
            else:
                for cell in flowPath:
                    basins[cell[0], cell[1]] = flag

    return basins
# F3 endorheic process
def process_endorheic(Venu):
    # mask_sin dir dem occ
    # 先把occ连通体寻找，保留occ average>50%的连通体
    # 以连通体内部最低点为出口点，试试PFD8策略

    baseName = os.path.basename(Venu)
    print(baseName)

    occfile = os.path.join(Venu, baseName + "_occ_.tif")
    fdirfile = os.path.join(Venu, baseName + "_dir_.tif")
    # sinkfile = os.path.join(Venu,baseName+"_occ_.tif")
    demfile = os.path.join(Venu, baseName + "_dem_.tif")
    outDirfile = os.path.join(Venu, baseName + "_outdir_1.tif")
    MASKFile = os.path.join(Venu, baseName + "_mask_.tif")
    MASK1 = Raster.get_raster(MASKFile)
    proj, geo, M_nodata = Raster.get_proj_geo_nodata(MASKFile)

    occ1 = Raster.get_raster(occfile)
    dem1 = Raster.get_raster(demfile)
    fdir1 = Raster.get_raster(fdirfile)
    proj, geo, o_nodata = Raster.get_proj_geo_nodata(occfile)
    proj, geo, f_nodata = Raster.get_proj_geo_nodata(fdirfile)
    _, _, dem_nodata = Raster.get_proj_geo_nodata(demfile)
    row, col = occ1.shape

    MASK = MASK1.copy()
    MASK[:, :] = M_nodata
    dem = dem1.copy()
    # dem[:, :] = dem_nodata
    occ = occ1.copy()
    occ[:, :] = o_nodata
    fdir = fdir1.copy()
    fdir[:, :] = f_nodata



    occmask1 = MASK1 > 0
    # occmask1 = occ1[MASK1 > 0].copy()
    # occmask1 = occ1.copy()

    # 仅保留average>50%的连通体 并且面积大于1km2
    print("开始识别连通体")
    numbers1, labledOcc1 = check_disjoint_regions(occmask1)
    # 检查四至：保证每个四至都在边界上的sink才是被处理的sink
    sink_extent = {}
    # 左列
    for i in range(row):
        sink_extent.setdefault(labledOcc1[i, 1], [0, 0, 0, 0])
        sink_extent[labledOcc1[i, 1]][2] = 1
    # 右边列
    for i in range(row):
        sink_extent.setdefault(labledOcc1[i, col - 2], [0, 0, 0, 0])
        sink_extent[labledOcc1[i, col - 2]][3] = 1
    # 上行
    for i in range(col):
        sink_extent.setdefault(labledOcc1[1, i], [0, 0, 0, 0])
        sink_extent[labledOcc1[1, i]][0] = 1
    # 下行
    for i in range(col):
        sink_extent.setdefault(labledOcc1[row - 2, i], [0, 0, 0, 0])
        sink_extent[labledOcc1[row - 2, i]][1] = 1

    print(sink_extent)
    for sink_extent_id in sink_extent:
        if sink_extent_id == 0:
            continue
        if sum(sink_extent[sink_extent_id]) == 4:
            final_mask = labledOcc1 == sink_extent_id

    MASK[final_mask] = MASK1[final_mask]
    occ[final_mask] = occ1[final_mask]
    dem[final_mask] = dem1[final_mask]
    fdir[final_mask] = fdir1[final_mask]

    # occmask = occ > 0
    occ1[occ1 != o_nodata] = 1
    occ1[occ1 != 1] = 0
    occmask = occ1.copy()
    numbers, labledOcc = check_disjoint_regions(occmask)

    print(numbers)

    savedOccCollectionIds = []
    for number in range(1, numbers + 1):
        labledmask = labledOcc == number
        if np.average(occ[labledmask]) < 50 or len(np.argwhere(labledOcc == number)) < 1111:
            labledOcc[labledmask] = 0
            continue
        savedOccCollectionIds.append(number)
    # 给这片洼地划分初始sink以及他们的上游
    print("连通体识别成功,有{:d}个连通体".format(len(savedOccCollectionIds)))

    print("开始划分初始洼地")
    initialbasins = initial_sink_upstream_extract(fdir,f_nodata)
    print("初始洼地划分数量{:d}".format(len(np.unique(initialbasins))))
    # Raster.save_raster(os.path.join(r'F:\青藏高原水体数据集\DATA\TBP_endorheic\173623', "sink_1111" +  ".tif"),initialbasins, proj, geo, gdal.GDT_Float32, -9999)

    # 处理连通体：
    tempIds = [-9999]
    print("开始处理连通体")
    outDEM = dem.copy()
    # outDEM[MASK==0] = 8888
    if len(savedOccCollectionIds) <= 1:
        # 识别出的水体数量小于等于1，则留向最低点
        # 掩膜sink

        maskSink = MASK.copy()
        maskSink[maskSink != 1] = 0  # 写成0 是因为好用：否则后面查最小值的时候会错误
        maskrow, maskcol = maskSink.shape

        maskDem = dem.copy()
        maskDem[maskSink == 0] = -9999
        maskDem[maskDem < -100] = 100000  # 先升高再赋空值，不知道为啥np.argmin识别不了-9999

        # 获取该区域的最小值的索引
        min_index_in_region = np.unravel_index(np.argmin(maskDem), maskDem.shape)
        outCell = (min_index_in_region[0], min_index_in_region[1])
        print("最低高程:{:f}".format(maskDem[outCell[0], outCell[1]]))

        # 记录出口点数据
        maskDir = fdir.copy()  # np.zeros((row,col))
        maskDir[maskSink == 1] = 0
        maskDir[outCell[0], outCell[1]] = 0
        maskDem[maskDem == 100000] = dem_nodata

        # 计算淹没内的流向
        # A = priorityFlood(maskSink, 0, maskDem, (outCell[0], outCell[1]), maskDir)  # Python
        A = queuePriorityFlood(maskSink, 0, maskDem, (outCell[0], outCell[1]), maskDir)  # Python

        # 写入流向
        for ti in range(maskrow):
            for tj in range(maskcol):
                if A[ti, tj] in [0, 255]:
                    # print(A[ti,tj])
                    continue
                fdir[ti, tj] = A[ti, tj]
                # print(A[ti,tj])
    else:
        for soid in savedOccCollectionIds:
            print(soid)
            bigCollection = np.zeros((row,col))   # 待处理的洼地掩膜
            soidmask = labledOcc==soid
            initialsinkIds = np.unique(initialbasins[soidmask])  # 记录occ连通体经过的洼地集水区，这些集水区最后合并成一个大的集合进行处理
            tempIds += list(initialsinkIds)
            # 定义需要重新赋值的固定值和新值的映射
            values_to_replace = initialsinkIds  # 需要替换的值
            new_value = 1  # 替换后的值
            # 使用 np.isin 找到需要替换的值并赋值
            bigCollection[np.isin(initialbasins,values_to_replace)] = new_value
            bigCollection[fdir==f_nodata] = 0
            outDEM[bigCollection==1] = dem_nodata
            # Raster.save_raster(os.path.join(r'F:\青藏高原水体数据集\DATA\TBP_endorheic\173623',"sink_"+str(soid)+".tif"),bigCollection,proj,geo,gdal.GDT_Byte,0)
            # break


            # 处理bigcollection
            cells = np.argwhere(bigCollection==1)
            minxy = np.amin(cells, axis=0)
            maxxy = np.amax(cells, axis=0)
            minX = minxy[0]
            maxX = maxxy[0]
            minY = minxy[1]
            maxY = maxxy[1]
            # print(minX,minY,maxX,maxY)

            # 掩膜sink
            maskSink = bigCollection[minX:maxX + 1, minY:maxY + 1].copy()
            maskSink[maskSink != 1] = 0  # 写成0 是因为好用：否则后面查最小值的时候会错误
            maskrow, maskcol = maskSink.shape

            maskDem = dem[minX:maxX + 1, minY:maxY + 1].copy()
            maskDem[maskSink == 0] = -9999
            maskDem[maskDem < -100] = 100000  # 先升高再赋空值，不知道为啥np.argmin识别不了-9999

            # 获取该区域的最小值的索引
            min_index_in_region = np.unravel_index(np.argmin(maskDem), maskDem.shape)
            outCell = (min_index_in_region[0],min_index_in_region[1])
            print("最低高程:{:f}".format(maskDem[outCell[0],outCell[1]]))

            # 记录出口点数据
            maskDir = fdir[minX:maxX + 1, minY:maxY + 1].copy()  # np.zeros((row,col))
            maskDir[maskSink == 1] = 0
            maskDir[outCell[0] , outCell[1]] = 0

            maskDem[maskDem == 100000] = dem_nodata
            # print(maskDir.shape,maskSink.shape,maskDem.shape)

            # 计算淹没内的流向
            # A = priorityFlood(maskSink, 0, maskDem, (outCell[0], outCell[1]), maskDir)  # Python
            A = queuePriorityFlood(maskSink, 0, maskDem, (outCell[0], outCell[1]), maskDir)  # Python

            # print(A.min(),A.max())
            # 写入流向
            for ti in range(maskrow):
                for tj in range(maskcol):
                    if A[ti, tj] in [0,255]:
                        # print(A[ti,tj])
                        continue
                    fdir[minX + ti, minY + tj] = A[ti, tj]
                    # print(A[ti,tj])

        print("连通体处理成功")

        temp_outDEMfile = os.path.join(Venu,baseName+"_tempoutdem_.tif")
        temp_outFillfile = os.path.join(Venu, baseName + "_tempoutFill_.tif")
        temp_outtdirfile = os.path.join(Venu, baseName + "_tempouttdir_.tif")

        Raster.save_raster(temp_outDEMfile, outDEM, proj, geo, gdal.GDT_Float32, dem_nodata)

        # 处理外流的洼地:搜索外边界上最低的点流到洼地以外的斜率最大的栅格
        # 正常填洼
        import whitebox
        wbt = whitebox.WhiteboxTools()
        # wbt.fill_depressions_planchon_and_darboux(temp_outDEMfile,temp_outFillfile,fix_flats=True)  # 会填充内部的无效值
        wbt.fill_depressions(temp_outDEMfile,temp_outFillfile,fix_flats=True)
        wbt.d8_pointer(temp_outFillfile,temp_outtdirfile,esri_pntr=True)
        # sink.Cal_filled(temp_outDEMfile,temp_outFillfile)
        # sink.Cal_dir(temp_outFillfile,temp_outtdirfile)
        newd8 = Raster.get_raster(temp_outtdirfile)
        FillDem = Raster.get_raster(temp_outFillfile)
        # newd8[MASK == 0] = 255
        cells0 = np.argwhere(newd8 == 0)
        modifiedcells = []
        for temp_cell in cells0:
            outdir = [0,math.inf]
            for k in range(8):
                nextCell = (temp_cell[0]+dmove[k][0],temp_cell[1]+dmove[k][1])
                if not check_boundary(nextCell[0],nextCell[1],row,col):
                    continue
                if dem[nextCell[0],nextCell[1]] == dem_nodata:
                    continue
                if outDEM[nextCell[0],nextCell[1]] != dem_nodata:
                    continue
                nextCellH = dem[nextCell[0],nextCell[1]]
                if nextCellH > outdir[1]:
                    continue
                outdir = [2**k,nextCellH]
            if outdir[0] == 0:
                modifiedcells.append([temp_cell[0],temp_cell[1]])
            else:
                newd8[temp_cell[0],temp_cell[1]] = outdir[0]
        fdir[outDEM != dem_nodata] = newd8[outDEM != dem_nodata]

        # # 计算洼地上游，从边界上最低的与其他已有流向的栅格流出
        for temp_cell in modifiedcells:
            popCells = [(temp_cell[0], temp_cell[1])]
            sinkCells = [(temp_cell[0], temp_cell[1],FillDem[temp_cell[0],temp_cell[1]])]
            minX = math.inf
            minY = math.inf
            maxX = -1
            maxY = -1
            while popCells:
                popCell = popCells.pop()
                minX = min(minX,popCell[0])
                minY = min(minY,popCell[1])
                maxX = max(maxX,popCell[0])
                maxY = max(maxY,popCell[1])

                for kk in range(8):
                    nextCell = (popCell[0] + dmove[kk][0], popCell[1] + dmove[kk][1])
                    if not check_boundary(nextCell[0], nextCell[1], row, col):
                        continue
                    if outDEM[nextCell[0], nextCell[1]] == dem_nodata:
                        continue
                    nextDir = fdir[nextCell[0],nextCell[1]]
                    if 2**((kk+4)%8) != nextDir:
                        continue
                    popCells.append(nextCell)
                    sinkCells.append((nextCell[0],nextCell[1],FillDem[nextCell[0],nextCell[1]]))
            # 构建掩膜计算矩阵，从边界上最低的与其他已有流向的栅格流出
            maxHCell = [-1, -1, -100]  # x ,y, fdir
            sinkCells.sort(key=lambda x:x[2])
            flag = False
            for cell in sinkCells:
                for kk in range(8):
                    nextCell = [cell[0] + dmove[kk][0], cell[1] + dmove[kk][1]]
                    if not check_boundary(nextCell[0], nextCell[1], row, col):
                        continue
                    if dem[nextCell[0], nextCell[1]] == dem_nodata:
                        continue
                    temp_nextCell = [nextCell[0],nextCell[1],FillDem[nextCell[0],nextCell[1]]]
                    if temp_nextCell in sinkCells:
                        continue

                    flag = True
                    maxHCell = [cell[0],cell[1],2**kk]
                    break
                if flag:
                    break

            # 掩膜sink
            maskSink = np.zeros((maxX-minX+1,maxY-minY+1))
            for sinkCell in sinkCells:
                maskSink[sinkCell[0]-minX,sinkCell[1]-minY] = 1

            maskSink[maskSink != 1] = 0  # 写成0 是因为好用：否则后面查最小值的时候会错误
            maskrow, maskcol = maskSink.shape
            maskDem = dem[minX:maxX + 1, minY:maxY + 1].copy()
            maskDem[maskSink == 0] = dem_nodata
            # 记录出口点数据
            maskDir = fdir[minX:maxX + 1, minY:maxY + 1].copy()  # np.zeros((row,col))
            maskDir[maskSink == 1] = 0
            outCell = [maxHCell[0]-minX,maxHCell[1]-minY]
            maskDir[outCell[0], outCell[1]] = maxHCell[2]

            # 计算淹没内的流向
            B = queuePriorityFlood(maskSink, 0, maskDem, (outCell[0], outCell[1]), maskDir)  # Python

            # 写入流向
            for ti in range(maskrow):
                for tj in range(maskcol):
                    if B[ti, tj] in [0, 255]:
                        # print(A[ti,tj])
                        continue
                    fdir[minX + ti, minY + tj] = B[ti, tj]
                    # print(A[ti,tj])


        # if os.path.exists(temp_outDEMfile):
        #     os.remove(temp_outDEMfile)
        # if os.path.exists(temp_outFillfile):
        #     os.remove(temp_outFillfile)
        # if os.path.exists(temp_outtdirfile):
        #     os.remove(temp_outtdirfile)
    fdir[MASK1 == 0] = f_nodata
    Raster.save_raster(outDirfile, fdir, proj, geo, gdal.GDT_Byte, f_nodata)
    return fdir


def sbatch_process_endorheic(Venu):

    venudirs = os.listdir(Venu)
    runVenus = []
    for venudir in venudirs:
        venu = os.path.join(Venu,venudir)  # 二级根目录：存放具体文件_run_datas
        runVenus.append(venu)

    po = Pool(15)
    for info in runVenus:
        po.apply_async(process_endorheic,(info,))

    po.close()
    po.join()

    pass

def make_endorheic_mask_fdir(Venu):
    """
    给这个文件夹下的数据生产掩膜_mask_.tif和初始流向_burndir，用于处理特殊的内流区流向
    :param Venu:
    :return:
    """
    baseName = os.path.basename(Venu)
    out_mask_path = os.path.join(Venu,baseName + '_mask_.tif')
    baseTifPath = os.path.join(Venu,baseName + "_burnDEM.tif")


    baseTif = Raster.get_raster(baseTifPath)
    proj,geo,nodata = Raster.get_proj_geo_nodata(baseTifPath)
    row,col = baseTif.shape

    MASK = np.zeros((row,col),dtype = np.int8)
    MASK[baseTif != nodata] = 1

    Raster.save_raster(out_mask_path,MASK,proj,geo,gdal.GDT_Byte,0)


    out_fdir_path = os.path.join(Venu, baseName + '_burndir.tif')
    sink.Cal_dir(baseTifPath,out_fdir_path)



def process_endorheic1(Venu):
    # mask_sin dir dem occ
    # 先把occ连通体寻找，保留occ average>50%的连通体
    # 以连通体内部最低点为出口点，试试PFD8策略

    baseName = os.path.basename(Venu)
    print(baseName)

    occfile = os.path.join(Venu, baseName + "_Occ.tif")
    fdirfile = os.path.join(Venu, baseName + "_burndir.tif")
    # sinkfile = os.path.join(Venu,baseName+"_occ_.tif")
    demfile = os.path.join(Venu, baseName + "_burnDEM.tif")
    outDirfile = os.path.join(Venu, baseName + "_ModifiedDir.tif")
    MASKFile = os.path.join(Venu, baseName + "_mask_.tif")
    MASK1 = Raster.get_raster(MASKFile)
    proj, geo, M_nodata = Raster.get_proj_geo_nodata(MASKFile)

    occ1 = Raster.get_raster(occfile)
    dem1 = Raster.get_raster(demfile)
    fdir1 = Raster.get_raster(fdirfile)
    proj, geo, o_nodata = Raster.get_proj_geo_nodata(occfile)
    proj, geo, f_nodata = Raster.get_proj_geo_nodata(fdirfile)
    _, _, dem_nodata = Raster.get_proj_geo_nodata(demfile)
    row, col = occ1.shape

    MASK = MASK1.copy()
    MASK[:, :] = M_nodata
    dem = dem1.copy()
    dem[:, :] = dem_nodata
    occ = occ1.copy()
    occ[:, :] = o_nodata
    fdir = fdir1.copy()
    fdir[:, :] = f_nodata

    occmask1 = MASK1 > 0

    # 仅保留average>50%的连通体 并且面积大于10km2
    print("开始识别连通体")
    numbers1, labledOcc1 = check_disjoint_regions(occmask1)
    # 检查四至：保证每个四至都在边界上的sink才是被处理的sink
    sink_extent = {}
    # 左列
    for i in range(row):
        sink_extent.setdefault(labledOcc1[i, 1], [0, 0, 0, 0])
        sink_extent[labledOcc1[i, 1]][2] = 1
    # 右边列
    for i in range(row):
        sink_extent.setdefault(labledOcc1[i, col - 2], [0, 0, 0, 0])
        sink_extent[labledOcc1[i, col - 2]][3] = 1
    # 上行
    for i in range(col):
        sink_extent.setdefault(labledOcc1[1, i], [0, 0, 0, 0])
        sink_extent[labledOcc1[1, i]][0] = 1
    # 下行
    for i in range(col):
        sink_extent.setdefault(labledOcc1[row - 2, i], [0, 0, 0, 0])
        sink_extent[labledOcc1[row - 2, i]][1] = 1

    print(sink_extent)
    for sink_extent_id in sink_extent:
        if sink_extent_id == 0:
            continue
        if sum(sink_extent[sink_extent_id]) == 4:
            final_mask = labledOcc1 == sink_extent_id

    MASK[final_mask] = MASK1[final_mask]
    occ[final_mask] = occ1[final_mask]
    dem[final_mask] = dem1[final_mask]
    fdir[final_mask] = fdir1[final_mask]

    # occmask = occ > 0
    occ1[occ1 != o_nodata] = 1
    occ1[occ1 != 1] = 0
    occmask = occ1.copy()
    numbers, labledOcc = check_disjoint_regions(occmask)

    print(numbers)

    savedOccCollectionIds = []
    for number in range(1, numbers + 1):
        labledmask = labledOcc == number
        if np.average(occ[labledmask]) < 50 or len(np.argwhere(labledOcc == number)) < 1111:
            labledOcc[labledmask] = 0
            continue
        savedOccCollectionIds.append(number)
    # 给这片洼地划分初始sink以及他们的上游
    print("连通体识别成功,有{:d}个连通体".format(len(savedOccCollectionIds)))

    print("开始划分初始洼地")
    initialbasins = initial_sink_upstream_extract(fdir, f_nodata)
    print("初始洼地划分数量{:d}".format(len(np.unique(initialbasins))))
    # Raster.save_raster(os.path.join(r'F:\青藏高原水体数据集\DATA\TBP_endorheic\173623', "sink_1111" +  ".tif"),initialbasins, proj, geo, gdal.GDT_Float32, -9999)

    # 处理连通体：
    tempIds = [-9999]
    print("开始处理连通体")
    outDEM = dem.copy()
    # outDEM[MASK==0] = 8888
    if len(savedOccCollectionIds) <= 1:
        # 识别出的水体数量小于等于1，则留向最低点
        # 掩膜sink

        maskSink = MASK.copy()
        maskSink[maskSink != 1] = 0  # 写成0 是因为好用：否则后面查最小值的时候会错误
        maskrow, maskcol = maskSink.shape

        maskDem = dem.copy()
        maskDem[maskSink == 0] = -9999
        maskDem[maskDem < -100] = 100000  # 先升高再赋空值，不知道为啥np.argmin识别不了-9999

        # 获取该区域的最小值的索引
        min_index_in_region = np.unravel_index(np.argmin(maskDem), maskDem.shape)
        outCell = (min_index_in_region[0], min_index_in_region[1])
        print("最低高程:{:f}".format(maskDem[outCell[0], outCell[1]]))

        # 记录出口点数据
        maskDir = fdir.copy()  # np.zeros((row,col))
        maskDir[maskSink == 1] = 0
        maskDir[outCell[0], outCell[1]] = 0
        maskDem[maskDem == 100000] = dem_nodata

        # 计算淹没内的流向
        # A = priorityFlood(maskSink, 0, maskDem, (outCell[0], outCell[1]), maskDir)  # Python
        A = queuePriorityFlood(maskSink, 0, maskDem, (outCell[0], outCell[1]), maskDir)  # Python

        # 写入流向
        for ti in range(maskrow):
            for tj in range(maskcol):
                if A[ti, tj] in [0, 255]:
                    # print(A[ti,tj])
                    continue
                fdir[ti, tj] = A[ti, tj]
                # print(A[ti,tj])
    else:
        for soid in savedOccCollectionIds:
            print(soid)
            bigCollection = np.zeros((row, col))  # 待处理的洼地掩膜
            soidmask = labledOcc == soid
            initialsinkIds = np.unique(initialbasins[soidmask])  # 记录occ连通体经过的洼地集水区，这些集水区最后合并成一个大的集合进行处理
            tempIds += list(initialsinkIds)
            # 定义需要重新赋值的固定值和新值的映射
            values_to_replace = initialsinkIds  # 需要替换的值
            new_value = 1  # 替换后的值
            # 使用 np.isin 找到需要替换的值并赋值
            bigCollection[np.isin(initialbasins, values_to_replace)] = new_value
            bigCollection[fdir == f_nodata] = 0
            outDEM[bigCollection == 1] = dem_nodata
            # Raster.save_raster(os.path.join(r'F:\青藏高原水体数据集\DATA\TBP_endorheic\173623',"sink_"+str(soid)+".tif"),bigCollection,proj,geo,gdal.GDT_Byte,0)
            # break

            # 处理bigcollection
            cells = np.argwhere(bigCollection == 1)
            minxy = np.amin(cells, axis=0)
            maxxy = np.amax(cells, axis=0)
            minX = minxy[0]
            maxX = maxxy[0]
            minY = minxy[1]
            maxY = maxxy[1]
            # print(minX,minY,maxX,maxY)

            # 掩膜sink
            maskSink = bigCollection[minX:maxX + 1, minY:maxY + 1].copy()
            maskSink[maskSink != 1] = 0  # 写成0 是因为好用：否则后面查最小值的时候会错误
            maskrow, maskcol = maskSink.shape

            maskDem = dem[minX:maxX + 1, minY:maxY + 1].copy()
            maskDem[maskSink == 0] = -9999
            maskDem[maskDem < -100] = 100000  # 先升高再赋空值，不知道为啥np.argmin识别不了-9999

            # 获取该区域的最小值的索引
            min_index_in_region = np.unravel_index(np.argmin(maskDem), maskDem.shape)
            outCell = (min_index_in_region[0], min_index_in_region[1])
            print("最低高程:{:f}".format(maskDem[outCell[0], outCell[1]]))

            # 记录出口点数据
            maskDir = fdir[minX:maxX + 1, minY:maxY + 1].copy()  # np.zeros((row,col))
            maskDir[maskSink == 1] = 0
            maskDir[outCell[0], outCell[1]] = 0

            maskDem[maskDem == 100000] = dem_nodata
            # print(maskDir.shape,maskSink.shape,maskDem.shape)

            # 计算淹没内的流向
            # A = priorityFlood(maskSink, 0, maskDem, (outCell[0], outCell[1]), maskDir)  # Python
            A = queuePriorityFlood(maskSink, 0, maskDem, (outCell[0], outCell[1]), maskDir)  # Python

            # print(A.min(),A.max())
            # 写入流向
            for ti in range(maskrow):
                for tj in range(maskcol):
                    if A[ti, tj] in [0, 255]:
                        # print(A[ti,tj])
                        continue
                    fdir[minX + ti, minY + tj] = A[ti, tj]
                    # print(A[ti,tj])

        print("连通体处理成功")

        temp_outDEMfile = os.path.join(Venu, baseName + "_tempoutdem_.tif")
        temp_outFillfile = os.path.join(Venu, baseName + "_tempoutFill_.tif")
        temp_outtdirfile = os.path.join(Venu, baseName + "_tempouttdir_.tif")

        Raster.save_raster(temp_outDEMfile, outDEM, proj, geo, gdal.GDT_Float32, dem_nodata)

        # 处理外流的洼地:搜索外边界上最低的点流到洼地以外的斜率最大的栅格
        # 正常填洼
        import whitebox
        wbt = whitebox.WhiteboxTools()
        wbt.fill_depressions_planchon_and_darboux(temp_outDEMfile, temp_outFillfile, fix_flats=True)
        wbt.d8_pointer(temp_outFillfile, temp_outtdirfile, esri_pntr=True)
        # sink.Cal_filled(temp_outDEMfile,temp_outFillfile)
        # sink.Cal_dir(temp_outFillfile,temp_outtdirfile)
        newd8 = Raster.get_raster(temp_outtdirfile)
        FillDem = Raster.get_raster(temp_outFillfile)
        # newd8[MASK == 0] = 255
        cells0 = np.argwhere(newd8 == 0)
        modifiedcells = []
        for temp_cell in cells0:
            outdir = [0, math.inf]
            for k in range(8):
                nextCell = (temp_cell[0] + dmove[k][0], temp_cell[1] + dmove[k][1])
                if not check_boundary(nextCell[0], nextCell[1], row, col):
                    continue
                if dem[nextCell[0], nextCell[1]] == dem_nodata:
                    continue
                if outDEM[nextCell[0], nextCell[1]] != dem_nodata:
                    continue
                nextCellH = dem[nextCell[0], nextCell[1]]
                if nextCellH > outdir[1]:
                    continue
                outdir = [2 ** k, nextCellH]
            if outdir[0] == 0:
                modifiedcells.append([temp_cell[0], temp_cell[1]])
            else:
                newd8[temp_cell[0], temp_cell[1]] = outdir[0]
        fdir[outDEM != dem_nodata] = newd8[outDEM != dem_nodata]

        # # 计算洼地上游，从边界上最低的与其他已有流向的栅格流出
        for temp_cell in modifiedcells:
            popCells = [(temp_cell[0], temp_cell[1])]
            sinkCells = [(temp_cell[0], temp_cell[1], FillDem[temp_cell[0], temp_cell[1]])]
            minX = math.inf
            minY = math.inf
            maxX = -1
            maxY = -1
            while popCells:
                popCell = popCells.pop()
                minX = min(minX, popCell[0])
                minY = min(minY, popCell[1])
                maxX = max(maxX, popCell[0])
                maxY = max(maxY, popCell[1])

                for kk in range(8):
                    nextCell = (popCell[0] + dmove[kk][0], popCell[1] + dmove[kk][1])
                    if not check_boundary(nextCell[0], nextCell[1], row, col):
                        continue
                    if outDEM[nextCell[0], nextCell[1]] == dem_nodata:
                        continue
                    nextDir = fdir[nextCell[0], nextCell[1]]
                    if 2 ** ((kk + 4) % 8) != nextDir:
                        continue
                    popCells.append(nextCell)
                    sinkCells.append((nextCell[0], nextCell[1], FillDem[nextCell[0], nextCell[1]]))
            # 构建掩膜计算矩阵，从边界上最低的与其他已有流向的栅格流出
            maxHCell = [-1, -1, -100]  # x ,y, fdir
            sinkCells.sort(key=lambda x: x[2])
            flag = False
            for cell in sinkCells:
                for kk in range(8):
                    nextCell = [cell[0] + dmove[kk][0], cell[1] + dmove[kk][1]]
                    if not check_boundary(nextCell[0], nextCell[1], row, col):
                        continue
                    if dem[nextCell[0], nextCell[1]] == dem_nodata:
                        continue
                    temp_nextCell = [nextCell[0], nextCell[1], FillDem[nextCell[0], nextCell[1]]]
                    if temp_nextCell in sinkCells:
                        continue

                    flag = True
                    maxHCell = [cell[0], cell[1], 2 ** kk]
                    break
                if flag:
                    break

            # 掩膜sink
            maskSink = np.zeros((maxX - minX + 1, maxY - minY + 1))
            for sinkCell in sinkCells:
                maskSink[sinkCell[0] - minX, sinkCell[1] - minY] = 1

            maskSink[maskSink != 1] = 0  # 写成0 是因为好用：否则后面查最小值的时候会错误
            maskrow, maskcol = maskSink.shape
            maskDem = dem[minX:maxX + 1, minY:maxY + 1].copy()
            maskDem[maskSink == 0] = dem_nodata
            # 记录出口点数据
            maskDir = fdir[minX:maxX + 1, minY:maxY + 1].copy()  # np.zeros((row,col))
            maskDir[maskSink == 1] = 0
            outCell = [maxHCell[0] - minX, maxHCell[1] - minY]
            maskDir[outCell[0], outCell[1]] = maxHCell[2]

            # 计算淹没内的流向
            B = queuePriorityFlood(maskSink, 0, maskDem, (outCell[0], outCell[1]), maskDir)  # Python

            # 写入流向
            for ti in range(maskrow):
                for tj in range(maskcol):
                    if B[ti, tj] in [0, 255]:
                        # print(A[ti,tj])
                        continue
                    fdir[minX + ti, minY + tj] = B[ti, tj]
                    # print(A[ti,tj])

        if os.path.exists(temp_outDEMfile):
            os.remove(temp_outDEMfile)
        if os.path.exists(temp_outFillfile):
            os.remove(temp_outFillfile)
        if os.path.exists(temp_outtdirfile):
            os.remove(temp_outtdirfile)
    Raster.save_raster(outDirfile, fdir, proj, geo, gdal.GDT_Byte, f_nodata)
    return fdir

def process_special_endorheic_fdir(Venu):
    """
    计算特殊内流区的流向:
    1、计算掩膜
    2、计算初始流向
    :param Venu:
    :return:
    """
    make_endorheic_mask_fdir(Venu)
    process_endorheic1(Venu)





def sbatch_buchong_endorheic(Venu):

    endorheic_dirs = os.listdir(Venu)
    for endorheic_dir in endorheic_dirs:
        # sbatch_process_endorheic(os.path.join(Venu, endorheic_dir))
        if 5<int(float(endorheic_dir.split('_')[1])) <400:
            sbatch_process_endorheic(os.path.join(Venu,endorheic_dir))
def Check1(venu):
    """
    检查没有改名字的sink，复制到指定文件，以便后续处理
    :param venu:
    :return:
    """
    fillDir = os.path.join(venu,'buchong.1')
    if not os.path.exists(fillDir):
        os.mkdir(fillDir)
        os.chmod(fillDir,0o777)

    result = []
    dirNames = os.listdir(venu)
    popPaths = [os.path.join(venu,dirName) for dirName in dirNames]
    n = 1
    while popPaths:
        popPath = popPaths.pop()
        dirName = os.path.basename(popPath)
        if len(dirName.split('.')) > 1:
            continue
        if len(dirName.split('_')) == 0:
            continue

        dirNames = os.listdir(popPath)
        for dirName1 in dirNames:
            popPaths.append(os.path.join(popPath,dirName1))
        # 后半段代码：复制所有sink的计算文件
        if dirName == 'endorheic':
            shutil.copytree(popPath,
                            os.path.join(fillDir, 'endorheic'+'_'+str(n)))
            n += 1

            # 前半段代码：检查是否没改名1_dem_.tif
            # filepaths = os.listdir(popPath)
            # for filepath in filepaths:
            #     file1 = os.path.join(popPath,filepath,filepath+'_dem.tif')
            #     file2 = os.path.join(popPath, filepath, filepath + '_dir.tif')
            #     file3 = os.path.join(popPath, filepath, filepath + '_mask.tif')
            #     file4 = os.path.join(popPath, filepath, filepath + '_occ.tif')
            #     if os.path.exists(file1) or os.path.exists(file2) or os.path.exists(file3) or os.path.exists(file4):
            #         shutil.copytree(os.path.join(popPath,filepath),os.path.join(fillDir,filepath))



if __name__ == '__main__':
    # lower_water(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\Global_20241114\MergeDEM\N5353W133132_FABDEM_V1-0.tif')
    # r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\Global_20241114\sink\sink_N5353W133132_FABDEM_V1-0.tif'
    # initial_sink_upstream_extract(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\Global_20241114\dir\dir_N5353W133132_FABDEM_V1-0.tif',
    #                               r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\Global_20241114\outsink.tif',
    #                               r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\Global_20241114\outbasin.tif')
    # process_endorheic(r"F:\青藏高原水体数据集\DATA\TBP_endorheic\6402257")
    # delete(r'F:\青藏高原水体数据集\DATA\TBP_endorheic\109524\109524_mask_.tif')

    process_endorheic(r'F:\青藏高原水体数据集\New_DATA\endorheic_point\129150')

    # occ1 = Raster.get_raster(r'F:\青藏高原水体数据集\New_DATA\endorheic_point\129150\129150_occ_.tif')
    # proj,geo,o_nodata = Raster.get_proj_geo_nodata(r'F:\青藏高原水体数据集\New_DATA\endorheic_point\129150\129150_occ_.tif')
    # occ1[occ1 != o_nodata] = 1
    # occ1[occ1 != 1] = 0
    #
    # numbers1, labledOcc1 = check_disjoint_regions(occ1)
    #
    #
    # print(numbers1,labledOcc1)

    pass