# -*- coding: utf-8 -*-
"""
@Time ： 2024/12/28 11:38
@Auth ：
@File ：postprocess.py
@IDE ：PyCharm
"""
import ast
import math
import shutil
import sqlite3
import csv

import numpy as np
import pandas as pd
from scipy import ndimage

import check_counterflow
import genral_functions
import heap_PF_D8
import valid
from genral_functions import *
import Raster
import db
import os
import process_sink
import endorheic_process
import Merge

def VisualCheck(Venu):
    basename = os.path.basename(Venu)
    dbPath = os.path.join(Venu, basename + "_burndb.db")
    if not os.path.exists(dbPath):
        return
    endor2exFile = os.path.join(Venu, basename + "_update_.csv")
    if not os.path.exists(endor2exFile):
        return
    if os.path.exists(os.path.join(Venu, basename + "_ModifiedDir.tif")):
        return
    endorheicPath = os.path.join(Venu, "endorheic")
    if not os.path.exists(endorheicPath):
        return
    # dbPath = r'F:\青藏高原水体数据集\DATA\TBP_endorheic\Asia_4020015090_burndb_.db'
    # endor2exFile = r'F:\青藏高原水体数据集\DATA\visiualCheck\Asia\090.csv'


    # Step1、修改数据库 \删除endorheic文件夹下的文件夹
    conn = sqlite3.connect(dbPath)
    cursor = conn.cursor()

    with open(endor2exFile, "r",encoding = "gbk") as file:
        reader = csv.DictReader(file)
        for row in reader:
            dirPath = os.path.join(endorheicPath,row["record"])
            if os.path.exists(dirPath):
                shutil.rmtree(dirPath)
            record_id = int(row["record"])
            new_age = int(row["newtype"])
            cursor.execute("""
            UPDATE sink
            SET sinkType = ?
            WHERE SinkId = ?
            """, (new_age, record_id))

    # 提交更改
    conn.commit()
    print(dbPath, "修改成功")
    # 关闭连接
    conn.close()

    # Step2、\删除endorheic文件夹下的文件夹

def sbatch_VisualCheck(venu):
    """
    批量处理每个子文件夹
    :param venu:
    :return:
    """

    dirNames = os.listdir(venu)
    for dirName in dirNames:

        dirPath = os.path.join(venu,dirName)
        if len(dirName.split('.')) > 1:
            continue
        # if dirName != 'Africa_1020027430':
        #     continue

        # if os.path.exists(os.path.join(dirPath, dirName + "_ModifiedDir.tif")):
        #     continue

        try:
            VisualCheck(dirPath)
        except Exception as e:
            print(dirName,e)

        try:
            process_sink.exorheic_process(dirPath)
        except Exception as e:
            print(e)

        try:
            if os.path.exists(os.path.join(dirPath, dirName+"_ModifiedDir.tif")):
                continue
            if os.path.exists(os.path.join(dirPath,"endorheic")):
                endorheic_process.sbatch_process_endorheic(os.path.join(dirPath,"endorheic"))
        except Exception as e:
            print(e)

        try:
            # if os.path.exists(os.path.join(dirPath, dirName+"_ModifiedDir.tif")):
            #     continue
            if os.path.exists(os.path.join(dirPath,"endorheic")):
                Merge.merge_endorheic(dirPath)
                # pass
            elif os.path.exists(os.path.join(dirPath,dirName+"_burndb.db")):
                # pass
                # shutil.copyfile(os.path.join(dirPath, dirName+"_outDIR_.tif"),os.path.join(dirPath, dirName+"_ModifiedDir_.tif"))
                A = Raster.get_raster(os.path.join(dirPath, dirName+"_outDIR_.tif"))
                proj,geo,nodata = Raster.get_proj_geo_nodata(os.path.join(dirPath, dirName+"_outDIR_.tif"))
                A[A==nodata] = 255
                Raster.save_raster(os.path.join(dirPath, dirName+"_ModifiedDir.tif"),A,proj,geo,gdal.GDT_Byte,255)
            else:
                A = Raster.get_raster(os.path.join(dirPath, dirName + "_ModifiedDir.tif"))
                proj, geo, nodata = Raster.get_proj_geo_nodata(os.path.join(dirPath, dirName + "_ModifiedDir.tif"))
                A[A == nodata] = 255
                Raster.save_raster(os.path.join(dirPath, dirName + "_ModifiedDir.tif"), A, proj, geo, gdal.GDT_Byte,
                                   255)
        except Exception as e:
            print(e)



# 处理对向流的代码 #
def counterflow_process(fdir_file,outfile):
    """
    处理由于对向流造成的洼地：
    1）找到对向流的栅格并追溯这两个格子的全部上游，全部赋值为1，找到最大连通体
    2）单独储存上游掩膜，流向，dem，
    3）实施breach策略消除：以边界最低点为出口点实施PF+D8算法

    :param fdir_file:
    :return:
    """
    conuterflow = {1:16,2:32,4:64,8:128,16:1,32:2,64:4,128:8}
    fdir = Raster.get_raster(fdir_file)
    proj,geo,f_nodata = Raster.get_proj_geo_nodata(fdir_file)

    row,col = fdir.shape
    vis = np.zeros((row,col),dtype = np.int8)
    result_mask = np.zeros((row,col),dtype = np.int8)
    nn = 1

    baseDir = os.path.dirname(fdir_file)
    venu_name = os.path.basename(fdir_file).split('_')[1]


    for i in range(row):
        for j in range(col):
            if fdir[i,j] == f_nodata:
                continue
            if vis[i,j] == 1:
                continue
            now_fdir = fdir[i,j]
            if now_fdir not in dmove_dic:
                continue
            nextCell = (i+dmove_dic[now_fdir][0],j+dmove_dic[now_fdir][1])
            if not check_boundary(nextCell[0],nextCell[1],row,col):
                continue

            nextfdir = fdir[nextCell[0],nextCell[1]]
            if nextfdir == conuterflow[now_fdir]:
                # 对向流，进行处理
                vis[i,j] = 1
                vis[nextCell[0],nextCell[1]] = 1
                print('存在对向流')
                result_mask[i,j] = nn
                result_mask[nextCell[0],nextCell[1]] = nn
                nn += 1

                minrow = math.inf
                mincol = math.inf
                maxrow = -1
                maxcol = -1

                # popCells = [(i, j)]
                temppopCells = [(i, j),nextCell]
                popCells = [(i, j), nextCell]
                vis[i, j] = 1
                while popCells:
                    popCell = popCells.pop()
                    minrow = min(minrow, popCell[0])
                    mincol = min(mincol, popCell[1])
                    maxrow = max(maxrow, popCell[0])
                    maxcol = max(maxcol, popCell[1])
                    result_mask[popCell[0], popCell[1]] = 1
                    upCells = get_upstreamCells(popCell[0], popCell[1], fdir, row, col)
                    for upCell in upCells:
                        if vis[upCell[0],upCell[1]] == 1:
                            continue
                        popCells.insert(0,upCell)
                        temppopCells.append(upCell)

                # 生成sink掩膜
                mask = np.zeros((maxrow - minrow + 1 + 2, maxcol - mincol + 1 + 2), dtype=np.int8)  # 做1个栅格的缓冲区
                # maskdem = np.zeros((maxrow-minrow+1+2,maxcol-mincol+1+2)) # 做1个栅格的缓冲区
                # maskdem[:,:] = d_nodata
                # maskdem = dem[minrow - 1:maxrow + 2, mincol - 1:maxcol + 2].copy()
                maskdir = np.zeros((maxrow - minrow + 1 + 2, maxcol - mincol + 1 + 2))  # 做1个栅格的缓冲区
                maskdir[:, :] = f_nodata
                for cell in temppopCells:
                    mask[cell[0] - minrow + 1, cell[1] - mincol + 1] = 1
                    # maskdem[cell[0]-minrow+1,cell[1]-mincol+1] = dem[cell[0],cell[1]]
                    maskdir[cell[0] - minrow + 1, cell[1] - mincol + 1] = fdir[cell[0], cell[1]]
                temp_geo = (geo[0] + geo[1] * (mincol - 1), geo[1], geo[2], geo[3] + geo[5] * (minrow - 1), geo[4], geo[5])

                if not os.path.exists(os.path.join(baseDir,venu_name)):
                    os.mkdir(os.path.join(baseDir,venu_name))
                Raster.save_raster(os.path.join(baseDir, venu_name, 'mask_'+str(nn) + '.tif'), mask, proj, temp_geo,gdal.GDT_Byte, 0)
                Raster.save_raster(os.path.join(baseDir,venu_name,'mask_dir_'+str(nn)+'.tif'),maskdir,proj,temp_geo,gdal.GDT_Byte,f_nodata)




                # while popCells:
                #     popCell = popCells.pop()
                #     result_mask[popCell[0],popCell[1]] = 1
                #     upCells = get_upstreamCells(popCell[0],popCell[1],fdir,row,col)
                #     for upCell in upCells:
                #         if vis[upCell[0],upCell[1]] == 1:
                #             continue
                #         popCells.append(upCell)

    Raster.save_raster(outfile,result_mask,proj,geo,gdal.GDT_Byte,0)
def counterflow_process_check(fdir_file):
    """
    处理由于对向流造成的洼地：
    1）找到对向流的栅格并追溯这两个格子的全部上游，全部赋值为1，找到最大连通体
    2）单独储存上游掩膜，流向，dem，
    3）实施breach策略消除：以边界最低点为出口点实施PF+D8算法

    :param fdir_file:
    :return:
    """
    conuterflow = {1:16,2:32,4:64,8:128,16:1,32:2,64:4,128:8}
    fdir = Raster.get_raster(fdir_file)
    proj,geo,f_nodata = Raster.get_proj_geo_nodata(fdir_file)

    row,col = fdir.shape
    vis = np.zeros((row,col),dtype = np.int8)
    for i in range(row):
        for j in range(col):
            if fdir[i,j] == f_nodata:
                continue
            if vis[i,j] == 1:
                continue
            now_fdir = fdir[i,j]
            if now_fdir not in dmove_dic:
                continue
            nextCell = (i+dmove_dic[now_fdir][0],j+dmove_dic[now_fdir][1])
            if not check_boundary(nextCell[0],nextCell[1],row,col):
                continue

            nextfdir = fdir[nextCell[0],nextCell[1]]
            if nextfdir == conuterflow[now_fdir]:
                # 对向流，进行处理
                vis[i,j] = 1
                vis[nextCell[0],nextCell[1]] = 1
                print('存在对向流',nextfdir , now_fdir)
def counterflow_get_upstream(csv_file,fdr_file):
    A = pd.read_excel(csv_file)

    paras = []
    for text1 in A['Name']:
        coords_str = text1.split(".tif")[-1].strip()
        coords = ast.literal_eval(coords_str)

        paras += coords

    fdr = Raster.get_raster(fdr_file)
    proj,geo,f_nodata = Raster.get_proj_geo_nodata(fdr_file)
    row,col = fdr.shape
    huan_cell = []
    for cell in paras:
        # 检查每个环对应的流向，剔除内流终点
        if fdr[cell[0],cell[1]] == 0:
            continue
        huan_cell.append(cell)
    print(huan_cell)
    temp_position = [['lon','lat']]
    for start_cell in huan_cell:
        lon = geo[0]+start_cell[1]*geo[1]
        lat = geo[3]+start_cell[0]*geo[5]
        temp_position.append([lon,lat])
    with open("/datanode05/zhangbin/FAB_hydrography/CODE/SA_points.csv",'w') as f:
        writer = csv.writer(f)
        writer.writerows(temp_position)
        f.close()
    # basin_counter = np.zeros((row,col),dtype = np.int8)
    # vis = np.zeros((row,col),dtype = np.int8)
    # for start_cell in huan_cell:
    #
    #     pop_cells = [start_cell]
    #     #
    #     # minrow = math.inf
    #     # mincol = math.inf
    #     # maxrow = -1
    #     # maxcol = -1
    #     temppopCells = []
    #     while pop_cells:
    #
    #         pop_cell = pop_cells.pop()
    #         # minrow = min(minrow, pop_cell[0])
    #         # mincol = min(mincol, pop_cell[1])
    #         # maxrow = max(maxrow, pop_cell[0])
    #         # maxcol = max(maxcol, pop_cell[1])
    #         temppopCells.append(pop_cell)
    #         # pops.append(pop_cell)
    #         now_dir = fdr[pop_cell[0], pop_cell[1]]
    #
    #         if vis[pop_cell[0], pop_cell[1]] != 0:
    #             continue
    #         if now_dir not in dmove_dic:
    #             break
    #         vis[pop_cell[0], pop_cell[1]] = 1
    #         basin_counter[pop_cell[0], pop_cell[1]] = 1
    #
    #         upcells = genral_functions.get_upstreamCells(pop_cell[0], pop_cell[1], fdr, row, col)
    #         for upcell in upcells:
    #             if vis[upcell[0], upcell[1]] == 0:
    #                 pop_cells.append(upcell)
    # Raster.save_raster("/datanode05/zhangbin/FAB_hydrography/run_data/SouthAmerica/huan_counterflow.tif",basin_counter,proj,geo,gdal.GDT_Byte,0)
def divide_conuterflow_region(counterflow_file,dem_file,fdir_file,venu):
    """
    拆分对流栅格的上游，按照每个文件夹“dem、dir、sink”的方式组织保存
    :param counterflow_file:
    :param venu:
    :return:
    """
    if not os.path.exists(venu):
        os.mkdir(venu)
        # os.chmod(venu,0o777)
    conuterflow = Raster.get_raster(counterflow_file)
    proj,geo,c_nodata = Raster.get_proj_geo_nodata(counterflow_file)

    dem = Raster.get_raster(dem_file)
    proj,geo,d_nodata = Raster.get_proj_geo_nodata(dem_file)

    fdir = Raster.get_raster(fdir_file)
    proj,geo,f_nodata = Raster.get_proj_geo_nodata(fdir_file)

    row,col = conuterflow.shape

    vis = np.zeros((row,col),dtype = np.int8)
    ids = 1
    for i in range(row):
        for j in range(col):
            if conuterflow[i,j] == c_nodata:
                continue
            if vis[i,j] == 1:
                continue
            # 寻找最大连通体
            # 首先找到同id的最大单元连通体
            # 记录四至
            minrow = math.inf
            mincol = math.inf
            maxrow = -1
            maxcol = -1

            popCells = [(i, j)]
            temppopCells = [(i, j)]
            vis[i, j] = 1
            while popCells:
                popCell = popCells.pop()
                minrow = min(minrow,popCell[0])
                mincol = min(mincol,popCell[1])
                maxrow = max(maxrow,popCell[0])
                maxcol = max(maxcol,popCell[1])

                for k in range(8):
                    nextCell = (popCell[0] + dmove[k][0], popCell[1] + dmove[k][1])
                    if not (0 <= nextCell[0] < row and 0 <= nextCell[1] < col):
                        continue
                    if vis[nextCell[0], nextCell[1]] == 1:
                        continue
                    if conuterflow[nextCell[0], nextCell[1]] == c_nodata:
                        continue
                    if conuterflow[nextCell[0], nextCell[1]] != conuterflow[popCell[0], popCell[1]]:
                        continue
                    popCells.append(nextCell)
                    temppopCells.append(nextCell)
                    vis[nextCell[0], nextCell[1]] = 1
            # 生成sink掩膜
            mask = np.zeros((maxrow-minrow+1+2,maxcol-mincol+1+2),dtype = np.int8) # 做1个栅格的缓冲区
            # maskdem = np.zeros((maxrow-minrow+1+2,maxcol-mincol+1+2)) # 做1个栅格的缓冲区
            # maskdem[:,:] = d_nodata
            maskdem = dem[minrow-1:maxrow+2,mincol-1:maxcol+2].copy()
            maskdir = np.zeros((maxrow-minrow+1+2,maxcol-mincol+1+2)) # 做1个栅格的缓冲区
            maskdir[:,:] = f_nodata
            for cell in temppopCells:
                mask[cell[0]-minrow+1,cell[1]-mincol+1] = 1
                # maskdem[cell[0]-minrow+1,cell[1]-mincol+1] = dem[cell[0],cell[1]]
                maskdir[cell[0]-minrow+1,cell[1]-mincol+1] = fdir[cell[0],cell[1]]
            temp_geo = (geo[0] + geo[1] * (mincol-1), geo[1], geo[2], geo[3] + geo[5] * (minrow-1), geo[4], geo[5])

            temp_path = os.path.join(venu,str(ids))
            if not os.path.exists(temp_path):
                os.mkdir(temp_path)
                os.chmod(temp_path,0o777)

            outsinkfile = os.path.join(temp_path,'sink_'+str(ids)+'_.tif')
            outdemfile = os.path.join(temp_path,'dem_'+str(ids)+'_.tif')
            outdirfile = os.path.join(temp_path,'dir_'+str(ids)+'_.tif')
            Raster.save_raster(outsinkfile,mask,proj,temp_geo,gdal.GDT_Byte,0)
            Raster.save_raster(outdemfile,maskdem,proj,temp_geo,gdal.GDT_Float32,d_nodata)
            Raster.save_raster(outdirfile,maskdir,proj,temp_geo,gdal.GDT_Byte,f_nodata)
            # os.chmod(outsinkfile,0o777)
            # os.chmod(outdemfile,0o777)
            # os.chmod(outdirfile,0o777)
            ids += 1

def sbatch_PFD8_exorheic(mask_dir):
    import special_process

    paras = []
    big_files = []
    dirPaths = os.listdir(mask_dir)
    for dirPath in dirPaths:
        if len(dirPath.split('.')) > 1:
            continue
        sinkfile = os.path.join(mask_dir,dirPath,'sink_'+dirPath+'_.tif')
        demfile = os.path.join(mask_dir,dirPath,'dem_'+dirPath+'_.tif')
        dirfile = os.path.join(mask_dir,dirPath,'dir_'+dirPath+'_.tif')
        outpath = os.path.join(mask_dir,dirPath,'outdir_'+dirPath+'_.tif')
        size_bytes = os.path.getsize(demfile)  # 文件大小（字节）
        size_mb = size_bytes / (1024 * 1024)  # 转换为 MB
        if size_mb > 500:
            print(demfile)
            big_files.append(demfile)
            continue
        paras.append([demfile,sinkfile,dirfile,outpath])
        # special_process.PFD8_exorheic(demfile,sinkfile,dirfile,outpath)

    po = Pool(10)
    for para in paras:
        # po.apply_async(special_process.PFD8_exorheic,(para[0],para[1],para[2],para[3],))
        po.apply_async(special_process.modify_conuter_flow, (para[0], para[1], para[2], para[3],))
        # special_process.modify_conuter_flow(para[0], para[1], para[2], para[3])


    po.close()
    po.join()
    print(big_files)
def check_all_0(mask_dir):
    """
    检查生成的流向是不是0，如果是，则删除，不是则保留
    :param mask_dir:
    :return:
    """
    dirPaths = os.listdir(mask_dir)
    for dirPath in dirPaths:
        outpath = os.path.join(mask_dir, dirPath, 'outdir_' + dirPath + '_.tif')
        if not os.path.exists(outpath):
            print(outpath)
            continue
        # temp = Raster.get_raster(outpath)
        # proj,geo,nodata = Raster.get_proj_geo_nodata(outpath)
        # temp[temp == 0] = 255
        # Raster.save_raster(outpath,temp,proj,geo,gdal.GDT_Byte,nodata)




def merge_counterflow_outdir(mask_dir,outfdir_file,out_fdr_file):
    """
    将修复后的对向流栅格上游流向合并至原始流向进行修正。
    :param mask_dir:
    :param outfdir_file:
    :return:
    """

    paras = [outfdir_file]
    dirPaths = os.listdir(mask_dir)
    for dirPath in dirPaths:
        if len(dirPath.split('.')) > 1:
            continue
        dirfile = os.path.join(mask_dir,dirPath,'modified_dir_'+dirPath+'_.tif')
        if os.path.exists(dirfile):
            paras.append(dirfile)
    #
    # paras = [outfdir_file,"/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正缺失的海/1/modified_dir_1.tif",
    #          "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正缺失的海/2/modified_dir_2.tif","/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正缺失的海/3/modified_dir_3.tif",
    #          "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正缺失的海/2/modified_dir_21_.tif","/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正缺失的海/3/modified_dir_31_.tif"]
    venu = os.path.dirname(outfdir_file)


    # 构造 BuildVRTOptions（正确方式）
    vrt_options = gdal.BuildVRTOptions(
        separate=False,
        resolution='highest',
        VRTNodata=255  # 设定 NoData 值，可根据你的影像情况修改
    )

    # 构建 VRT
    vrt_path = os.path.join(venu, "temp_merged6.vrt")
    if os.path.exists(vrt_path):
        os.remove(vrt_path)
    vrt_ds = gdal.BuildVRT(vrt_path, paras, options=vrt_options)
    vrt_ds = None  # 保存并关闭

    # 输出 GeoTIFF + 压缩
    # os.path.join(venu, "Final_FDIR.tif")
    gdal.Translate(
        out_fdr_file,
        vrt_path,
        format="GTiff",
        creationOptions=[
            "COMPRESS=LZW",
            "TILED=YES",
            "BIGTIFF=IF_SAFER"
        ]
    )
    os.remove(vrt_path)
def F2(sink_file,dir_file,dem_file,outdir_file):
    """
    简易算法，把对流的栅格纠正为流向不为上游的最大坡度下游，
    :param sink_file:
    :param dir_file:
    :param dem_file:
    :param outdir_file:
    :return:
    """
    sink = Raster.get_raster(sink_file)

    fdir = Raster.get_raster(dir_file)

    dem = Raster.get_raster(dem_file)

    conuterflow = {1: 16, 2: 32, 4: 64, 8: 128, 16: 1, 32: 2, 64: 4, 128: 8}

    proj, geo, f_nodata = Raster.get_proj_geo_nodata(dir_file)
    correct_cells = []
    row, col = fdir.shape
    vis = np.zeros((row, col), dtype=np.int8)

    new_fdir = fdir.copy()
    for i in range(row):
        for j in range(col):
            if fdir[i, j] == f_nodata:
                continue
            if vis[i, j] == 1:
                continue
            now_fdir = fdir[i, j]
            if now_fdir not in dmove_dic:
                continue
            nextCell = (i + dmove_dic[now_fdir][0], j + dmove_dic[now_fdir][1])
            if not check_boundary(nextCell[0], nextCell[1], row, col):
                continue

            nextfdir = fdir[nextCell[0], nextCell[1]]
            if nextfdir == conuterflow[now_fdir]:
                # 对向流，进行处理
                vis[i, j] = 1
                vis[nextCell[0], nextCell[1]] = 1
                print('存在对向流', nextfdir, now_fdir)
                correct_cells.append((i,j))

    for cell in correct_cells:
        now_dir = fdir[cell[0],cell[1]]
        flag = True
        for k in range(8):
            # 流向一个不流向他的
            next_cell = (cell[0] + dmove[k][0], cell[1] + dmove[k][1])
            next_dir = fdir[next_cell[0],next_cell[1]]
            if next_dir != conuterflow[now_dir]:
                new_fdir[cell[0], cell[1]] = 2 ** k
                flag = False
                break
        if flag:
            print(cell,'被包围了')

    Raster.save_raster(outdir_file,new_fdir,proj,geo,gdal.GDT_Byte,f_nodata)


def process_huan(fdir_file):

    fdir = Raster.get_raster(fdir_file)
    proj,geo,f_nodata = Raster.get_proj_geo_nodata(fdir_file)

    fdir[4168,283] = 1

    Raster.save_raster(fdir_file,fdir,proj,geo,gdal.GDT_Byte,f_nodata)
def modify_huan(fdir_file):
    """
    构建树索引图，来更新汇流累积量时加快速度
    1) 初始化路径
    2） 从边界出发，寻找下游：如果没有下游了停止，则更新此条路径，赋予新id；如果中途汇入有id的栅格，则停止赋予此id
    :param fdir_file:
    :param output_file:
    :return:
    """


    def cal_next_code(now_cell,now_fdir,now_fdir_name,row_num,col_num):
        """
        构建表间汇流。
        1）判断next_cell的位置（上下左右栅格以及所在区域）
        2)返回信息
        :param next_cell:
        :param now_fdir_file:
        :return:
        """
        infos = now_fdir_name.split('_')  # [fdir,row_num,col_num,region] # 区域编码的行列号，不同于参数中的栅格行列号
        now_row = int(float(infos[3]) // float(infos[2]))
        now_col = int(float(infos[3]) % float(infos[2]))
        new_row,new_col = int(float(infos[1])+1),int(float(infos[2])+1)

        if now_cell[0] == 0 and now_cell[1] == col_num-1:
            # 上边界
            if now_fdir in [1, 2]:
                # 右
                if now_col == float(infos[2])-1:
                    # 区域边界条件，如果超出区域范围，则不记算新的区域id
                    return new_row,new_col,int(float(infos[3]) + 1)
                new_id = int(float(infos[3]) + 1)
                return int(float(infos[1])),int(float(infos[2])),new_id
            if now_fdir in [32, 64]:
                # 上
                if now_row == 0:
                    return new_row,new_col,int((now_row - 1)*float(infos[2]) + now_col)
                new_id = int((now_row - 1)*float(infos[2]) + now_col)
                return int(float(infos[1])),int(float(infos[2])),new_id
            if now_fdir in [128]:
                # 右上
                if ((now_row == 0) or (now_col == float(infos[2])-1)):
                    return new_row,new_col,int((now_row - 1)*float(infos[2]) + now_col + 1)
                new_id = int((now_row - 1)*float(infos[2]) + now_col + 1)
                return int(float(infos[1])),int(float(infos[2])),new_id

        if now_cell[0] == 0 and now_cell[1] == 0:
            # 上边界
            if now_fdir in [8, 16]:
                # 左
                if now_col == 0:
                    return new_row,new_col,int(float(infos[3]) - 1)
                new_id = int(float(infos[3]) - 1)
                return int(float(infos[1])),int(float(infos[2])),new_id
            if now_fdir in [128, 64]:
                # 上
                if now_row == 0:
                    return new_row,new_col,int((now_row - 1)*float(infos[2]) + now_col)
                new_id = int((now_row - 1)*float(infos[2]) + now_col)
                return int(float(infos[1])),int(float(infos[2])),new_id
            if now_fdir in [32]:
                # 左上
                if ((now_col == 0) or (now_row == 0)):
                    return new_row,new_col,int((now_row - 1)*float(infos[2]) + now_col - 1)
                new_id = int((now_row - 1)*float(infos[2]) + now_col - 1)
                return int(float(infos[1])),int(float(infos[2])),new_id

        if now_cell[0] == row_num-1 and now_cell[1] == 0:
            # 下边界
            if now_fdir in [32, 16]:
                # 左
                if now_col == 0:
                    return new_row,new_col,int(float(infos[3]) - 1)
                new_id = int(float(infos[3]) - 1)
                return int(float(infos[1])),int(float(infos[2])),new_id
            if now_fdir in [2, 4]:
                # 下
                if now_row == float(infos[1])-1:
                    return new_row,new_col,int((now_row + 1)*float(infos[2]) + now_col)
                new_id = int((now_row + 1)*float(infos[2]) + now_col)
                return int(float(infos[1])),int(float(infos[2])),new_id
            if now_fdir in [8]:
                # 左下
                if ((now_col == 0) or (now_row == float(infos[1])-1)):
                    return new_row,new_col,int((now_row + 1)*float(infos[2]) + now_col - 1)
                new_id = int((now_row + 1)*float(infos[2]) + now_col - 1)
                return int(float(infos[1])),int(float(infos[2])),new_id

        if now_cell[0] == row_num-1 and now_cell[1] == col-1:
            # 下边界
            if now_fdir in [1, 128]:
                # 右
                if now_col == float(infos[2])-1:
                    return new_row,new_col,int(float(infos[3]) + 1)
                new_id = int(float(infos[3]) + 1)
                return int(float(infos[1])),int(float(infos[2])),new_id
            if now_fdir in [8, 4]:
                # 下
                if now_row == float(infos[1])-1:
                    return new_row,new_col,int((now_row + 1)*float(infos[2]) + now_col)
                new_id = int((now_row + 1)*float(infos[2]) + now_col)
                return int(float(infos[1])),int(float(infos[2])),new_id
            if now_fdir in [2]:
                # 右下
                if ((now_col == float(infos[2])-1) or (now_row == float(infos[1])-1)):
                    return new_row,new_col,int((now_row + 1)*float(infos[2]) + now_col + 1)
                new_id = int((now_row + 1)*float(infos[2]) + now_col + 1)
                return int(float(infos[1])),int(float(infos[2])),new_id

        if now_cell[0] == 0 :
            # 上边界
            if now_fdir in [32, 64, 128]:
                # 上
                if now_row == 0:
                    return new_row,new_col,int((now_row - 1)*float(infos[2]) + now_col)
                new_id = int((now_row - 1)*float(infos[2]) + now_col)
                return int(float(infos[1])),int(float(infos[2])),new_id
        if now_cell[0] == row_num-1 :
            # 下边界
            if now_fdir in [2,4,8]:
                # 下
                if now_row == float(infos[1])-1:
                    return new_row,new_col,int((now_row + 1)*float(infos[2]) + now_col)
                new_id = int((now_row + 1)*float(infos[2]) + now_col)
                return int(float(infos[1])),int(float(infos[2])),new_id
        if now_cell[1] == 0 :
            # 左边界
            if now_fdir in [16, 32,8]:
                # 左
                if now_col == 0:
                    return new_row,new_col,int((now_row)*float(infos[2]) + now_col-1)
                new_id = int((now_row)*float(infos[2]) + now_col-1)
                return int(float(infos[1])),int(float(infos[2])),new_id
        if now_cell[1] == col-1 :
            # 右边界
            if now_fdir in [1,2,128]:
                # 右
                if now_col == float(infos[2])-1:
                    return new_row,new_col,int((now_row)*float(infos[2]) + now_col+1)
                new_id = int((now_row)*float(infos[2]) + now_col+1)
                return int(float(infos[1])),int(float(infos[2])),new_id

    baseName = os.path.basename(fdir_file)
    infos1 = baseName.split('_')
    fdir = Raster.get_raster(fdir_file)
    proj,geo,f_nodata = Raster.get_proj_geo_nodata(fdir_file)
    row,col = fdir.shape

    mask = np.zeros((row,col), dtype=bool)
    mask[fdir != f_nodata] = True
    # 记录边界shange
    # 找边缘：腐蚀后与原始掩膜做异或（边界为 True）
    eroded = ndimage.binary_erosion(mask)
    boundary = mask ^ eroded  # XOR 得到边界像元
    # print(boundary)
    print(baseName)

    # 获取边界像元的行列号
    rows, cols = np.where(boundary)
    # 输出

    tree_table = {}  # 记录表内汇流
    id_cell = {}

    table_cell_id = {}
    table_cell = {}  # 记录表间连接关系

    vis = np.zeros((row, col))
    vis[:, :] = -9999
    road_id = 1


    upvis = np.zeros((row,col))
    up_id = 1
    for r, c in zip(rows, cols):
        huan_vis = np.zeros((row, col))  # 用来检测环
        # print(f"边界像元：行 {r}, 列 {c}")
        r = int(r)
        c = int(c)
        if vis[r,c] != -9999:

            continue
        pop_cells = [(r, c)]
        road_cells = [(r, c)]

        case = 0
        while pop_cells:
            pop_cell = pop_cells.pop()
            huan_vis[pop_cell[0],pop_cell[1]] = 1
            # print(pop_cell)
            now_fdir = fdir[pop_cell[0], pop_cell[1]]
            if now_fdir not in dmove_dic:
                case = 1
                continue
            next_cell = (pop_cell[0] + dmove_dic[now_fdir][0], pop_cell[1] + dmove_dic[now_fdir][1])
            if not (0 <= next_cell[0] < row and 0 <= next_cell[1] < col):
                # 超出边界时的计算，要计算相邻的区域id
                new_row, new_col, next_id = cal_next_code(pop_cell, now_fdir, baseName, row, col)
                case = 1
                continue

            if vis[next_cell[0],next_cell[1]] != -9999:
                # 中止寻找下游，赋予前路栅格该id
                temp_id = vis[next_cell[0],next_cell[1]]
                for temp_cell in road_cells:
                    vis[temp_cell[0],temp_cell[1]] = temp_id
                continue

            if fdir[next_cell[0],next_cell[1]] == f_nodata:
                # 没超出边界，但是nextCEll是空值，则记录当前file，和其他的信息
                id_cell.setdefault(road_id, pop_cell)
                tree_table.setdefault((baseName, r, c), (baseName, pop_cell[0], pop_cell[1]))
                # new_row, new_col, next_id = cal_next_code(pop_cell, now_fdir, baseName, row, col)
                # if next_id != -9999:
                next_file = infos1[0] + '_' + infos1[1] + '_' + infos1[2] + '_' + infos1[3] + '_.tif'
                # next_file = infos1[0] + '_' + infos1[1] + '_' + infos1[2] + '_' + str(next_id) + '_.tif'
                case = 1
                continue
            if huan_vis[next_cell[0],next_cell[1]] == 1:
                print(baseName,'有环',"需要修正流向栅格坐标为",next_cell,geo)
                road_cells.append(next_cell)
                # print(road_cells)

                # # 使用 np.where 找到所有下标
                indices = [i for i, x in enumerate(road_cells) if x == next_cell]
                #
                # # print(indices)  # 输出 [1 3 5]
                print(road_cells[indices[0]:indices[-1]+1])
                # correct_cells = road_cells[indices[0]:indices[-1]].copy()
                # # 寻找他们的上游
                # new_popcells = correct_cells.copy()
                # while new_popcells:
                #     new_pop_cell = new_popcells.pop()
                #     upvis[new_pop_cell[0],new_pop_cell[1]] = up_id
                #     temp_results = get_upstreamCells(new_pop_cell[0],new_pop_cell[1],fdir,row,col)
                #     for ii in temp_results:
                #         if upvis[ii[0],ii[1]] != up_id:
                #             new_popcells.append(ii)
                #
                # # 查找该循环里非本up_id的流向
                # ttflag = True
                # for jj in correct_cells:
                #     for kkk in range(8):
                #         next_jk = (jj[0]+dmove[kkk][0],jj[1]+dmove[kkk][1])
                #         if upvis[next_jk[0],next_jk[1]] != up_id:
                #             fdir[jj[0],jj[1]] = 2**kkk
                #             ttflag = False
                #             break
                # if ttflag:
                #     print('环依然存在,失败')
                # up_id += 1

                # return
                break
            pop_cells.append(next_cell)
            road_cells.append(next_cell)

        if case == 1:
            for temp_cell in road_cells:
                vis[temp_cell[0], temp_cell[1]] = road_id
            road_id += 1

        # Raster.save_raster(fdir_file,fdir,proj,geo,gdal.GDT_Byte,f_nodata)



def test_local_huan(start_cells,fdr_file,outfile):
    """
    构建两个流向图层，一个为原始的流向，另一个为修正的流向

    1）找到环的上游集水区；
    2）以环内的cell为起始点，按bfs寻找最低的边界点，记录这个点
    3）将最近的上游点指向这个点，并沿原来的流向依次迭代，将流向指向每次迭代的上游，直到到达环内像元。

    所需数据：流向、环掩膜、上游集水区
    :param fdr_file:
    :return:
    """

    # start_cells = [(102, 3753), (3, 3678), (5, 3688), (7, 3696), (21, 3756), (32, 3747), (128, 3912), (164, 3988), (209, 4098), (260, 4275), (265, 4304), (297, 4558), (309, 4756), (314, 4901), (1, 3671), (16, 3733), (5, 3685), (9, 3705), (11, 3714), (13, 3722), (10, 3708), (12, 3718), (18, 3738), (16, 3732), (57, 3810), (83, 3841), (54, 3756), (58, 3755), (44, 3754), (45, 3755), (52, 3756), (102, 3867), (58, 3812), (56, 3755), (63, 3753), (65, 3820), (154, 3963), (110, 3879), (73, 3750), (101, 3865), (78, 3750), (149, 3952), (85, 3766), (85, 3746), (83, 3747), (120, 3897), (90, 3763), (91, 3746), (102, 3754), (95, 3760), (100, 3861), (96, 3744), (93, 3746), (108, 3876), (101, 3862), (102, 3757), (137, 3929), (119, 3895), (116, 3887), (116, 3889), (129, 3914), (132, 3919), (122, 3901), (123, 3903), (126, 3908), (173, 4006), (131, 3917), (203, 4077), (136, 3925), (139, 3932), (145, 31774), (183, 4028), (162, 3981), (140, 3934), (160, 3975), (161, 3978), (151, 3956), (153, 3961), (152, 3959), (170, 3998), (172, 4003), (225, 4145), (167, 3994), (178, 4016), (251, 4237), (173, 4005), (178, 4017), (193, 4052), (190, 4042), (177, 4013), (195, 4057), (233, 4171), (186, 4034), (187, 4037), (191, 4045), (288, 4456), (207, 4090), (242, 4201), (204, 4083), (197, 4062), (196, 4060), (206, 4089), (246, 4221), (203, 4080), (306, 4679), (219, 4127), (231, 4163), (214, 4113), (217, 4122), (218, 4125), (222, 4135), (232, 4167), (224, 4138), (231, 4165), (231, 4164), (307, 4696), (239, 4193), (239, 4194), (245, 4210), (283, 4408), (240, 4196), (248, 4223), (243, 4207), (255, 4257), (243, 4208), (243, 4209), (278, 4379), (254, 4250), (305, 4664), (249, 4232), (252, 4244), (251, 4233), (285, 4430), (252, 4247), (257, 4263), (313, 4855), (314, 4910), (280, 4390), (284, 4419), (257, 4268), (275, 4362), (261, 4286), (295, 4538), (303, 4644), (265, 4301), (314, 4933), (267, 4318), (269, 4330), (273, 4352), (273, 4355), (276, 4370), (304, 4651), (277, 4376), (279, 4387), (291, 4490), (309, 4747), (279, 4386), (312, 4835), (281, 4397), (294, 4516), (299, 4580), (302, 4620), (281, 4399), (285, 4436), (285, 4437), (285, 4438), (301, 4599), (290, 4482), (289, 4470), (288, 4463), (313, 4840), (287, 4447), (293, 4517), (298, 4569), (302, 4626), (292, 4503), (290, 4473), (293, 4520), (296, 4533), (293, 4506), (313, 4864), (295, 4535), (301, 4603), (308, 4719), (309, 4733), (314, 4944), (302, 4631), (302, 4615), (303, 4647), (311, 4826), (302, 4613), (307, 4685), (309, 4755), (314, 4884), (314, 4896), (308, 4728), (308, 4739), (308, 4725), (314, 4892), (313, 4839), (311, 4775), (311, 4797), (311, 4805), (314, 4940), (312, 4829), (314, 4849), (313, 4863), (313, 4921), (314, 4879), (314, 4937), (314, 4942)]
    fdr = Raster.get_raster(fdr_file)
    proj,geo,nodata = Raster.get_proj_geo_nodata(fdr_file)

    row,col = fdr.shape
    new_fdr = fdr.copy()
    # while True:
    #     if len(start_cells) == 0:
    #         break
    new_fdr = fdr.copy()
    # -------------------------------------------------------------------------- 1、提取流域 -------------------------------------------------------------------------------
    vis = np.zeros((row,col))
    n = 1
    id_start_cell = {}
    for start_cell in start_cells:

        pop_cells = [start_cell]
        # print(start_cell)
        id_start_cell.setdefault(n,start_cell)
        pops = []
        while pop_cells:

            pop_cell = pop_cells.pop()
            pops.append(pop_cell)
            now_dir = fdr[pop_cell[0],pop_cell[1]]
            if vis[pop_cell[0],pop_cell[1]] != 0:
                continue
            if now_dir not in dmove_dic:
                break
            vis[pop_cell[0], pop_cell[1]] = n

            upcells = genral_functions.get_upstreamCells(pop_cell[0],pop_cell[1],fdr,row,col)
            for upcell in upcells:
                if vis[upcell[0],upcell[1]] == 0:
                    pop_cells.append(upcell)
        n += 1

    print('1')
    # Raster.save_raster(r'F:\青藏高原水体数据集\New_DATA\example\NA\fdir_21_31_134__basin1.tif',vis,proj,geo,gdal.GDT_Float32,0)
    # -------------------------------------------------------------------------- 1、提取流域 -------------------------------------------------------------------------------

    # -------------------------------------------------------------------------- 2、寻找环 -------------------------------------------------------------------------------

    vis1 = np.zeros((row, col))  # 环
    all_huan = []
    for start_cell in start_cells:
        temp_huan = []
        pop_cells = [start_cell]
        n = 1
        while pop_cells:
            pop_cell = pop_cells.pop()

            now_dir = fdr[pop_cell[0], pop_cell[1]]
            if vis1[pop_cell[0],pop_cell[1]] != 0:
                continue

            if now_dir not in dmove_dic:
                break
            temp_huan.append(pop_cell)
            vis1[pop_cell[0], pop_cell[1]] = 1
            next_cell = (pop_cell[0]+dmove_dic[now_dir][0],pop_cell[1]+dmove_dic[now_dir][1])
            if vis1[next_cell[0],next_cell[1]] == 0:
                pop_cells.append(next_cell)
        if len(temp_huan) != 0:
            all_huan.append(temp_huan)
    # Raster.save_raster(r'F:\青藏高原水体数据集\New_DATA\example\NA\fdir_21_31_134_mask_basin1.tif', vis1, proj, geo,
    #                    gdal.GDT_Float32, 0)
    print('2')
    # -------------------------------------------------------------------------- 2、寻找环 -------------------------------------------------------------------------------

    id_up_sinks = {}  # 记录流向本sink的上游sink，不让生成新的环路
    id_down_sinks = {}  # 记录流向本sink的下游sink，不让生成新的环路

    for huan in all_huan:
        # -------------------------------------------------------------------------- 3、寻找最近出口点 -------------------------------------------------------------------------------

        now_id = vis[huan[0][0],huan[0][1]]
        id_up_sinks.setdefault(now_id,[])  # 生成空表：初始化
        id_down_sinks.setdefault(now_id, [-1])  # 生成空表：初始化
        vis2 = np.zeros((row, col), dtype=np.int8)  # 环
        flag = False
        outlet = (-1,-1)
        first_cell = [-1,-1]
        # print(huan)
        # print(now_id)
        while huan:

            pop_cell = huan.pop()
            vis2[pop_cell[0],pop_cell[1]] = 1
            # print(pop_cell)
            for k in range(8):

                next_cell = (pop_cell[0]+dmove[k][0],pop_cell[1]+dmove[k][1])
                next_id = vis[next_cell[0],next_cell[1]]

                # 如果是本id
                if next_id == now_id:
                    if vis2[next_cell[0],next_cell[1]] == 0:
                        huan.insert(0,next_cell)
                        vis2[next_cell[0], next_cell[1]] = 1
                else:
                    # 如果这个id不是本id
                    id_up_sinks.setdefault(next_id,[])
                    id_down_sinks.setdefault(next_id,[-1])
                    # 如果id的下游不是本id，
                    if id_down_sinks[next_id][-1] != now_id:
                        # print(id_down_sinks[next_id][-1] , now_id)
                        # 并且也不是本id的上游：停止
                        up_flag = True
                        # --------------------- 所有的上游 --------------------------
                        for temp_id in id_up_sinks[now_id]:
                            if next_id == temp_id:
                                # 假的
                                up_flag = False
                                break
                            for all_up_id in id_up_sinks[temp_id]:
                                if next_id == all_up_id:
                                    # 假的
                                    up_flag = False
                                    break
                            if not up_flag:
                                break
                        # --------------------- 所有的上游 --------------------------
                        if up_flag:
                            # ok了，结束
                            outlet = next_cell
                            new_fdr[pop_cell[0], pop_cell[1]] = 2 ** k
                            first_cell = pop_cell
                            flag = True
                            # 更新上下游关系
                            # print(now_id,next_id,id_down_sinks[next_id],id_up_sinks[next_id],id_up_sinks[now_id],pop_cell,next_cell)
                            id_down_sinks[now_id].append(next_id)
                            id_up_sinks[next_id].append(now_id)
                            id_up_sinks[next_id] += id_up_sinks[now_id]
                            # print(now_id, next_id, id_down_sinks[next_id], id_up_sinks[next_id],
                            #       id_up_sinks[now_id], pop_cell, next_cell)
                            id_start_cell.pop(now_id)
                            break
                    # outlet = next_cell
                    # new_fdr[pop_cell[0],pop_cell[1]] = 2**k
                    # first_cell = pop_cell
                    # flag = True
                    # break
            if flag:
                break
        # -------------------------------------------------------------------------- 3、寻找最近出口点 -------------------------------------------------------------------------------

        # -------------------------------------------------------------------------- 4、回溯修正流向 -------------------------------------------------------------------------------
        reverse_fdr = {1:16,2:32,4:64,8:128,16:1,32:2,64:4,128:8}  # 逆向流
        pop_cells = [first_cell]
        while pop_cells:

            pop_cell = pop_cells.pop()
            if pop_cell == [-1,-1]:
                print("-1 -1 :",pop_cell,huan,now_id,fdr[pop_cell[0],pop_cell[1]])
            now_dir = fdr[pop_cell[0],pop_cell[1]]
            if now_dir not in dmove_dic:
                break
            next_cell = [pop_cell[0]+dmove_dic[now_dir][0],pop_cell[1]+dmove_dic[now_dir][1]]
            new_fdr[next_cell[0],next_cell[1]] = reverse_fdr[now_dir]
            pop_cells.append(next_cell)
            if vis1[next_cell[0],next_cell[1]] != 0:
                break
    start_cells = list(id_start_cell.values())
    fdr = new_fdr.copy()
    print(start_cells)
        # -------------------------------------------------------------------------- 4、回溯修正流向 -------------------------------------------------------------------------------
    Raster.save_raster(outfile, new_fdr, proj, geo,
                               gdal.GDT_Byte, nodata)

def process_huan_global(start_cells,dem_file,fdr_file,venu):
    """
    根据环上的点，生成文件夹，把环的上游的DEM和流向存到文件夹里
    :param start_cells:
    :param dem_file:
    :param fdr_file:
    :param venu:
    :return:
    """

    # dem = Raster.get_raster(dem_file)
    # proj, geo, d_nodata = Raster.get_proj_geo_nodata(dem_file)

    fdr = Raster.get_raster(fdr_file)
    proj, geo, f_nodata = Raster.get_proj_geo_nodata(fdr_file)

    row, col = fdr.shape

    vis = np.zeros((row, col),dtype = np.uint16)
    n = 1
    # id_start_cell = {}
    n_huan_dic = [] # 记录n和环点的
    for start_cell in start_cells:

        pop_cells = [start_cell]
        # id_start_cell.setdefault(n, start_cell)
        temp_n_cell = [n,start_cell]
        minrow = math.inf
        mincol = math.inf
        maxrow = -1
        maxcol = -1
        temppopCells = []
        while pop_cells:

            pop_cell = pop_cells.pop()
            minrow = min(minrow, pop_cell[0])
            mincol = min(mincol, pop_cell[1])
            maxrow = max(maxrow, pop_cell[0])
            maxcol = max(maxcol, pop_cell[1])
            temppopCells.append(pop_cell)
            # pops.append(pop_cell)
            now_dir = fdr[pop_cell[0], pop_cell[1]]

            if vis[pop_cell[0], pop_cell[1]] != 0:
                continue
            if now_dir not in dmove_dic:
                break
            vis[pop_cell[0], pop_cell[1]] = n

            upcells = genral_functions.get_upstreamCells(pop_cell[0], pop_cell[1], fdr, row, col)
            for upcell in upcells:
                if vis[upcell[0], upcell[1]] == 0:
                    pop_cells.append(upcell)

        # 生成sink掩膜
        temp_n_cell += [minrow,mincol,maxrow,maxcol,start_cell[0] - minrow + 1, start_cell[1] - mincol + 1]
        n_huan_dic.append(temp_n_cell)
        mask = np.zeros((maxrow - minrow + 1 + 2, maxcol - mincol + 1 + 2), dtype=np.int8)  # 做1个栅格的缓冲区
        # maskdem = np.zeros((maxrow-minrow+1+2,maxcol-mincol+1+2)) # 做1个栅格的缓冲区
        # maskdem[:,:] = d_nodata
        # maskdem = dem[minrow - 1:maxrow + 2, mincol - 1:maxcol + 2].copy()

        maskdir = np.zeros((maxrow - minrow + 1 + 2, maxcol - mincol + 1 + 2),dtype = np.uint8)  # 做1个栅格的缓冲区
        maskdir[:, :] = f_nodata

        for cell in temppopCells:
            mask[cell[0] - minrow + 1, cell[1] - mincol + 1] = 1
            # maskdem[cell[0]-minrow+1,cell[1]-mincol+1] = dem[cell[0],cell[1]]
            maskdir[cell[0] - minrow + 1, cell[1] - mincol + 1] = fdr[cell[0], cell[1]]
        temp_geo = (geo[0] + geo[1] * (mincol - 1), geo[1], geo[2], geo[3] + geo[5] * (minrow - 1), geo[4], geo[5])

        temp_path = os.path.join(venu, str(n))
        if not os.path.exists(temp_path):
            os.mkdir(temp_path)
            # os.chmod(temp_path, 0o777)

        outsinkfile = os.path.join(temp_path, 'sink_' + str(n) + '_.tif')
        outdemfile = os.path.join(temp_path, 'dem_' + str(n) + '_.tif')
        outdirfile = os.path.join(temp_path, 'dir_' + str(n) + '_.tif')

        Raster.save_raster(outsinkfile, mask, proj, temp_geo, gdal.GDT_Byte, 0)
        # Raster.save_raster(outdemfile, maskdem, proj, temp_geo, gdal.GDT_Float32, d_nodata)
        Raster.save_raster(outdirfile, maskdir, proj, temp_geo, gdal.GDT_Byte, f_nodata)
        # os.chmod(outsinkfile, 0o777)
        # os.chmod(outdemfile, 0o777)
        # os.chmod(outdirfile, 0o777)

        n += 1
    outcsvfile = os.path.join(venu, 'lookup_' + str(n) + '_.csv')
    with open(outcsvfile,'w') as f:
        writer = csv.writer(f)
        writer.writerows(n_huan_dic)


def sbatch_process_huan(csv_file,dem_file,fdr_file,outvenu):
    """
    修正csv里提供的每个环的起始点环

    :param venu:
    :param csv_file:
    :return:
    """
    import re
    import ast

    if not os.path.exists(outvenu):
        os.mkdir(outvenu)


    A = pd.read_excel(csv_file)

    paras = []
    for text1 in A['Name']:
        coords_str = text1.split(".tif")[-1].strip()
        coords = ast.literal_eval(coords_str)

        paras += coords

    fdr = Raster.get_raster(fdr_file)
    proj,geo,f_nodata = Raster.get_proj_geo_nodata(fdr_file)

    huan_cell = []
    for cell in paras:
        # 检查每个环对应的流向，剔除内流终点
        if fdr[cell[0],cell[1]] == 0:
            continue
        huan_cell.append(cell)
    print(huan_cell)
    process_huan_global(huan_cell,dem_file,fdr_file,outvenu)   # 记录环到每个单独的文件夹，运行成功后可注释掉


    # po = Pool(30)
    # for para in paras:
    #
    #     po.apply_async(test_local_huan,(para[1],para[0],para[2],))
    #
    # po.close()
    # po.join()

def sbatch_test_local_huan(csv_file,dem_file,fdr_file,outvenu,lookup_csv):
    """
    单独处理大型环，可手动迭代使用
    :param csv_file:
    :param dem_file:
    :param fdr_file:
    :param outvenu:
    :return:
    """

    import re
    import ast

    # if not os.path.exists(outvenu):
    #     os.mkdir(outvenu)

    A = pd.read_excel(csv_file)

    paras = []
    for text1 in A['Name']:
        # print(text1)
        coords = ast.literal_eval(text1)

        paras += coords
    print(type(paras),paras)

    for file in paras:
        baseName = os.path.basename(file)



def merge_all_processed_fdr(huan_fdr_venu,region_venu,out_fdr_tif):
    """
    将所有的正常流向和修正后的流向复制到一个文件夹里，用修正后的流向替换掉原来有环的流向
    :param huan_fdr_venu:
    :param region_venu:
    :param out_venu:
    :return:
    """

    paras = []
    for name in os.listdir(region_venu):
        if os.path.exists(os.path.join(huan_fdr_venu,name)):
            paras.append(os.path.join(huan_fdr_venu,name))
        else:
            paras.append(os.path.join(region_venu, name))

    # 构造 BuildVRTOptions（正确方式）
    vrt_options = gdal.BuildVRTOptions(
        separate=False,
        resolution='highest',
        VRTNodata=255  # 设定 NoData 值，可根据你的影像情况修改
    )

    # 构建 VRT

    vrt_path = os.path.join(os.path.dirname(out_fdr_tif), 'mergedem11.vrt')
    vrt_ds = gdal.BuildVRT(vrt_path, paras, options=vrt_options)
    vrt_ds = None  # 保存并关闭

    # 输出 GeoTIFF + 压缩
    gdal.Translate(
        out_fdr_tif,
        vrt_path,
        format="GTiff",
        creationOptions=[
            "COMPRESS=LZW",
            "TILED=YES",
            "BIGTIFF=IF_SAFER"
        ]
    )
    os.remove(vrt_path)
    # gdal.Warp(outfile, files, options=[
    #     "-co", "COMPRESS=LZW",  # 选择LZW压缩
    #     "-co", "BIGTIFF=YES"])  # 允许大文件])
    print("合并完成，输出文件：", out_fdr_tif)

def temp_NA(huan_fdr_venu,region_venu):

    A = ['fdir_21_31_480_.tif','fdir_21_31_511_.tif','fdir_21_31_511_.tif','fdir_21_31_480_.tif']

    paras = []
    for name in A:
        if os.path.exists(os.path.join(huan_fdr_venu,name)):
            paras.append(os.path.join(huan_fdr_venu,name))
        else:
            paras.append(os.path.join(region_venu, name))

    # 构造 BuildVRTOptions（正确方式）
    vrt_options = gdal.BuildVRTOptions(
        separate=False,
        resolution='highest',
        VRTNodata=255  # 设定 NoData 值，可根据你的影像情况修改
    )

    # 构建 VRT
    out_fdr_tif = "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/NorthAmerica/补丁/NA4.tif"
    vrt_path = os.path.join(os.path.dirname(out_fdr_tif), 'mergedem11.vrt')
    vrt_ds = gdal.BuildVRT(vrt_path, paras, options=vrt_options)
    vrt_ds = None  # 保存并关闭

    # 输出 GeoTIFF + 压缩
    gdal.Translate(
        out_fdr_tif,
        vrt_path,
        format="GTiff",
        creationOptions=[
            "COMPRESS=LZW",
            "TILED=YES",
            "BIGTIFF=IF_SAFER"
        ]
    )
    os.remove(vrt_path)
    # gdal.Warp(outfile, files, options=[
    #     "-co", "COMPRESS=LZW",  # 选择LZW压缩
    #     "-co", "BIGTIFF=YES"])  # 允许大文件])
    print("合并完成，输出文件：", out_fdr_tif)

def huan_based_on_dem(sink_file,dem_file,out_fdr):
    """
    根据dem找到边界上的最低高程点作为流出点
    :param sink_file:
    :param dem_file:
    :return:
    """

    sink = Raster.get_raster(sink_file)
    proj,geo,s_nodata = Raster.get_proj_geo_nodata(sink_file)
    sink[sink!=1] = 0

    row,col = sink.shape

    dem = Raster.get_raster(dem_file)
    proj,geo,d_nodata = Raster.get_proj_geo_nodata(dem_file)

    # 找到最低点
    min_cell = [0,0,10000]
    for i in range(row):
        for j in range(col):
            if sink[i,j] == 1:
                continue
            flag = False
            for k in range(8):
                next_cell = (i+dmove[k][0],j+dmove[k][1])
                if not check_boundary(next_cell[0],next_cell[1],row,col):
                    continue
                if sink[next_cell[0],next_cell[1]] == 1:
                    flag = True
                    break

            if flag:
                if dem[i,j] < min_cell[2]:
                    min_cell = [i,j,dem[i,j]]

    print(min_cell)
    # for k in range(8):
    #     next_cell = (min_cell[0] + dmove[k][0], min_cell[1] + dmove[k][1])
    #     if sink[next_cell[0],next_cell[1]] == 1:
    #         min_cell = [next_cell[0],next_cell[1]]
    #         break


    print(valid.pixel_to_lonlat(geo,min_cell[0],min_cell[1]))
    maskDEM = dem.copy()
    maskDEM[sink != 1] = d_nodata
    A = heap_PF_D8.optimized_flow_repair(maskDEM,sink,0,(min_cell[0],min_cell[1]))   # 堆栈优化后的版本
    A[sink != 1] = 255
    Raster.save_raster(out_fdr,A,proj,geo,gdal.GDT_Byte,255)

def sbatch_huan_based_on_dem(venu):
    for file in os.listdir(venu):
        if len(file.split('.')) > 1:
            continue

        sink_file = os.path.join(venu, file, 'sink_' + str(file) + '_.tif')
        dem_file = os.path.join(venu, file, 'dem_' + str(file) + '_.tif')
        out_fdr = os.path.join(venu, file, 'modified_dir_' + str(file) + '_.tif')
        if os.path.exists(os.path.join(venu, file, 'modified_dir_' + str(file) + '_.tif')):
            continue

        # huan_based_on_dem(sink_file,dem_file,out_fdr)
        A = Raster.get_raster(os.path.join(venu, file, 'dir_' + str(file) + '_.tif'))
        proj,geo,f_nodata = Raster.get_proj_geo_nodata(os.path.join(venu, file, 'dir_' + str(file) + '_.tif'))

        A[A!=f_nodata] = 2
        A[A == f_nodata] = 255
        Raster.save_raster(out_fdr,A,proj,geo,gdal.GDT_Byte,255)



if __name__ == "__main__":

    import special_process

    huan_based_on_dem(r'F:\青藏高原水体数据集\New_DATA\全局环\EU\177\sink_177_.tif',
                      r'F:\青藏高原水体数据集\New_DATA\全局环\EU\177\dem_177_.tif',
                      r'F:\青藏高原水体数据集\New_DATA\全局环\EU\177\modified_fdr.tif')

    # file = r'F:\青藏高原水体数据集\3\EU\3\outdir_3_.tif'
    # A = Raster.get_raster(file)
    # proj,geo,nodata = Raster.get_proj_geo_nodata(file)
    # lon = 49.4902408#-58.9836726
    # lat = 53.4588929#-3.2430080
    #
    # row,col = valid.lonlat_to_pixel(geo,lon,lat)
    # print(A[row,col])
    # A[row,col] = 2
    #
    # A[A == 0] = nodata
    # Raster.save_raster(r'F:\青藏高原水体数据集\3\EU\3\outdir_3_new.tif',A,proj,geo,gdal.GDT_Byte,nodata)
    pass