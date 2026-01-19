# -*- coding: utf-8 -*-
"""
@Time ： 2025/3/11 11:28
@Auth ：
@File ：split_cal_acc.py
@IDE ：PyCharm
"""
import csv
import os
import pickle
import time
from scipy import ndimage
import numpy as np
from whitebox import WhiteboxTools
from osgeo import gdal
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
import multiprocessing

import postprocess
import sink
from genral_functions import *
import Raster

from multiprocessing import Pool
import numpy as np
import os

def Check():

    Po = Pool(10)
    fileNames = os.listdir("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/region/")
    paras = []
    for fileName in fileNames:
        paras.append(os.path.join("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/region/",fileName))

    for para in paras:
        Po.apply_async(postprocess.counterflow_process_check,(para,))

    Po.close()
    Po.join()


def split_raster(input_path, output_dir, block_height=2000):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        # os.chmod(output_dir,0o777)
    block_width = block_height
    fdir = Raster.get_raster(input_path)
    proj,geo,f_nodata = Raster.get_proj_geo_nodata(input_path)

    row,col = fdir.shape
    row_num = (row + block_height - 1) // block_height
    col_num = (col + block_width - 1) // block_width

    block_num = row_num * col_num

    for i in range(row_num):
        for j in range(col_num):

            single_id = i * col_num + j

            start_row = i * block_height
            end_row = start_row + block_height
            start_col = j * block_width
            end_col = start_col + block_width

            temp_arr = fdir[start_row:end_row,start_col:end_col].copy()

            if np.all(temp_arr== f_nodata) :
                print(row_num,col_num,single_id)
                continue

            new_geo = (geo[0] + geo[1] * start_col, geo[1], geo[2], geo[3] + geo[5] * start_row, geo[4], geo[5])

            out_file = os.path.join(output_dir,"fdir_" + str(row_num) + "_" + str(col_num) + "_" + str(single_id) + "_.tif")
            Raster.save_raster(out_file,temp_arr,proj,new_geo,gdal.GDT_Byte,f_nodata)

    print("共生成了{:d}块流向单元".format(block_num))

def sbatch_wht_acc(fdir_dir,output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        os.chmod(output_dir,0o777)
    paras = [[os.path.join(fdir_dir, fileName), " ", os.path.join(output_dir, fileName.split('.')[0] + '.tif')] for fileName in os.listdir(fdir_dir)]
    Po = Pool(15)
    for para in paras:
        # print(para)
        Po.apply_async(sink.Cal_acc, (para[0], para[2],))
        pass
    Po.close()
    Po.join()
def construct_tree_table(fdir_file,output_file,out_tree_file,out_table_file,height = 2000):
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
    # print(baseName)

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


    for r, c in zip(rows, cols):
        huan_vis = np.zeros((row, col))  # 用来检测环
        # print(f"边界像元：行 {r}, 列 {c}")
        r = int(r)
        c = int(c)
        if vis[r,c] != -9999:
            tree_table.setdefault((baseName,r,c), (baseName,id_cell[vis[r,c]][0],id_cell[vis[r,c]][1]))
            table_cell.setdefault((baseName, id_cell[vis[r,c]][0],id_cell[vis[r,c]][1]), table_cell_id[vis[r,c]])
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
                id_cell.setdefault(road_id,pop_cell)
                tree_table.setdefault((baseName,r, c),(baseName,pop_cell[0],pop_cell[1]))
                table_cell.setdefault((baseName, pop_cell[0],pop_cell[1]), (baseName, -1, -1)) # 记录非表间连接关系的序列，如流向内流或海洋
                table_cell_id.setdefault(road_id, (baseName, -1, -1))
                continue
            next_cell = (pop_cell[0] + dmove_dic[now_fdir][0], pop_cell[1] + dmove_dic[now_fdir][1])
            if not (0 <= next_cell[0] < row and 0 <= next_cell[1] < col):
                # 超出边界时的计算，要计算相邻的区域id
                id_cell.setdefault(road_id, pop_cell)
                tree_table.setdefault((baseName, r, c), (baseName, pop_cell[0], pop_cell[1]))
                new_row, new_col, next_id = cal_next_code(pop_cell, now_fdir, baseName, row, col)
                # if next_id != -9999:
                next_file = infos1[0] + '_' + str(new_row) + '_' + str(new_col) + '_' + str(next_id) + '_.tif'
                # next_file = infos1[0] + '_' + infos1[1] + '_' + infos1[2] + '_' + str(next_id) + '_.tif'
                table_cell.setdefault((baseName, pop_cell[0], pop_cell[1]), (next_file, (next_cell[0] + height) % height,(next_cell[1] + height) % height))  # (next_cell[0]+height)%height,(next_cell[1]+height)%height)
                table_cell_id.setdefault(road_id, (next_file, (next_cell[0] + height) % height, (next_cell[1] + height) % height))
                case = 1
                continue

            if vis[next_cell[0],next_cell[1]] != -9999:
                # 中止寻找下游，赋予前路栅格该id
                temp_id = vis[next_cell[0],next_cell[1]]
                for temp_cell in road_cells:
                    vis[temp_cell[0],temp_cell[1]] = temp_id
                tree_table.setdefault((baseName,r, c),(baseName,id_cell[temp_id][0],id_cell[temp_id][1]))
                table_cell.setdefault((baseName, id_cell[temp_id][0],id_cell[temp_id][1]), table_cell_id[temp_id])
                continue

            if fdir[next_cell[0],next_cell[1]] == f_nodata:
                # 没超出边界，但是nextCEll是空值，则记录当前file，和其他的信息
                id_cell.setdefault(road_id, pop_cell)
                tree_table.setdefault((baseName, r, c), (baseName, pop_cell[0], pop_cell[1]))

                # new_row, new_col, next_id = cal_next_code(pop_cell, now_fdir, baseName, row, col)

                # if next_id != -9999:
                next_file = infos1[0] + '_' + infos1[1] + '_' + infos1[2] + '_' + infos1[3] + '_.tif'
                # next_file = infos1[0] + '_' + infos1[1] + '_' + infos1[2] + '_' + str(next_id) + '_.tif'
                table_cell.setdefault((baseName, pop_cell[0], pop_cell[1]), (next_file, (next_cell[0] + height) % height,(next_cell[1] + height) % height))  # (next_cell[0]+height)%height,(next_cell[1]+height)%height)
                table_cell_id.setdefault(road_id, (next_file, (next_cell[0] + height) % height, (next_cell[1] + height) % height))
                case = 1
                continue
            if huan_vis[next_cell[0],next_cell[1]] == 1:
                print(baseName,'有环（内流区）',"需要修正流向栅格坐标为",next_cell,geo)
                break
            pop_cells.append(next_cell)
            road_cells.append(next_cell)
        if case == 1:
            for temp_cell in road_cells:
                vis[temp_cell[0], temp_cell[1]] = road_id
            road_id += 1


    # tree_table = {}  # 记录表内汇流
    # id_cell = {}
    #
    # table_cell_id = {}
    # table_cell = {} # 记录表间连接关系
    #
    # vis = np.zeros((row,col))
    # vis[:,:] = -9999
    # road_id = 1
    # # 上边界
    # for j in range(col):
    #     if fdir[0,j] == f_nodata:
    #         continue
    #     if vis[0,j] != -9999:
    #         tree_table.setdefault((baseName,0,j), (baseName,id_cell[vis[0,j]][0],id_cell[vis[0,j]][1]))
    #         table_cell.setdefault((baseName, id_cell[vis[0,j]][0],id_cell[vis[0,j]][1]), table_cell_id[vis[0,j]])
    #         continue
    #     pop_cells = [(0,j)]
    #     road_cells = [(0,j)]
    #     case = 0
    #     while pop_cells:
    #         pop_cell = pop_cells.pop()
    #         now_fdir = fdir[pop_cell[0],pop_cell[1]]
    #         if now_fdir not in dmove_dic:
    #             case = 1
    #             id_cell.setdefault(road_id,pop_cell)
    #             tree_table.setdefault((baseName,0,j),(baseName,pop_cell[0],pop_cell[1]))
    #             table_cell.setdefault((baseName, pop_cell[0],pop_cell[1]), (baseName, -1, -1)) # 记录非表间连接关系的序列，如流向内流或海洋
    #             table_cell_id.setdefault(road_id, (baseName, -1, -1))
    #             continue
    #         next_cell = (pop_cell[0] + dmove_dic[now_fdir][0], pop_cell[1]+dmove_dic[now_fdir][1])
    #         if not (0<=next_cell[0]<row and 0<=next_cell[1]<col):
    #             id_cell.setdefault(road_id, pop_cell)
    #             tree_table.setdefault((baseName,0,j),(baseName,pop_cell[0],pop_cell[1]))
    #             new_row,new_col,next_id = cal_next_code(pop_cell,now_fdir,baseName,row,col)
    #             # if next_id != -9999:
    #             next_file = infos1[0] + '_' + str(new_row) + '_' + str(new_col) + '_' + str(next_id) + '_.tif'
    #             # next_file = infos1[0] + '_' + infos1[1] + '_' + infos1[2] + '_' + str(next_id) + '_.tif'
    #             table_cell.setdefault((baseName,pop_cell[0],pop_cell[1]),(next_file,(next_cell[0]+height)%height,(next_cell[1]+height)%height))  #(next_cell[0]+height)%height,(next_cell[1]+height)%height)
    #             table_cell_id.setdefault(road_id,(next_file,(next_cell[0]+height)%height,(next_cell[1]+height)%height))
    #             case = 1
    #             continue
    #         if vis[next_cell[0],next_cell[1]] != -9999:
    #             # 中止寻找下游，赋予前路栅格该id
    #             temp_id = vis[next_cell[0],next_cell[1]]
    #             for temp_cell in road_cells:
    #                 vis[temp_cell[0],temp_cell[1]] = temp_id
    #             tree_table.setdefault((baseName,0,j),(baseName,id_cell[temp_id][0],id_cell[temp_id][1]))
    #             table_cell.setdefault((baseName, id_cell[temp_id][0],id_cell[temp_id][1]), table_cell_id[temp_id])
    #             continue
    #
    #         pop_cells.append(next_cell)
    #         road_cells.append(next_cell)
    #
    #     if case == 1:
    #         for temp_cell in road_cells:
    #             vis[temp_cell[0], temp_cell[1]] = road_id
    #         road_id += 1
    # # 下边界
    # for j in range(col):
    #     if fdir[row-1, j] == f_nodata:
    #         continue
    #     if vis[row-1, j] != -9999:
    #         tree_table.setdefault((baseName,row-1, j), (baseName,id_cell[vis[row-1, j]][0],id_cell[vis[row-1, j]][1]))
    #         table_cell.setdefault((baseName, id_cell[vis[row-1, j]][0],id_cell[vis[row-1, j]][1]), table_cell_id[vis[row-1, j]])
    #         continue
    #     pop_cells = [(row-1, j)]
    #     road_cells = [(row-1, j)]
    #     case = 0
    #     while pop_cells:
    #         pop_cell = pop_cells.pop()
    #         now_fdir = fdir[pop_cell[0], pop_cell[1]]
    #         if now_fdir not in dmove_dic:
    #             case = 1
    #             id_cell.setdefault(road_id, pop_cell)
    #             tree_table.setdefault((baseName,row-1, j), (baseName,pop_cell[0],pop_cell[1]))
    #             table_cell.setdefault((baseName, pop_cell[0],pop_cell[1]), (baseName, -1, -1))  # 记录非表间连接关系的序列，如流向内流或海洋
    #             table_cell_id.setdefault(road_id, (baseName, -1, -1))
    #             continue
    #         next_cell = (pop_cell[0] + dmove_dic[now_fdir][0], pop_cell[1] + dmove_dic[now_fdir][1])
    #         if not (0 <= next_cell[0] < row and 0 <= next_cell[1] < col):
    #             case = 1
    #             id_cell.setdefault(road_id, pop_cell)
    #             tree_table.setdefault((baseName,row - 1, j), (baseName,pop_cell[0],pop_cell[1]))
    #             new_row,new_col,next_id = cal_next_code(pop_cell, now_fdir, baseName, row, col)
    #             # if next_id!=-9999:
    #             next_file = infos1[0] + '_' + str(new_row) + '_' + str(new_col) + '_' + str(next_id) + '_.tif'
    #             # next_file = infos1[0] + '_' + infos1[1] + '_' + infos1[2] + '_' + str(next_id) + '_.tif'
    #             table_cell.setdefault((baseName, pop_cell[0],pop_cell[1]), (next_file,(next_cell[0]+height)%height,(next_cell[1]+height)%height))
    #             table_cell_id.setdefault(road_id, (next_file,(next_cell[0]+height)%height,(next_cell[1]+height)%height))
    #             continue
    #         if vis[next_cell[0], next_cell[1]] != -9999:
    #             # 中止寻找下游，赋予前路栅格该id
    #             temp_id = vis[next_cell[0], next_cell[1]]
    #             for temp_cell in road_cells:
    #                 vis[temp_cell[0], temp_cell[1]] = temp_id
    #             tree_table.setdefault((baseName,row - 1, j), (baseName, id_cell[temp_id][0],id_cell[temp_id][1]))
    #             table_cell.setdefault((baseName, id_cell[temp_id][0],id_cell[temp_id][1]), table_cell_id[temp_id])
    #             continue
    #
    #         pop_cells.append(next_cell)
    #         road_cells.append(next_cell)
    #
    #     if case == 1:
    #         for temp_cell in road_cells:
    #             vis[temp_cell[0], temp_cell[1]] = road_id
    #         road_id += 1
    # # 左边界
    # for i in range(row):
    #     if fdir[i,0] == f_nodata:
    #         continue
    #     if vis[i,0] != -9999:
    #         tree_table.setdefault((baseName,i,0), (baseName,id_cell[vis[i,0]][0],id_cell[vis[i,0]][1]))
    #         table_cell.setdefault((baseName,id_cell[vis[i,0]][0],id_cell[vis[i,0]][1]), table_cell_id[vis[i,0]])
    #         continue
    #     pop_cells = [(i,0)]
    #     road_cells = [(i,0)]
    #     case = 0
    #     while pop_cells:
    #         pop_cell = pop_cells.pop()
    #         now_fdir = fdir[pop_cell[0], pop_cell[1]]
    #         if now_fdir not in dmove_dic:
    #             case = 1
    #             id_cell.setdefault(road_id, pop_cell)
    #             tree_table.setdefault((baseName,i,0), (baseName, pop_cell[0],pop_cell[1]))
    #             table_cell.setdefault((baseName, pop_cell[0],pop_cell[1]), (baseName, -1, -1))  # 记录非表间连接关系的序列，如流向内流或海洋
    #             table_cell_id.setdefault(road_id, (baseName, -1, -1))
    #             continue
    #         next_cell = (pop_cell[0] + dmove_dic[now_fdir][0], pop_cell[1] + dmove_dic[now_fdir][1])
    #         if not (0 <= next_cell[0] < row and 0 <= next_cell[1] < col):
    #             case = 1
    #             id_cell.setdefault(road_id, pop_cell)
    #             tree_table.setdefault((baseName,i, 0), (baseName, pop_cell[0],pop_cell[1]))
    #             new_row,new_col,next_id = cal_next_code(pop_cell, now_fdir, baseName, row, col)
    #             # if next_id != -9999:
    #             next_file = infos1[0] + '_' + str(new_row) + '_' + str(new_col) + '_' + str(next_id) + '_.tif'
    #             # next_file = infos1[0] + '_' + infos1[1] + '_' + infos1[2] + '_' + str(next_id) + '_.tif'
    #             table_cell.setdefault((baseName, pop_cell[0],pop_cell[1]), (next_file,(next_cell[0]+height)%height,(next_cell[1]+height)%height))
    #             table_cell_id.setdefault(road_id, (next_file,(next_cell[0]+height)%height,(next_cell[1]+height)%height))
    #             continue
    #         if vis[next_cell[0], next_cell[1]] != -9999:
    #             # 中止寻找下游，赋予前路栅格该id
    #             temp_id = vis[next_cell[0], next_cell[1]]
    #             for temp_cell in road_cells:
    #                 vis[temp_cell[0], temp_cell[1]] = temp_id
    #             tree_table.setdefault((baseName,i, 0), (baseName, id_cell[temp_id][0],id_cell[temp_id][1]))
    #             table_cell.setdefault((baseName, id_cell[temp_id][0],id_cell[temp_id][1]), table_cell_id[temp_id])
    #             continue
    #
    #         pop_cells.append(next_cell)
    #         road_cells.append(next_cell)
    #
    #     if case == 1:
    #         for temp_cell in road_cells:
    #             vis[temp_cell[0], temp_cell[1]] = road_id
    #         road_id += 1
    # # 右边界
    # for i in range(row):
    #     if fdir[i,col-1] == f_nodata:
    #         continue
    #     if vis[i,col-1] != -9999:
    #         tree_table.setdefault((baseName,i, col - 1), (baseName, id_cell[vis[i,col-1]][0],id_cell[vis[i,col-1]][1]))
    #         table_cell.setdefault((baseName, id_cell[vis[i,col-1]][0],id_cell[vis[i,col-1]][1]), table_cell_id[vis[i, col - 1]])
    #         continue
    #     pop_cells = [(i,col-1)]
    #     road_cells = [(i,col-1)]
    #     case = 0
    #     while pop_cells:
    #         pop_cell = pop_cells.pop()
    #         now_fdir = fdir[pop_cell[0], pop_cell[1]]
    #         if now_fdir not in dmove_dic:
    #             case = 1
    #             id_cell.setdefault(road_id, pop_cell)
    #             tree_table.setdefault((baseName,i,col-1), (baseName, pop_cell[0],pop_cell[1]))
    #             table_cell.setdefault((baseName, pop_cell[0],pop_cell[1]), (baseName, -1, -1))  # 记录非表间连接关系的序列，如流向内流或海洋
    #             table_cell_id.setdefault(road_id, (baseName, -1, -1))
    #             continue
    #         next_cell = (pop_cell[0] + dmove_dic[now_fdir][0], pop_cell[1] + dmove_dic[now_fdir][1])
    #         if not (0 <= next_cell[0] < row and 0 <= next_cell[1] < col):
    #             case = 1
    #             id_cell.setdefault(road_id, pop_cell)
    #             tree_table.setdefault((baseName,i, col - 1), (baseName, pop_cell[0],pop_cell[1]))
    #             new_row,new_col,next_id = cal_next_code(pop_cell, now_fdir, baseName, row, col)
    #             # if next_id != -9999:
    #             next_file = infos1[0] + '_' + str(new_row) + '_' + str(new_col) + '_' + str(next_id) + '_.tif'
    #             # next_file = infos1[0] + '_' + infos1[1] + '_' + infos1[2] + '_' + str(next_id) + '_.tif'
    #             table_cell.setdefault((baseName, pop_cell[0],pop_cell[1]), (next_file,(next_cell[0]+height)%height,(next_cell[1]+height)%height))
    #             table_cell_id.setdefault(road_id, (next_file,(next_cell[0]+height)%height,(next_cell[1]+height)%height))
    #             continue
    #         if vis[next_cell[0], next_cell[1]] != -9999:
    #             # 中止寻找下游，赋予前路栅格该id
    #             temp_id = vis[next_cell[0], next_cell[1]]
    #             for temp_cell in road_cells:
    #                 vis[temp_cell[0], temp_cell[1]] = temp_id
    #             tree_table.setdefault((baseName,i, col - 1),(baseName, id_cell[temp_id][0],id_cell[temp_id][1]))
    #             table_cell.setdefault((baseName, id_cell[temp_id][0],id_cell[temp_id][1]), table_cell_id[temp_id])
    #             continue
    #
    #         pop_cells.append(next_cell)
    #         road_cells.append(next_cell)
    #
    #     if case == 1:
    #         for temp_cell in road_cells:
    #             vis[temp_cell[0], temp_cell[1]] = road_id
    #         road_id += 1
    # print(table_cell)
    # Raster.save_raster(output_file,vis,proj,geo,gdal.GDT_Float32,-9999)

    tree_table_csv = list(tree_table.items())
    with open(out_tree_file,'w',newline="") as f:
        writer = csv.writer(f)
        writer.writerows(tree_table_csv)
        f.close()
    os.chmod(out_tree_file,0o777)
    # for i in table_cell:
    #     print(i,table_cell[i])
    table_cell_csv = list(table_cell.items())
    with open(out_table_file,'w',newline="") as f:
        writer = csv.writer(f)
        writer.writerows(table_cell_csv)
        f.close()
    os.chmod(out_tree_file,0o777)
def construct_tree_masks(fdir_dir,output_dir,outtable_dir,height = 2000):


    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        # os.chmod(output_dir,0o777)
    if not os.path.exists(outtable_dir):
        os.makedirs(outtable_dir, exist_ok=True)
        # os.chmod(outtable_dir,0o777)

    # paras = [[os.path.join(fdir_dir,fileName)," ",os.path.join(output_dir,fileName.split('.')[0]+'tree_.csv'),os.path.join(outtable_dir,fileName.split('.')[0]+'table_.csv')] for fileName in os.listdir(fdir_dir)]
    paras = []
    for fileName in os.listdir(fdir_dir):
        if fileName.split('.')[-1]!='tif':
            continue
        paras.append([os.path.join(fdir_dir,fileName),os.path.join(output_dir,fileName.split('.')[0]+'tree_.tif'),os.path.join(output_dir,fileName.split('.')[0]+'tree_.csv'),os.path.join(outtable_dir,fileName.split('.')[0]+'table_.csv')])
    # print(paras)
    Po = Pool(5)
    for para in paras:
        # print(para)
        Po.apply_async(construct_tree_table,(para[0],para[1],para[2],para[3],height,))
        # construct_tree_table(para[0],para[1],para[2],para[3],height)
    Po.close()
    Po.join()
    # construct_tree_mask(r'F:\FABDEM\split_trail\region\fdir_4_4_0_.tif',r'F:\FABDEM\split_trail\region_mask\A.tif')

def merge_tree_table(tree_dir,table_dir,merge_tree,merge_table):
    """
    合并生成的表内汇流关系表和表间关系表。
    表间关系表用来查询更新的汇流累积量。
    :param output_dir:
    :param outtable_dir:
    :param merge_tree:
    :param merge_table:
    :return:
    """
    fileNames = [os.path.join(tree_dir,fileName) for fileName in os.listdir(tree_dir)]
    merge_tree_csv = []
    for fileName in fileNames:
        with open(fileName,'r') as f:
            reader = csv.reader(f)
            for i in reader:
                # print(i)
                merge_tree_csv.append(i)
            f.close()
    with open(merge_tree,'w',newline="") as f:
        writer = csv.writer(f)
        writer.writerows(merge_tree_csv)
        f.close()

    fileNames = [os.path.join(table_dir, fileName) for fileName in os.listdir(table_dir)]
    merge_tree_csv = []
    for fileName in fileNames:
        with open(fileName, 'r') as f:
            reader = csv.reader(f)
            for i in reader:
                # print(i[1].split(',')[0].split('_')[-2] == '-9999')
                #
                # if i[1].split(',')[0].split('_')[-2] == '-9999':
                #     continue
                merge_tree_csv.append(i)
            f.close()

    with open(merge_table, 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerows(merge_tree_csv)
        f.close()

def get_upstream(node,table,tree,results):
    """
    根据节点追溯上游
    :param node:
    :param table:
    :param tree:
    :return:
    """
    upstream = set()
    # print(table[start])
    pop_cells = [node]
    while pop_cells:
        pop_cell = pop_cells.pop()
        # if pop_cell in results:
        #     upstream += results[pop_cell]
        #     break
        # print(pop_cell)
        if pop_cell not in table:
            continue
        upstreams = table[pop_cell]
        # upstream += table[pop_cell]
        for temp_upstream in upstreams:
            pop_cells += tree[temp_upstream]
            upstream.add(temp_upstream)
            # upstream += tree[temp_upstream]
            # upstream.append(temp_upstream)  # 要加上直接连接的栅格，还有这些栅格的上游

    return node,list(upstream)
def cal_conflunce(tree_file,table_file,out_acc_file):
    """
    构建汇流树、计算每块栅格应该更新的汇流累积量
    :param tree_file:
    :param table_file:
    :return:
    """
    tree = {} # 表内关系  下游:[上游]
    table = {} # 表间关系  下游:[上游]
    with open(tree_file,'r') as f:
        reader = csv.reader(f)
        for i in reader:
            tree.setdefault(i[1],[]).append(i[0])
        f.close()

    with open(table_file,'r') as f:
        reader = csv.reader(f)
        for i in reader:
            table.setdefault(i[1],[]).append(i[0])
        f.close()
    result = {}
    acc_result = []
    for node in table:
        node,upstreams = get_upstream(node,table,tree,result)
        result.setdefault(node,upstreams)
        # print(node,result[node])
        acc_result.append([[node]+upstreams])


    with open(out_acc_file,'wb') as  f:
        pickle.dump(result,f)
        f.close()

def make_acc_table(acc_dir,out_file):
    """
    将每块单元的边界像元的汇流值记录下来，{(文件，行，列)：acc}
    :param acc_dir:
    :param out_file:
    :return:
    """


    file_acc = {}
    for fileName in os.listdir(acc_dir):
        if fileName.split('.')[-1] != 'tif':
            continue
        filePath = os.path.join(acc_dir,fileName)
        acc = Raster.get_raster(filePath)
        proj, geo, a_nodata = Raster.get_proj_geo_nodata(filePath)
        row, col = acc.shape

        mask = np.zeros((row, col), dtype=bool)
        mask[acc != a_nodata] = True
        # 记录边界shange
        # 找边缘：腐蚀后与原始掩膜做异或（边界为 True）
        eroded = ndimage.binary_erosion(mask)
        boundary = mask ^ eroded  # XOR 得到边界像元

        # 获取边界像元的行列号
        rows, cols = np.where(boundary)
        # 输出
        for i, j in zip(rows,cols):
            file_acc.setdefault((fileName, i, j), acc[i, j])

        # for j in range(col):
        #     # 第一行
        #     if acc[0,j] != a_nodata:
        #         file_acc.setdefault((fileName,0,j),acc[0,j])
        # for j in range(col):
        #     # 最后一行
        #     if acc[row-1,j] != a_nodata:
        #         file_acc.setdefault((fileName,row-1,j),acc[row-1,j])
        # for i in range(row):
        #     # 第一列
        #     if acc[i,0] != a_nodata:
        #         file_acc.setdefault((fileName,i,0),acc[i,0])
        # for i in range(row):
        #     # 最后一列
        #     if acc[i,col-1] != a_nodata:
        #         file_acc.setdefault((fileName,i,col-1),acc[i,col-1])

    with open(out_file,'wb') as f:
        pickle.dump(file_acc,f)
        f.close()
def cal_acc_update(acc_upstream_file,file_acc_file,out_acc_update_dir):
    """
    计算每个传递点应该更新的上游汇流累积量
    :param acc_upstream_file:
    :param acc_dir:
    :param out_acc_update_file:
    :return:
    """
    if not os.path.exists(out_acc_update_dir):
        os.mkdir(out_acc_update_dir)
        os.chmod(out_acc_update_dir,0o777)
    with open(file_acc_file,'rb') as f:
        file_acc = pickle.load(f)
        for i in file_acc:
            print(i)
            break
        f.close()

    with open(acc_upstream_file,'rb') as f:
        cell_upstream = pickle.load(f)
        for i in cell_upstream:
            print(i)
            break
        f.close()
    # print(cell_upstream)

    result = {}
    for cell in cell_upstream:
        upstreams = cell_upstream[cell]


        acc_num = 0
        for upstream in upstreams:
            infos = upstream.split(',')
            fileName = infos[0][2:-1]
            # print(infos)
            # row = infos[1][1:]
            # col = infos[2][1:-1]
            row = int(float(infos[1][1:]))
            col = int(float(infos[2][1:-1]))

            # if fileName == 'fdir_4_4_7_.tif':
            #     print(cell, upstreams)
            #     break
            if (fileName, row, col) in file_acc:
                acc_num += float(file_acc[(fileName, row, col)]) + 1
        result.setdefault(cell.split(',')[0][2:-1],[]).append((int(float(cell.split(',')[1][1:])),int(float(cell.split(',')[2][1:-1])),acc_num))

    for i in result:
        # print(i)
        filpath = os.path.join(out_acc_update_dir,i[:-4]+'.pkl')
        with open(filpath,'wb') as f:
            pickle.dump(result[i],f)
            f.close()
    #
    # print(result)

def update_acc(acc_table_dir,acc_dir,fdir_dir,out_acc_dir):
    """
    更新汇流累积量
    :param acc_table_dir:
    :param acc_dir:
    :param out_acc_dir:
    :return:
    """
    if not os.path.exists(out_acc_dir):
        os.mkdir(out_acc_dir)
        os.chmod(out_acc_dir,0o777)
    # print(os.listdir(acc_table_dir))
    for fileName in os.listdir(acc_table_dir):

        acc_file = os.path.join(acc_dir,fileName[:-4]+'.tif')
        print(fileName,acc_file)

        # if fileName[:-4]+'.tif' != 'fdir_20_34_192_.tif':
        #     continue

        if not os.path.exists(acc_file):
            continue
        fdir_file = os.path.join(fdir_dir,fileName[:-4]+'.tif')
        acc_table_path = os.path.join(acc_table_dir, fileName)
        fdir = Raster.get_raster(fdir_file)
        acc = Raster.get_raster(acc_file)
        CCC = acc.copy().astype(np.float64)


        proj,geo,a_nodata = Raster.get_proj_geo_nodata(acc_file)
        row,col = acc.shape

        # vis = np.zeros((10000,10000))

        with open(acc_table_path,'rb') as f:
            con = pickle.load(f)
            for i in con:
                # print(int(float(i[0])),int(float(i[1])),float(i[2]))
                # print(i)
                update_acc = float(i[2])
                x = int(float(i[0]))
                y = int(float(i[1]))
                # 强制修正边界。因为计算边界传递时是按照固定值height计算的，所以接到的像元可能并不是刚好的height，所以需要匹配修正一下
                if int(float(i[0])) >= row:
                    x = row - 1
                if int(float(i[1])) >= col:
                    y = col -1
                if int(float(i[0])) < 0:
                    x = 0
                if int(float(i[1])) < 0:
                    y = 0
                pop_cells = [[x,y]]

                vis = np.zeros((10000, 10000))   # 这一行改过，原本在L829,改对了

                while pop_cells:
                    pop_cell = pop_cells.pop()
                    vis[pop_cell[0],pop_cell[1]] = 1

                    CCC[pop_cell[0],pop_cell[1]] += float(update_acc)

                    now_fdir = fdir[pop_cell[0], pop_cell[1]]
                    if now_fdir not in dmove_dic:
                        continue
                    next_cell = [pop_cell[0] + dmove_dic[now_fdir][0], pop_cell[1] + dmove_dic[now_fdir][1]]
                    if not (0 <= next_cell[0] < row and 0 <= next_cell[1] < col):
                        continue
                    if vis[next_cell[0],next_cell[1]] == 1:
                        print(pop_cell,now_fdir,next_cell,fdir[next_cell[0],next_cell[1]])
                        continue
                    pop_cells.append(next_cell)
            f.close()

        out_new_acc = os.path.join(out_acc_dir,fileName[:-4]+'acc.tif')
        Raster.save_raster(out_new_acc,CCC,proj,geo,gdal.GDT_Float64,a_nodata)
        # break

def merge_update_acc(new_acc_dir,out_final_acc):
    """
    合并最终更新的acc。
    :param new_acc_dir:
    :param out_final_acc:
    :return:
    """
    files = [os.path.join(new_acc_dir,fileName) for fileName in os.listdir(new_acc_dir)]
    # files = ["/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final1/special/PF/special_outdir.tif","/datanode05/zhangbin/FAB_hydrography/Result/SouthAmerica_ACC_final1/region/SouthAmerica_fdr.tif"]
    proj,geo,nodata = Raster.get_proj_geo_nodata(files[0])
    # 构造 BuildVRTOptions（正确方式）
    vrt_options = gdal.BuildVRTOptions(
        separate=False,
        resolution='highest',
        VRTNodata=nodata  # 设定 NoData 值，可根据你的影像情况修改
    )
    baseDir = os.path.dirname(out_final_acc)
    # 构建 VRT
    vrt_path = os.path.join(baseDir, "temp_merged3.vrt")
    vrt_ds = gdal.BuildVRT(vrt_path, files, options=vrt_options)
    vrt_ds = None  # 保存并关闭


    # 输出 GeoTIFF + 压缩
    gdal.Translate(
        out_final_acc,
        vrt_path,
        format="GTiff",
        creationOptions=[
            "COMPRESS=LZW",
            "TILED=YES",
            "BIGTIFF=IF_SAFER"
        ]
    )
    os.chmod(vrt_path,0o777)
    os.remove(vrt_path)
    A = Raster.get_raster(out_final_acc)
    proj,geo,nodata = Raster.get_proj_geo_nodata(out_final_acc)
    A[A==nodata] = -9999
    Raster.save_raster(out_final_acc,A,proj,geo,gdal.GDT_Float64,-9999)

def cal_acc_work_flow_example():
    """
    测试数据及代码.
    note：尚未加入自动计算汇流累积量的代码，需要自己手动准备acc目录
    :return:
    """
    fdir_file = r'F:\FABDEM\split_trail\测试数据\DEM\run_fdir.tif'
    region = r'F:\FABDEM\split_trail\测试数据\region'  # 存放分块流向的目录
    acc = r'F:\FABDEM\split_trail\测试数据\acc'  # 存放分块初始汇流累积量的目录
    region_tree = r'F:\FABDEM\split_trail\测试数据\region_tree'  # 存放分块内汇流的目录：汇入点-流出点
    region_table = r'F:\FABDEM\split_trail\测试数据\region_table'  # 存放分块间汇流的目录：
    merge_tree = r'F:\FABDEM\split_trail\测试数据\merge_tree.csv' #  记录总体的表内汇流
    merge_table = r'F:\FABDEM\split_trail\测试数据\merge_table.csv' # 记录总体的表间汇流
    acc_upstream = r'F:\FABDEM\split_trail\测试数据\acc_upstream.pkl' # 记录每个汇入点的上游
    file_acc = r'F:\FABDEM\split_trail\测试数据\file_acc_origin.pkl' # 记录每个交界点的初始汇流累积量
    acc_table = r'F:\FABDEM\split_trail\测试数据\acc_table' # 记录每分块的传输汇流值
    new_acc = r'F:\FABDEM\split_trail\测试数据\new_acc'  # 存放更新后分块汇流累积量的目录
    merge_acc = r'F:\FABDEM\split_trail\测试数据\FINAL_ACC.tif'  # 存放汇流累积量的文件

    t1= time.time()
    split_raster(fdir_file, region, block_height= 2000)
    t2 = time.time()
    print("Spliting region consumes {:f}s".format((t2-t1)))
    #
    construct_tree_masks(region, region_tree,region_table,height=2000)
    t3 = time.time()
    print("constructing region confluence consumes {:f}s".format((t3 - t2)))

    merge_tree_table(region_tree, region_table,merge_tree,merge_table)
    t4 = time.time()
    print("Merging confluence consumes {:f}s".format((t4 - t3)))

    cal_conflunce(merge_tree, merge_table,acc_upstream)
    t5 = time.time()
    print("Calculating confluence consumes {:f}s".format((t5 - t4)))

    make_acc_table(acc, file_acc)
    t6 = time.time()
    print("Making acc table consumes {:f}s".format((t6 - t5)))

    cal_acc_update(acc_upstream, file_acc,acc_table)
    t7 = time.time()
    print("Getting upstreams consumes {:f}s".format((t7 - t6)))

    update_acc(acc_table, acc, region,new_acc)
    t8 = time.time()
    print("Updating acc consumes {:f}s".format((t8 - t7)))
    #
    merge_update_acc(new_acc, merge_acc)
    t9 = time.time()
    print("Merge acc consumes {:f}s".format((t9 - t8)))

    print("All work flow consumes {:f}s".format((t9 - t1)))

def cal_acc_work_flow_SA():
    """
    测试数据及代码.
    note：尚未加入自动计算汇流累积量的代码，需要自己手动准备acc目录
    :return:
    """
    fdir_file = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/Final_FDIR.tif"
    region = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/region"  # 存放分块流向的目录
    acc = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/acc"  # 存放分块初始汇流累积量的目录
    region_tree = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/region_tree"  # 存放分块内汇流的目录：汇入点-流出点
    region_table = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/region_table"  # 存放分块间汇流的目录：
    merge_tree = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/merge_tree.csv" #  记录总体的表内汇流
    merge_table = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/merge_table.csv" # 记录总体的表间汇流
    acc_upstream = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/acc_upstream.pkl" # 记录每个汇入点的上游
    file_acc = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/file_acc_origin.pkl" # 记录每个交界点的初始汇流累积量
    acc_table = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/acc_table" # 记录每分块的传输汇流值
    new_acc = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/new_acc"  # 存放更新后分块汇流累积量的目录
    merge_acc = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/FINAL_ACC.tif"  # 存放汇流累积量的文件


    t1= time.time()
    # split_raster(fdir_file, region, block_height= 10000)
    # sbatch_wht_acc("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/region/","/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/acc/")
    t2 = time.time()
    print("Spliting region consumes {:f}s".format((t2-t1)))
    # # #
    construct_tree_masks(region, region_tree,region_table,height=10000)
    t3 = time.time()
    print("constructing region confluence consumes {:f}s".format((t3 - t2)))

    merge_tree_table(region_tree, region_table,merge_tree,merge_table)
    t4 = time.time()
    print("Merging confluence consumes {:f}s".format((t4 - t3)))

    cal_conflunce(merge_tree, merge_table,acc_upstream)
    t5 = time.time()
    print("Calculating confluence consumes {:f}s".format((t5 - t4)))

    make_acc_table(acc, file_acc)
    t6 = time.time()
    print("Making acc table consumes {:f}s".format((t6 - t5)))

    cal_acc_update(acc_upstream, file_acc,acc_table)
    t7 = time.time()
    print("Getting upstreams consumes {:f}s".format((t7 - t6)))
    #
    update_acc(acc_table, acc, region,new_acc)
    t8 = time.time()
    print("Updating acc consumes {:f}s".format((t8 - t7)))
    #
    merge_update_acc(new_acc, merge_acc)
    t9 = time.time()
    print("Merge acc consumes {:f}s".format((t9 - t8)))
    #
    print("All work flow consumes {:f}s".format((t9 - t1)))

def cal_normal_acc_work_flow(fdir_file,acc_venu):
    """
    测试数据及代码.
    note：尚未加入自动计算汇流累积量的代码，需要自己手动准备acc目录
    :return:
    """
    # fdir_file = fdir_file
    # region = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/region"  # 存放分块流向的目录
    # acc = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/acc"  # 存放分块初始汇流累积量的目录
    # region_tree = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/region_tree"  # 存放分块内汇流的目录：汇入点-流出点
    # region_table = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/region_table"  # 存放分块间汇流的目录：
    # merge_tree = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/merge_tree.csv" #  记录总体的表内汇流
    # merge_table = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/merge_table.csv" # 记录总体的表间汇流
    # acc_upstream = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/acc_upstream.pkl" # 记录每个汇入点的上游
    # file_acc = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/file_acc_origin.pkl" # 记录每个交界点的初始汇流累积量
    # acc_table = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/acc_table" # 记录每分块的传输汇流值
    # new_acc = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/new_acc"  # 存放更新后分块汇流累积量的目录
    # merge_acc = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/FINAL_ACC.tif"  # 存放汇流累积量的文件
    if not os.path.exists(acc_venu):
        os.mkdir(acc_venu)
    region = os.path.join(acc_venu,'region')
    acc = os.path.join(acc_venu,'acc')
    region_tree = os.path.join(acc_venu,'region_tree')
    region_table = os.path.join(acc_venu,'region_table')
    merge_tree = os.path.join(acc_venu,'merge_tree.csv')
    merge_table = os.path.join(acc_venu,'merge_table.csv')
    acc_upstream = os.path.join(acc_venu,'acc_upstream.pkl')
    file_acc = os.path.join(acc_venu,'file_acc_origin.pkl')
    acc_table = os.path.join(acc_venu,'acc_table')
    new_acc = os.path.join(acc_venu,'new_acc')
    merge_acc = os.path.join(acc_venu,'FINAL_ACC.tif')


    t1= time.time()
    split_raster(fdir_file, region, block_height= 10000)
    sbatch_wht_acc(region,acc)
    t2 = time.time()
    print("Spliting region consumes {:f}s".format((t2-t1)))
    # # #
    construct_tree_masks(region, region_tree,region_table,height=10000)
    t3 = time.time()
    print("constructing region confluence consumes {:f}s".format((t3 - t2)))

    merge_tree_table(region_tree, region_table,merge_tree,merge_table)
    t4 = time.time()
    print("Merging confluence consumes {:f}s".format((t4 - t3)))
    #
    cal_conflunce(merge_tree, merge_table,acc_upstream)
    t5 = time.time()
    print("Calculating confluence consumes {:f}s".format((t5 - t4)))
    #
    make_acc_table(acc, file_acc)
    t6 = time.time()
    print("Making acc table consumes {:f}s".format((t6 - t5)))

    cal_acc_update(acc_upstream, file_acc,acc_table)
    t7 = time.time()
    print("Getting upstreams consumes {:f}s".format((t7 - t6)))
    #
    update_acc(acc_table, acc, region,new_acc)
    t8 = time.time()
    print("Updating acc consumes {:f}s".format((t8 - t7)))
    #
    merge_update_acc(new_acc, merge_acc)
    t9 = time.time()
    print("Merge acc consumes {:f}s".format((t9 - t8)))
    #
    print("All work flow consumes {:f}s".format((t9 - t1)))

if __name__ == "__main__":

    # 测试
    # cal_acc_work_flow()

    # 应用SA
    # cal_acc_work_flow_SA()

    # construct_tree_table(r'F:\FABDEM\split_trail\region\fdir_4_4_0_.tif',r'',r'F:\FABDEM\split_trail\region\s.csv')

    # split_raster(r'F:\FABDEM\split_trail\fdir.tif',r'F:\FABDEM\split_trail\region')
    # construct_tree_masks(r'F:\FABDEM\split_trail\region',r'F:\FABDEM\split_trail\region_mask')
    # fdir_dir = r'F:\FABDEM\split_trail\region'
    # output_dir = r'F:\FABDEM\split_trail\region_mask'
    # paras = [[os.path.join(fdir_dir, fileName), os.path.join(output_dir, fileName.split('.')[0] + 'mask_.tif')] for
    #          fileName in os.listdir(fdir_dir)]
    # print(paras)
    # Po = Pool(6)
    # for para in paras:
    #     print(para)
    #     Po.apply_async(construct_tree_mask,(para[0],para[1]))
    #
    # Po.close()
    # Po.join()
    #
    # construct_tree_masks(r'F:\FABDEM\split_trail\region',r'F:\FABDEM\split_trail\region_tree',r'F:\FABDEM\split_trail\region_table')
    # merge_tree_table(r'F:\FABDEM\split_trail\region_tree',r'F:\FABDEM\split_trail\region_table',r'F:\FABDEM\split_trail\merge_tree.csv',r'F:\FABDEM\split_trail\merge_table.csv')
    # cal_conflunce(r'F:\FABDEM\split_trail\merge_tree.csv',r'F:\FABDEM\split_trail\merge_table.csv',r'F:\FABDEM\split_trail\acc_upstream.pkl')
    # make_acc_table(r'F:\FABDEM\split_trail\acc',r'F:\FABDEM\split_trail\file_acc_origin.pkl')
    # cal_acc_update(r'F:\FABDEM\split_trail\acc_upstream.pkl',r'F:\FABDEM\split_trail\file_acc_origin.pkl',r'F:\FABDEM\split_trail\acc_table')
    # update_acc(r'F:\FABDEM\split_trail\acc_table',r'F:\FABDEM\split_trail\acc',r'F:\FABDEM\split_trail\region',r'F:\FABDEM\split_trail\new_acc')
    # merge_update_acc(r'F:\FABDEM\split_trail\new_acc',r'F:\FABDEM\split_trail\updated_acc.tif')

    # A = Raster.get_raster( r'F:\FABDEM\split_trail\测试数据\FINAL_ACC.tif')
    # B = Raster.get_raster( r'F:\FABDEM\split_trail\测试数据\DEM\acc_arcmap.tif')
    # print(np.all(A==B))



    construct_tree_table(r'F:\青藏高原水体数据集\DATA\Acc\fdir_27_23_379_f.tif',
                         "",
                         r'F:\青藏高原水体数据集\DATA\Acc\fdir_27_23_379_tree.csv',
                         r'F:\青藏高原水体数据集\DATA\Acc\fdir_27_23_379_table.csv',height=10000)

    pass
