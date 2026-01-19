# -*- coding: utf-8 -*-
"""
计算GRDC提供的站点在FABDir中的上游，做面积误差计算和空间重叠验证
1）要计算两套站点的上游，分别存储在不同的文件夹，转矢量，投影计算面积
2)记录GRDC的提供的面积和矢量位置
3)验证

@Time ： 2025/3/8 12:52
@Auth ：
@File ：valid.py
@IDE ：PyCharm
"""
import csv
import math
import os
import pickle
import shutil

import numpy as np
from osgeo import gdal,ogr,osr
import pandas as pd

import Find
import Raster
from multiprocessing import Pool
from pyproj import Geod
import geopandas as gpd

import compare_NHD
import genral_functions
import sink
import matplotlib.pyplot as plt


dmove = [(0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1), (-1, 0), (-1, 1)]
dmove_dic = {1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1), 16: (0, -1), 32: (-1, -1), 64: (-1, 0), 128: (-1, 1)}

def get_rever_D8(dir, row, col):
    """
    查询输入栅格的上游栅格
    :param dir: array of dir
    :param row: row of the cell
    :param col:
    :return: [(i,j),(),]
    """
    up_cell = []
    row_num, col_num = dir.shape

    for i in range(8):
        now_loc = (row + dmove[i][0], col + dmove[i][1])
        # print(now_loc)
        if 0<=now_loc[0]<row_num and 0<=now_loc[1]<col_num:
            if dir[now_loc[0], now_loc[1]] == 2 ** ((i + 4) % 8):
                up_cell.append(now_loc)

    return up_cell

def lonlat_to_pixel(geo, lon, lat):
    """根据经纬度坐标获取栅格数据中的像素位置"""

    # 获取仿射变换参数
    origin_x, pixel_width, _, origin_y, _, pixel_height = geo

    # 计算像素坐标（行、列）
    col = int((lon - origin_x) / pixel_width)
    row = int((lat - origin_y) / pixel_height)

    return row, col

def pixel_to_lonlat(geo, row, col):
    """根据行列号获取经纬度坐标"""

    # 取仿射变换参数
    origin_x, pixel_width, _, origin_y, _, pixel_height = geo

    # 计算经纬度
    lon = origin_x + col * pixel_width
    lat = origin_y + row * pixel_height  # 注意 pixel_height 通常是负值

    return lon, lat



def callback(result):
    print(result)

def get_basin_FAB(fdir,lon,lat,geo,proj,outfile):
    """ """
    row,col = lonlat_to_pixel(geo,lon,lat)
    vis = np.zeros((fdir.shape[0],fdir.shape[1]))
    max_row = row
    min_row = row
    max_col = col
    min_col = col
    # basinCells = [(row,col)]
    # popCells = [(row,col)]
    basinCells = [(row, col), (row + 1, col + 1), (row, col + 1), (row - 1, col - 1), (row - 1, col),
                  (row - 1, col + 1), (row, col - 1), (row + 1, col - 1), (row + 1, col)]
    popCells = [(row, col), (row + 1, col + 1), (row, col + 1), (row - 1, col - 1), (row - 1, col), (row - 1, col + 1),
                (row, col - 1), (row + 1, col - 1), (row + 1, col)]
    for basinCell in basinCells:
        vis[basinCell[0],basinCell[1]] = 1

    while popCells:
        popCell = popCells.pop()
        vis[popCell[0], popCell[1]] = 1
        upCells = get_rever_D8(fdir,popCell[0],popCell[1])

        # basinCells += upCells
        # popCells += upCells

        for upCell in upCells:
            if vis[upCell[0],upCell[1]] == 1:
                continue
            basinCells.append(upCell)
            popCells.append(upCell)
            max_row = max(max_row,upCell[0])
            max_col = max(max_col,upCell[1])
            min_row = min(min_row, upCell[0])
            min_col = min(min_col, upCell[1])

    # 坐标转换至小坐标系
    new_row = max_row - min_row + 1
    new_col = max_col - min_col + 1
    basin = np.zeros((new_row,new_col))
    for basinCell in basinCells:
        basin[basinCell[0]-min_row,basinCell[1]-min_col] = 1

    new_lon,new_lat = pixel_to_lonlat(geo,min_row,min_col)
    new_geo = (new_lon,geo[1],geo[2],new_lat,geo[4],geo[5])

    Raster.save_raster(outfile,basin,proj,new_geo,gdal.GDT_Byte,0)
def get_basin_FAB1(basearea,fdir,acc,lon,lat,geo,proj,outfile):
    """ """
    row,col = lonlat_to_pixel(geo,lon,lat)
    # vis = np.zeros((fdir.shape[0],fdir.shape[1]),dtype = np.int8)
    max_row = row
    min_row = row
    max_col = col
    min_col = col
    basinCells = []
    # popCells = [(row,col)]
    outletCells = [(row,col)]
    initial_area = acc[row,col] * 30 * 30 / 1000000

    # basinCells = [(row, col), (row + 1, col + 1), (row, col + 1), (row - 1, col - 1), (row - 1, col),
    #               (row - 1, col + 1), (row, col - 1), (row + 1, col - 1), (row + 1, col)]
    # popCells = [(row, col), (row + 1, col + 1), (row, col + 1), (row - 1, col - 1), (row - 1, col), (row - 1, col + 1),
    #             (row, col - 1), (row + 1, col - 1), (row + 1, col)]
    # for basinCell in basinCells:
    #     vis[basinCell[0],basinCell[1]] = 1
    n = 0
    # 追溯下游，直到面积合适了
    outlet_cell = [row,col]


    n = 0
    while abs((initial_area - basearea) / basearea) > 0.2:# abs(((len(basinCells)*30*30/1000000) - basearea)/basearea) > 0.2:
        if initial_area > basearea:
            break
        nowDir = fdir[outlet_cell[0], outlet_cell[1]]
        if nowDir not in dmove_dic:
            break

        outlet_cell[0] += dmove_dic[nowDir][0]
        outlet_cell[1] += dmove_dic[nowDir][1]
        initial_area = acc[outlet_cell[0],outlet_cell[1]] * 30 * 30 / 1000000
        n += 1
        if n > 200:
            return [lon,lat,-1,-1,-1]

        # ----------------------------------------------------------- 校正提取流域一体化的代码 -----------------------------------------
        # temp_Cell = outletCells.pop()
        # popCells = [temp_Cell]
        # basinCells.append(temp_Cell)
        # if vis[temp_Cell[0],temp_Cell[1]] == 1:
        #     break
        # vis[temp_Cell[0],temp_Cell[1]] = 1
        #
        # while popCells:
        #     popCell = popCells.pop()
        #     max_row = max(max_row, popCell[0])
        #     max_col = max(max_col, popCell[1])
        #     min_row = min(min_row, popCell[0])
        #     min_col = min(min_col, popCell[1])
        #     vis[popCell[0], popCell[1]] = 1
        #     upCells = get_rever_D8(fdir,popCell[0],popCell[1])
        #
        #     # basinCells += upCells
        #     # popCells += upCells
        #
        #     for upCell in upCells:
        #         if vis[upCell[0],upCell[1]] == 1:
        #             continue
        #         basinCells.append(upCell)
        #         popCells.append(upCell)
        #
        # nowDir = fdir[temp_Cell[0], temp_Cell[1]]
        # if nowDir not in dmove_dic:
        #     break
        # outletCells.append((temp_Cell[0] + dmove_dic[nowDir][0], temp_Cell[1] + dmove_dic[nowDir][1]))
        # n += 1
        # ----------------------------------------------------------- 校正提取流域一体化的代码 -----------------------------------------

    new_lon, new_lat = pixel_to_lonlat(geo, outlet_cell[0], outlet_cell[1])
    return [new_lon,new_lat,initial_area,outlet_cell[0],outlet_cell[1]]


    # o_lon,o_lat = pixel_to_lonlat(geo,temp_Cell[0],temp_Cell[1])

    # 坐标转换至小坐标系
    # new_row = max_row - min_row + 1
    # new_col = max_col - min_col + 1
    # basin = np.zeros((new_row,new_col),dtype = np.int8)
    # for basinCell in basinCells:
    #     basin[basinCell[0]-min_row,basinCell[1]-min_col] = 1
    #
    # new_lon,new_lat = pixel_to_lonlat(geo,min_row,min_col)
    # new_geo = (new_lon,geo[1],geo[2],new_lat,geo[4],geo[5])
    #
    # Raster.save_raster(outfile,basin,proj,new_geo,gdal.GDT_Byte,0)
    # return [new_lon,new_lat]
def get_basin_FAB2(basearea,stream,Stream_nodata,fdir,lon,lat,geo,proj,outfile):

    """
    追溯到下游最近河网，如果面积大于实测面积，就沿河网上移
    :param basearea:
    :param stream:
    :param fdir:
    :param lon:
    :param lat:
    :param geo:
    :param proj:
    :param outfile:
    :return:
    """
    row,col = lonlat_to_pixel(geo,lon,lat)

    max_row = row
    min_row = row
    max_col = col
    min_col = col
    basinCells = []
    # popCells = [(row,col)]

    # basinCells = [(row, col), (row + 1, col + 1), (row, col + 1), (row - 1, col - 1), (row - 1, col),
    #               (row - 1, col + 1), (row, col - 1), (row + 1, col - 1), (row + 1, col)]
    # popCells = [(row, col), (row + 1, col + 1), (row, col + 1), (row - 1, col - 1), (row - 1, col), (row - 1, col + 1),
    #             (row, col - 1), (row + 1, col - 1), (row + 1, col)]
    # for basinCell in basinCells:
    #     vis[basinCell[0],basinCell[1]] = 1

    # F1:按流向寻找下游最近河网
    pop_cells=[(row,col)]
    # Vis[row,col]=1
    # print(pop_cells)
    while pop_cells:
        pop_cell=pop_cells.pop()
        # print(pop_cell)
        if stream[pop_cell[0],pop_cell[1]]!=Stream_nodata:
            # print(pop_cell)
            row,col = pop_cell
            break

        else:
            now_dir=fdir[pop_cell[0],pop_cell[1]]
            if now_dir in dmove_dic:
                next_cell = (pop_cell[0] + dmove_dic[now_dir][0], pop_cell[1] + dmove_dic[now_dir][1])
                pop_cells.insert(0, next_cell)

    outletCells = [(row, col)]
    n = 0
    # 追溯下游，直到面积合适了
    while abs(((len(basinCells)*30*30/1000000) - basearea)/basearea) > 0.2:
        if n>10:
            break
        basinCells = []
        vis = np.zeros((fdir.shape[0], fdir.shape[1]), dtype=np.int8)
        temp_Cell = outletCells.pop()
        popCells = [temp_Cell]
        basinCells.append(temp_Cell)
        if vis[temp_Cell[0],temp_Cell[1]] == 1:
            break
        vis[temp_Cell[0],temp_Cell[1]] = 1

        while popCells:
            popCell = popCells.pop()
            max_row = max(max_row, popCell[0])
            max_col = max(max_col, popCell[1])
            min_row = min(min_row, popCell[0])
            min_col = min(min_col, popCell[1])
            vis[popCell[0], popCell[1]] = 1
            upCells = get_rever_D8(fdir,popCell[0],popCell[1])

            # basinCells += upCells
            # popCells += upCells

            for upCell in upCells:
                if vis[upCell[0],upCell[1]] == 1:
                    continue
                basinCells.append(upCell)
                popCells.append(upCell)

        upCells1 = get_rever_D8(fdir, temp_Cell[0], temp_Cell[1])
        for upCell1 in upCells1:
            if stream[upCell1[0],upCell1[1]] != Stream_nodata:
                outletCells.append(upCell1)
        n += 1

    # o_lon,o_lat = pixel_to_lonlat(geo,temp_Cell[0],temp_Cell[1])

    # 坐标转换至小坐标系
    new_row = max_row - min_row + 1
    new_col = max_col - min_col + 1
    basin = np.zeros((new_row,new_col))
    for basinCell in basinCells:
        basin[basinCell[0]-min_row,basinCell[1]-min_col] = 1

    new_lon,new_lat = pixel_to_lonlat(geo,min_row,min_col)
    new_geo = (new_lon,geo[1],geo[2],new_lat,geo[4],geo[5])

    Raster.save_raster(outfile,basin,proj,new_geo,gdal.GDT_Byte,0)
    # return o_lon,o_lat
def get_lon_lat_area(fn,GRDC_id,outvenu):

    '''
    :param fn: GRDC.shp (每个要素矢量)
    :return:
    '''

    ds = ogr.Open(fn)
    layer = ds.GetLayer()
    feature = layer.GetFeature(0)

    # lon = feature.GetField("LONG_NEW")
    # lat = feature.GetField("LAT_NEW")
    GRDC_NO = feature.GetField("GRDC_NO")
    lon = feature.GetField("LONG_ORG")
    lat = feature.GetField("LAT_ORG")
    rep_area = feature.GetField("AREA")
    hyd_area = feature.GetField("AREA_HYS")
    lon_New = feature.GetField("LONG_NEW")
    lat_New = feature.GetField("LAT_NEW")

    # if not os.path.exists(outvenu):
    #     os.mkdir(outvenu)
        # os.chmod(outvenu,0o777)
    outfile1 = os.path.join(outvenu,str(GRDC_id)+'_FAB_ORG.tif')
    # get_basin_FAB(fdir,lon,lat,geo,proj,outfile)

    outfile2 = os.path.join(outvenu, str(GRDC_id) + '_FAB_NEW.tif')
    # get_basin_FAB(fdir, lon_New, lat_New, geo, proj, outfile)

    return [GRDC_NO,lon,lat,lon_New,lat_New,rep_area,hyd_area,outfile1,outfile2]

def get_basin_main(folder,Dir_file,acc_file,venu,continent = 'AS'):
    """
    根据GRDC提供的面积校正流域出口，按照下游追溯，知道接近最近的面积
    :param folder:存放GRDC矢量的文件夹
    :param inFile:
    :param dicFile:
    :return:
    """
    folder = "/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/statbas_shp_zip/"  #'/datanode05/zhangbin/NWEI/data/GRDC'
    if not os.path.exists(venu):
        os.mkdir(venu)
        # os.chmod(venu,0o777)

    outvenu = os.path.join(venu,'FAB_basins1')
    if not os.path.exists(outvenu):
        os.mkdir(outvenu)
        # os.chmod(outvenu,0o777)
    #
    fdir = Raster.get_raster(Dir_file)
    proj,geo,f_nodata = Raster.get_proj_geo_nodata(Dir_file)

    acc = Raster.get_raster(acc_file)

    # 读取站点列表
    # spath = "/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/SouthAmerica.txt"
    # spath = "/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/NorthAmerica.txt"
    # spath = "/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/Australia.txt"
    # spath = "/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/Siberia.txt"
    # spath = "/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/Europe.txt"
    # spath = "/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/Africa.txt"
    # spath = "/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/Arctic.txt"
    # spath = "/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/Asia.txt"
    # spath = "/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/Greenland.txt"
    point_path = {'AS':"/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/Asia.txt",
                  'GR':"/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/Greenland.txt",
                  "AR":"/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/Arctic.txt",
                  "AF":"/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/Africa.txt",
                  "EU":"/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/Europe.txt",
                  "SI":"/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/Siberia.txt",
                  "AU":"/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/Australia.txt",
                  "NA":"/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/NorthAmerica.txt",
                  "SA":"/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/SouthAmerica.txt"}
    spath = point_path[continent]
    stations = np.loadtxt(spath, dtype=np.int32)
    stations = stations.ravel()

    result = []
    paras = []
    for station in stations[:]:
        # print(station)
        print(station)
        fn = os.path.join(folder, "grdc_basins_smoothed_md_no_%d.shp" % station)
        temp = get_lon_lat_area(fn,station,outvenu)   # [GRDC_NO,lon,lat,lon_New,lat_New,rep_area,hyd_area,outfile1,outfile2]
        try:
            temp1 = get_basin_FAB1(float(temp[5]),fdir,acc,temp[3],temp[4],geo,proj,temp[8])

        except Exception as e:
            print(station,e)
            temp1 = [-1,-1,-1,-1,-1]

        try:
            temp2 = get_basin_FAB1(float(temp[5]), fdir,acc, temp[1], temp[2], geo, proj, temp[7])

        except Exception as e:
            print(station, e)
            temp2 = [-1,-1,-1,-1,-1]
        r = [station] + temp + temp1 + temp2
        print(r)
        result.append(r)

        # result.append(temp)
        # paras.append(temp)
        # break
    # , columns=["GRDC_NO",'lon','lat','lon_New','lat_New', "GRDC_AREA", "HYS_AREA"]
    df = pd.DataFrame(result)   # ,'o_lon_new','o_lat_new','o_lon','o_lat'
    df.to_csv(os.path.join(venu, 'Europe_GRDC_station.csv'))
    # os.chmod(os.path.join(venu, 'GRC_area.csv'),0o777)

def extract_basin(Dir_file,split_csv):
    """
    根据校正的流域出口提取流域
    :param split_csv:
    :return:
    """

    with open(split_csv,'r') as f:
        con = []
        reader = csv.reader(f)
        n = 0
        for i in reader:
            if n==0:
                n += 1
                continue
            if os.path.exists(i[10]):
                continue
            if float(i[26]) < 0:
                continue
            con.append([int(float(i[24])),int(float(i[25])),i[10]])

        f.close()
    # print(con)
    fdir = Raster.get_raster(Dir_file)
    proj, geo, f_nodata = Raster.get_proj_geo_nodata(Dir_file)
    row,col = fdir.shape

    paras = []
    for lon_lat in con:
        if os.path.exists(lon_lat[2]):
            continue
        paras.append([lon_lat[0],lon_lat[1],lon_lat[2]])
        # try:
        #     vis = np.zeros((row, col), dtype=np.int8)
        #     max_row = 0
        #     min_row = row
        #     max_col = 0
        #     min_col = col
        #     # row1, col1 = lonlat_to_pixel(geo, lon_lat[0], lon_lat[1])
        #     row1 = lon_lat[0]
        #     col1 = lon_lat[1]
        #     # print(row1,col1)
        #     popCells = [(row1,col1)]
        #     basinCells = [(row1,col1)]
        #
        #     if vis[row1,col1] == 1:
        #         break
        #     vis[row1,col1] = 1
        #
        #     while popCells:
        #         popCell = popCells.pop()
        #         max_row = max(max_row, popCell[0])
        #         max_col = max(max_col, popCell[1])
        #         min_row = min(min_row, popCell[0])
        #         min_col = min(min_col, popCell[1])
        #         vis[popCell[0], popCell[1]] = 1
        #         upCells = get_rever_D8(fdir,popCell[0],popCell[1])
        #
        #         # basinCells += upCells
        #         # popCells += upCells
        #
        #         for upCell in upCells:
        #             if vis[upCell[0],upCell[1]] == 1:
        #                 continue
        #             basinCells.append(upCell)
        #             popCells.append(upCell)
        #
        #
        #     # 坐标转换至小坐标系
        #     new_row = max_row - min_row + 1
        #     new_col = max_col - min_col + 1
        #     basin = np.zeros((new_row,new_col),dtype = np.int8)
        #     for basinCell in basinCells:
        #         basin[basinCell[0]-min_row,basinCell[1]-min_col] = 1
        #
        #     new_lon,new_lat = pixel_to_lonlat(geo,min_row,min_col)
        #     new_geo = (new_lon,geo[1],geo[2],new_lat,geo[4],geo[5])
        #
        #     Raster.save_raster(lon_lat[2],basin,proj,new_geo,gdal.GDT_Byte,0)
        # except Exception as e:
        #     print(lon_lat,e)

    po = Pool(10)
    for para in paras:
        # fdir, row1, col1,row,col,geo,proj,outfile
        po.apply_async(compare_NHD.extravt_basin_for_po,(fdir,para[0],para[1],row,col,geo,proj,para[2],))
    po.close()
    po.join()

def one_time_extract_basin(fdr_file,csv_file):

    """
    计算每个出口点的直接上游，后续再构建流域间的汇流关系，进而避免重复计算，来提高效率
    :return:
    """
    import numpy as np
    from collections import deque

    # 假设 D8 流向编码：
    # 1=E, 2=SE, 3=S, 4=SW, 5=W, 6=NW, 7=N, 8=NE
    # 或根据你的 D8 格式修改 delta_row, delta_col

    # 8 邻域对应的行列偏移
    delta_row = [0, 1, 1, 1, 0, -1, -1, -1]
    delta_col = [1, 1, 0, -1, -1, -1, 0, 1]

    def get_upstream_neighbors(fdir, row, col):
        """返回流入 (row, col) 的上游像元列表"""
        neighbors = []
        rows, cols = fdir.shape
        for i in range(8):
            r = row + delta_row[i]
            c = col + delta_col[i]
            if 0 <= r < rows and 0 <= c < cols:
                # 判断 (r,c) 的流向是否指向 (row,col)
                flow_dir = fdir[r, c]
                if flow_dir == ((i + 4) % 8 + 1):  # 反向 D8 判断
                    neighbors.append((r, c))
        return neighbors

    def multi_exit_basin(fdir, export_points):
        """
        多出口批处理回溯
        fdir: 2D numpy array, D8 流向
        export_points: list of tuples [(row, col, id), ...]
        return: basin_id 栅格
        """
        rows, cols = fdir.shape
        basin_id = np.zeros((rows,cols),dtype = np.int16)
        basin_id[:,:] = -9999
        # visited = np.zeros((rows,cols),dtype = np.int8)


        # 初始化队列，将出口入队
        for r, c, bid in export_points:
            if genral_functions.check_boundary(r,c,rows,cols):
                basin_id[r, c] = bid
                # visited[r, c] = 1

        # 初始化队列，将出口入队
        for cell in export_points:
            if not genral_functions.check_boundary(cell[0], cell[1], rows, cols):
                continue
            pop_cells = [cell]
            bid = cell[2]

            while pop_cells:

                pop_cell = pop_cells.pop()

                upCells = genral_functions.get_upstreamCells(pop_cell[0],pop_cell[1],fdir,rows,cols)
                for upCell in upCells:
                    if basin_id[upCell[0], upCell[1]] != -9999:
                        continue
                    pop_cells.append(upCell)
                    basin_id[upCell[0],upCell[1]] = bid


        return basin_id

    def construct_confluence(fdr,basin_file,points):
        """
        构建汇流关系表
        :param fdr:
        :param points:
        :return:
        """

        confluence = {}
        basin = Raster.get_raster(basin_file)
        for r, c, bid in points:
            basin_id = basin[r,c]
            now_fdr = fdr[r,c]

            if now_fdr not in dmove_dic:
                # confluence.setdefault(now_fdr,-9999)
                continue
            next_cell = (r+dmove_dic[now_fdr][0],c+dmove_dic[now_fdr][1])
            next_id = basin[next_cell[0],next_cell[1]]

            confluence.setdefault(next_id,[]).append(basin_id)

        return confluence

    def raster2shp(raster,gpkg):
        from osgeo import gdal, ogr, osr

        layer_name = "basin_poly"

        ds = gdal.Open(raster)
        band = ds.GetRasterBand(1)

        nodata = band.GetNoDataValue()

        # 创建 GPKG
        driver = ogr.GetDriverByName("GPKG")
        out_ds = driver.CreateDataSource(gpkg)

        # 空间参考
        srs = osr.SpatialReference()
        srs.ImportFromWkt(ds.GetProjection())

        # 创建图层
        layer = out_ds.CreateLayer(
            layer_name,
            srs=srs,
            geom_type=ogr.wkbPolygon
        )

        # ⚠️ 关键：字段用于接收“栅格值”
        field_defn = ogr.FieldDefn("basin_id", ogr.OFTInteger)
        layer.CreateField(field_defn)

        # Polygonize
        gdal.Polygonize(
            band,  # 输入栅格
            None,  # mask（None = 忽略 NoData）
            layer,  # 输出图层
            0,  # 字段索引（basin_id）
            [],  # options
            callback=None
        )

        # 释放
        out_ds = None
        ds = None

    def get_upstream_from_pkl(tree, huc12):
        # with open(treeFile, 'rb') as f:
        #     con = pickle.load(f)
        #     f.close()

        ups = [huc12]
        result = []
        while ups:
            nowid = ups.pop()
            result.append(nowid)
            if nowid in tree:
                ups += tree[nowid]

        return result

    def create_output_shp_from_gpkg(basinList, filePath, outPath):
        """
        输出流域范围
        :param basinList: [HUC12]
        :param filePath: NHDPlus.gdb
        :param outPath: 输出路径
        :return:
        """

        # Step 1: 加载 NHDFlowline 图层
        # driver = ogr.GetDriverByName("OpenFileGDB")
        # gdb = driver.Open(filePath, 0)  # 只读模式
        huc12_layer = "basin_poly"
        ds = ogr.Open(filePath)
        # inLayer = gdb.GetLayer(huc12_layer)
        inLayer = ds.GetLayer(0)
        srs = inLayer.GetSpatialRef()
        featureDefn = inLayer.GetLayerDefn()

        shpDriver = ogr.GetDriverByName("ESRI Shapefile")
        outDs = shpDriver.CreateDataSource(outPath)
        outLayer = outDs.CreateLayer("data", srs=srs, geom_type=ogr.wkbMultiPolygon)
        # Add field values from input Layer
        for i in range(0, featureDefn.GetFieldCount()):
            fieldDefn = featureDefn.GetFieldDefn(i)
            outLayer.CreateField(fieldDefn)

        # id_filter = "huc12 IN ({})".format(", ".join(map(str, basinList)))
        # # Apply the filter to the layer
        # inLayer.SetAttributeFilter(id_filter)
        area = 0
        for feature in inLayer:
            if feature.GetField("basin_id") in basinList:
                outLayer.CreateFeature(feature)
                area += feature.GetField("area_km2")
        # print(Total_area)
        outLayer.SyncToDisk()
        outDs.Destroy()

        return area
    def cal_area(gpkg_file):
        import geopandas as gpd

        gdf = gpd.read_file(gpkg_file, layer="basin_poly")

        gdf_area = gdf.to_crs("EPSG:6933")  # World Cylindrical Equal Area

        gdf["area_km2"] = gdf_area.geometry.area / 1e6

        gdf.to_file(
            gpkg_file,
            layer="basin_poly",
            driver="GPKG"
        )
    def cal_area_shp(shp_file):
        import geopandas as gpd

        gdf = gpd.read_file(shp_file)

        gdf_area = gdf.to_crs("EPSG:6933")  # World Cylindrical Equal Area

        gdf["area_km2"] = gdf_area.geometry.area / 1e6

        return gdf["area_km2"].sum()

    def get_upstream_basin_area_shp(con,gpkg_file,treefile,folder,vectorPath,result_csv):
        """
        先根据汇流树把上游id构建出来
        再输出完整流域的shp
        再计算CSI和RE
        :param con:
        :param gpkg_file:
        :param treefile:
        :return:
        """

        if not os.path.exists(vectorPath):
            os.mkdir(vectorPath)
        with open(treefile, 'rb') as f:
            tree = pickle.load(f)
            f.close()

        result = []
        for i in con:
            info = [-1,i[3], i[1]]  # i[1] 是fid，0，1，2，3   i[3]是GRDC的id

            fn = os.path.join(folder, "grdc_basins_smoothed_md_no_%d.shp" % int(float(info[1])))
            temp = list(get_lon_lat_area(fn, info[1], ""))
            # FABbasinORG = os.path.join(FABbasinPath, "%d_FAB_ORG.tif" % int(float(info[1])))  #
            # FABbasinNEW = os.path.join(FABbasinPath, "%d_FAB_NEW.tif" % int(float(info[1])))
            FABbasinOGRvector = os.path.join(vectorPath, "%d_FAB_OGR.shp" % int(float(info[2])))
            FABbasinNEWvector = os.path.join(vectorPath, "%d_FAB_NEW.shp" % int(float(info[2])))
            # if not os.path.exists(FABbasinOGRvector):
                # paras.append([FABbasinORG,FABbasinOGRvector])
            ups = get_upstream_from_pkl(tree,float(info[2]))
            area = create_output_shp_from_gpkg(ups,gpkg_file,FABbasinOGRvector)
            basearea = float(temp[6])
            re = (area - basearea) / basearea


            if os.path.exists(FABbasinOGRvector):
                # paras.append([FABbasinORG, FABbasinOGRvector])

                try:
                    if os.path.exists(FABbasinOGRvector):
                        temp.append(area)
                        temp.append(CSI(fn, FABbasinOGRvector))
                        temp.append(re)

                    else:
                        temp.append(str(0))
                        temp.append(str(0))
                        temp.append(str(0))
                except Exception as e:
                    print(FABbasinOGRvector, e)
                    temp.append(str(0))
                    temp.append(str(0))
                    temp.append(str(0))
                result.append(temp)
            # break
        with open(result_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(result)
            f.close()


    def get_upstream_basin_area_shp_NHD(con,gpkg_file,treefile,vectorPath,result_csv):
        """
        先根据汇流树把上游id构建出来
        再输出完整流域的shp
        再计算CSI和RE
        :param con:
        :param gpkg_file:
        :param treefile:
        :return:
        """
        NHD_folder = "/datanode05/xiefy/zhangb/FAB_Hydro/NHD_valid/NHDPlusUPS/"
        if not os.path.exists(vectorPath):
            os.mkdir(vectorPath)
        with open(treefile, 'rb') as f:
            tree = pickle.load(f)
            f.close()
        area_hu = {}
        with open("/datanode05/xiefy/zhangb/FAB_Hydro/NHD_valid/NHD_lon_lat_HU.csv",'r') as f:
            reader = csv.reader(f)
            for i in reader:
                temp_infos = i[0].split('_')
                area_hu.setdefault((temp_infos[0],temp_infos[1]),temp_infos[2])
            f.close()


        result = []
        for i in con:
            info = [-1,i[3], i[0]]  # i[1] 是fid，0，1，2，3   i[3]是GRDC的id

            # old_HU = i[1].rstrip()
            # HU = old_HU.zfill(12)
            # print(old_HU,HU)
            HU = area_hu[(str(i[2]),str(i[3]))]


            fn = os.path.join(NHD_folder, str(HU)+'.shp')
            # print(fn)
            temp = [i[0],i[1],i[2],i[3],i[4],i[5],i[6]]


            FABbasinOGRvector = os.path.join(vectorPath, "%d_FAB_OGR.shp" % int(float(info[2])))


            if not os.path.exists(fn):
                print(HU,fn)
            try:
                NHDarea = cal_area_shp(fn)
                temp.append(NHDarea)
            except:
                temp.append(-1)
            result.append(temp)

            # ups = get_upstream_from_pkl(tree,float(info[2]))
            # try:
            #     # area = create_output_shp_from_gpkg(ups,gpkg_file,FABbasinOGRvector)
            #     area = cal_area_shp(FABbasinOGRvector)
            #     basearea = float(temp[6])
            #     re = (area - basearea) / basearea
            #     pass
            # except:
            #     continue
            # if not os.path.exists(fn):
            #     continue
            #
            # if os.path.exists(FABbasinOGRvector):
            #     # paras.append([FABbasinORG, FABbasinOGRvector])
            #
            #     try:
            #         if os.path.exists(FABbasinOGRvector):
            #             temp.append(area)
            #             temp.append(CSI(fn, FABbasinOGRvector))
            #             temp.append(re)
            #
            #         else:
            #             temp.append(str(0))
            #             temp.append(str(0))
            #             temp.append(str(0))
            #     except Exception as e:
            #         print(FABbasinOGRvector, e)
            #         temp.append(str(0))
            #         temp.append(str(0))
            #         temp.append(str(0))
            #     result.append(temp)
            # # break
        with open(result_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(result)
            f.close()


    folder = "/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/statbas_shp_zip"
    NHD_infos = []
    NHD_infos_all = []
    with open(csv_file,'r') as f:
        con = []
        all_infos = []
        reader = csv.reader(f)
        n = 0
        for i in reader:
            if n==0:
                n += 1
                continue
            # if os.path.exists(i[10]):
            #     continue
            # if float(i[26]) < 0:
            #     continue
            # con.append([int(float(i[24])),int(float(i[25])),float(i[1])])
            if float(i[7])>0:
                NHD_infos.append([int(float(i[7])),int(float(i[8])),float(i[0])])
                NHD_infos_all.append(i)
            all_infos.append(i)
            n += 1

    # fdr = Raster.get_raster(fdr_file)
    # proj,geo,nodata = Raster.get_proj_geo_nodata(fdr_file)

    basevenu = os.path.dirname(csv_file)
    out_one_time_basin = os.path.join(basevenu, 'out_basin_one_time.tif')
    out_one_time_basin_gpkg = os.path.join(basevenu, 'out_basin_one_time.gpkg')
    #
    out_one_time_basin_confluence = os.path.join(basevenu, 'basin_confluence.pkl')
    #
    vectorPath = os.path.join(basevenu, 'FAB_basins_vector')
    result_csv = os.path.join(basevenu, 'NHD_AREA.csv')

    # ------------------------------------ GRDC 验证 -------------------------------------------------------
    # basin_result = multi_exit_basin(fdr, con)
    # Raster.save_raster(out_one_time_basin,basin_result,proj,geo,gdal.GDT_Int16,-9999)
    #
    # conflluence = construct_confluence(fdr,out_one_time_basin,con)
    # with open(out_one_time_basin_confluence,'wb') as f:
    #     pickle.dump(conflluence,f)
    #
    #     f.close()
    # ------------------------------------ GRDC 验证 -------------------------------------------------------


    # # ------------------------------------ NHD 验证 -------------------------------------------------------
    # basin_result = multi_exit_basin(fdr, NHD_infos)
    # Raster.save_raster(out_one_time_basin, basin_result, proj, geo, gdal.GDT_Int16, -9999)
    #
    # conflluence = construct_confluence(fdr, out_one_time_basin, NHD_infos)
    # with open(out_one_time_basin_confluence, 'wb') as f:
    #     pickle.dump(conflluence, f)
    #
    #     f.close()
    # ------------------------------------ NHD 验证 -------------------------------------------------------
    # raster2shp(out_one_time_basin,out_one_time_basin_gpkg)
    # cal_area(out_one_time_basin_gpkg)

    # get_upstream_basin_area_shp(all_infos,out_one_time_basin_gpkg,out_one_time_basin_confluence,folder,vectorPath,result_csv)  # GRDC

    get_upstream_basin_area_shp_NHD(NHD_infos_all,out_one_time_basin_gpkg,out_one_time_basin_confluence,vectorPath,result_csv)  # NHD

    print("完成多出口流域回溯")
def CSI(GRDC_file,NWEI_file,PfafId = '1'):
    """
    计算单次流域面积CSI
    :param GRDC_file: shp
    :param NWEI_file: FAB_Basin.shp
    :return:
    """

    # 加载矢量文件
    file1 = gpd.read_file(GRDC_file)
    file2 = gpd.read_file(NWEI_file)

    # 投影到同一坐标系，确保单位一致
    # file1 = file1.to_crs("EPSG:3857")  # 3857
    # file2 = file2.to_crs("EPSG:3857")


    # NHD投影
    file1 = file1.to_crs("EPSG:5070")  # 3857 "EPSG:5070"
    file2 = file2.to_crs("EPSG:5070")

    # 计算交集
    intersection = gpd.overlay(file1, file2, how="intersection")

    # 计算面积
    intersection_area = intersection['geometry'].area.sum()
    file1_area = file1['geometry'].area.sum()
    file2_area = file2['geometry'].area.sum()

    # 计算重叠度
    overlap_ratio_1 = intersection_area / file1_area  # file1 的重叠度
    overlap_ratio_2 = intersection_area / file2_area  # file2 的重叠度

    # 并集面积和重叠系数（可选）
    union = gpd.overlay(file1, file2, how="union")
    union_area = union['geometry'].area.sum()
    overlap_coefficient = intersection_area / union_area
    # print(overlap_coefficient)
    # print(f"File1 的重叠度: {overlap_ratio_1:.2%}")
    # print(f"File2 的重叠度: {overlap_ratio_2:.2%}")
    print(f"重叠系数: {overlap_coefficient:.2%}")
    # 计算
    # return [str(PfafId),file1_area,file2_area,overlap_coefficient]
    return overlap_coefficient
def cal_area(raster_file,vector_file):

    # import whitebox
    # wbt = whitebox.WhiteboxTools()
    if not os.path.exists(vector_file):
        sink.wbt.raster_to_vector_polygons(raster_file,vector_file)
    # wbt.raster_to_vector_polygons(raster_file, vector_file)


    # 读取矢量数据
    gdf = gpd.read_file(vector_file)

    # 定义 WGS84 地理坐标系
    geod = Geod(ellps="WGS84")

    # 计算每个多边形的面积（单位：平方米）
    gdf['area_m2'] = gdf.geometry.apply(lambda geom: abs(geod.geometry_area_perimeter(geom)[0]))

    # 以平方千米为单位
    gdf['area_km2'] = gdf['area_m2'] / 1e6

    # 显示结果
    # print(gdf['area_m2'] / 1e6)

    return gdf['area_km2'].sum()

def raster2Feature(raster_file,vector_file):

    # sink.wbt.raster_to_vector_polygons(raster_file, vector_file)
    from osgeo import gdal, ogr, osr

    # 输入栅格
    raster_path = raster_file
    raster_ds = gdal.Open(raster_path)
    band = raster_ds.GetRasterBand(1)

    # 输出矢量
    vector_path = vector_file
    driver = ogr.GetDriverByName("ESRI Shapefile")
    driver.DeleteDataSource(vector_path)
    vector_ds = driver.CreateDataSource(vector_path)
    srs = osr.SpatialReference()
    srs.ImportFromWkt(raster_ds.GetProjection())
    layer = vector_ds.CreateLayer("polygonized", srs=srs)

    # 添加字段
    field_defn = ogr.FieldDefn("DN", ogr.OFTInteger)
    layer.CreateField(field_defn)

    # 栅格转矢量
    gdal.Polygonize(band, band, layer, 0, [], callback=None)

    vector_ds = None
    raster_ds = None
def sbatch_cal_area(venu,folder="/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/statbas_shp_zip"):
    """
    计算得到的流域交集CSI和面积
    :param venu:
    :param folder:
    :return:
    """
    baseName = os.path.basename(venu)
    stationfile = os.path.join(venu,'GRDC_station.csv')
    FABbasinPath = os.path.join(venu,'FAB_basins1')
    vectorPath =  os.path.join(venu,'FAB_basins_vector')
    result_csv = os.path.join(venu,'CSI.csv')
    if not os.path.exists(vectorPath):
        os.mkdir(vectorPath)
        # os.chmod(vectorPath,0o777)

    n = 0
    with open(stationfile,'r') as f:
        con = []
        reader = csv.reader(f)
        for i in reader:
            if n==0:
                n += 1
                continue
            con.append(i)
        f.close()

    # ------------------ 根据txt的站点 ------------------
    # spath = "/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/Asia.txt"
    # spath = "/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/Greenland.txt"
    # stations = np.loadtxt(spath, dtype=np.int32)
    # con = stations.ravel()
    # ------------------ 根据txt的站点 ------------------
    paras = []
    result = []
    for i in con:
        info = [-1,i[3]]

        fn = os.path.join(folder, "grdc_basins_smoothed_md_no_%d.shp" % int(float(info[1])))
        temp = list(get_lon_lat_area(fn, info[1], ""))
        FABbasinORG = os.path.join(FABbasinPath,"%d_FAB_ORG.tif" % int(float(info[1])))  #
        FABbasinNEW = os.path.join(FABbasinPath, "%d_FAB_NEW.tif" % int(float(info[1])))
        FABbasinOGRvector = os.path.join(vectorPath,"%d_FAB_OGR.shp" % int(float(info[1])))
        FABbasinNEWvector = os.path.join(vectorPath, "%d_FAB_NEW.shp" % int(float(info[1])))
        # if not os.path.exists(FABbasinOGRvector):
        #     paras.append([FABbasinORG,FABbasinOGRvector])
            # raster2Feature(FABbasinORG,FABbasinOGRvector)
        if os.path.exists(FABbasinORG):
            paras.append([FABbasinORG, FABbasinOGRvector])

            try:
                if os.path.exists(FABbasinORG):
                    temp.append(cal_area(FABbasinORG,FABbasinOGRvector))
                    temp.append(CSI(fn, FABbasinOGRvector))
                else:
                    temp.append(str(0))
                    temp.append(str(0))
            except Exception as e:
                print(FABbasinORG,e)
                temp.append(str(0))
                temp.append(str(0))
            result.append(temp)
    with open(result_csv,'w',newline='') as f:
        writer = csv.writer(f)
        writer.writerows(result)
        f.close()

    # po = Pool(20)
    # for para in paras:
    #     po.apply_async(raster2Feature,(para[0],para[1],))
    #
    # po.close()
    # po.join()
def post_GRDC_RE(csv_file):
    """
    处理split后的csv，计算最小的RE
    :param csv_file:
    :return:
    """

    df = pd.read_csv(csv_file)
    df['lon'] = 0
    df['lat'] = 0
    df['row'] = 0
    df['col'] = 0
    df['Area'] = 0
    df['Reliative_error'] = 0

    for index, row in df.iterrows():
        print(index,row)
        baseArea = row['7']
        area1 = row['12']
        area2 = row['17']
        re1 = (area1 - baseArea)/(area1 + baseArea)/2
        re2 = (area2 - baseArea)/(area2 + baseArea)/2

        if abs(re1) > abs(re2):
            df['lon'][index] = row['15']
            df['lat'][index] = row['16']
            df['Area'][index] = row['17']
            df['row'][index] = row['18']
            df['col'][index] = row['19']
            df['Reliative_error'][index] = re2
        else:
            df['lon'][index] = row['10']
            df['lat'][index] = row['11']
            df['Area'][index] = row['12']
            df['row'][index] = row['13']
            df['col'][index] = row['14']
            df['Reliative_error'][index] = re1

    # print(df)

    plt.scatter(df['7'],df['Area'])
    plt.show()
    #
    df.to_csv(csv_file)

def get_CSI_lower_05(csv_file,acc_file,shp_venu,outvenu):
    """
    根据CSI.csv找到CSI小于0.5的流域ID，来将shp复制到一个单独的文件夹，用于后续处理
    :param csv_file:
    :return:
    """
    if not os.path.exists(outvenu):
        os.mkdir(outvenu)
    with open(csv_file,'r') as f:
        reader = csv.reader(f)
        for i in reader:
            if float(i[10]) < 0.5:
                basePath = os.path.join(shp_venu,str(i[0])+"_FAB_OGR.shp")
                basePath1 = os.path.join(shp_venu, str(i[0]) + "_FAB_OGR.shx")
                basePath2 = os.path.join(shp_venu, str(i[0]) + "_FAB_OGR.dbf")
                if not os.path.exists(os.path.join(outvenu,str(i[0])+"_FAB_OGR.shp")):
                    shutil.copyfile(basePath,os.path.join(outvenu,str(i[0])+"_FAB_OGR.shp"))
                    shutil.copyfile(basePath1, os.path.join(outvenu, str(i[0]) + "_FAB_OGR.shx"))
                    shutil.copyfile(basePath2, os.path.join(outvenu, str(i[0]) + "_FAB_OGR.dbf"))
                Find.clip(acc_file,basePath,os.path.join(outvenu,str(i[0])+'acc.tif'))

def merge_continent_CSI_csv(venu):
    """
    把文件夹下各个大洲的数据进行合并，到global_CSI
    :param venu:
    :param out_global_csv:
    :return:
    """

    out_global_csv = os.path.join(venu,'Global_CSI.csv')

    if os.path.exists(out_global_csv):
        os.remove(out_global_csv)
    infos = []
    for file in os.listdir(venu):
        fila_path = os.path.join(venu,file)

        with open(fila_path,'r') as f:
            reader = csv.reader(f)
            for i in reader:
                if 0.1 < abs((float(i[6]) /float(i[9]))) < 10:
                    infos.append(i)
            f.close()

    with open(out_global_csv,'w',newline='') as f:
        writer = csv.writer(f)
        writer.writerows(infos)
        f.close()




if __name__ == "__main__":

    # cal_area(r'F:\青藏高原水体数据集\DATA\result\3846510_FAB_NEW.tif',
    #          r'F:\青藏高原水体数据集\DATA\result\3846510_FAB_NEW.shp')

    # CSI(r'F:\全球数据集\图表\素材\画图data\GRDC\statbas_shp_zip\grdc_basins_smoothed_md_no_3846510.shp',r'F:\青藏高原水体数据集\DATA\result\3846510_FAB_NEW.shp')
    # post_GRDC_RE(r'F:\青藏高原水体数据集\New_DATA\GRDC\new\Africa_GRDC_station.csv')
    # post_GRDC_RE(r'F:\青藏高原水体数据集\New_DATA\GRDC\new\Australia_GRDC_station.csv')
    # post_GRDC_RE(r'F:\青藏高原水体数据集\New_DATA\GRDC\new\Siberia_GRDC_station.csv')
    # post_GRDC_RE(r'F:\青藏高原水体数据集\New_DATA\GRDC\Asia\Asia_GRDC_station.csv')
    # post_GRDC_RE(r'F:\青藏高原水体数据集\New_DATA\GRDC\new\NorthAmerica_GRDC_station.csv')
    # post_GRDC_RE(r'F:\青藏高原水体数据集\New_DATA\GRDC\SouthAmerica\SouthAmerica_GRDC_station.csv')
    # post_GRDC_RE(r'F:\青藏高原水体数据集\New_DATA\GRDC\new\Arctic_GRDC_station.csv')
    post_GRDC_RE(r'F:\青藏高原水体数据集\New_DATA\GRDC\new\Europe_GRDC_station.csv')

    merge_continent_CSI_csv(r'F:\青藏高原水体数据集\New_DATA\GRDC\制图\各大洲GRDC验证csv')