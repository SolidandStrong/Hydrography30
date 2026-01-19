# -*- coding: utf-8 -*-
"""
@Time ： 2024/9/20 21:31
@Auth ：
@File ：genral_functions.py
@IDE ：PyCharm
"""
from osgeo import gdal
import numpy as np
import os
import time
import math
from multiprocessing import Pool
import rasterio
from shutil import copy
from rasterio.enums import Compression


dmove=[(0,1),(1,1),(1,0),(1,-1),(0,-1),(-1,-1),(-1,0),(-1,1)]
dmove_dic = {1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1), 16: (0, -1), 32: (-1, -1), 64: (-1, 0), 128: (-1, 1)}
cellSize = 30 # Size of cell

import heapq
import queue
class MaxPriorityQueue:
    def __init__(self):
        self.heap = []

    def push(self, x, y, H):
        # 将 H 取负，以便实现最大堆
        heapq.heappush(self.heap, (-H, x, y))

    def pop(self):
        # 弹出最大元素并恢复 H 的正值
        H, x, y = heapq.heappop(self.heap)
        return (x, y, -H)

    def peek(self):
        # 查看最大元素，不移除
        if self.heap:
            H, x, y = self.heap[0]
            return (x, y, -H)
        return None

    def is_empty(self):
        return len(self.heap) == 0

    def size(self):
        return len(self.heap)
class QueueUsingQueue:
    def __init__(self):
        self.queue = queue.Queue()
    def enqueue(self, x,y,H):
        # 入队操作：将元素放入队列
        self.queue.put((x,y,H))

    def dequeue(self):
        # 出队操作：从队列前端移除元素
        if not self.is_empty():
            return self.queue.get()
        else:
            return None

    def peek(self):
        # 查看队列前端元素
        if not self.is_empty():
            item = self.queue.get()
            self.enqueue(item[0],item[1],item[2])  # 将元素重新放回队列
            return item
        return None

    def is_empty(self):
        # 判断队列是否为空
        return self.queue.empty()

    def size(self):
        # 返回队列的大小
        return self.queue.qsize()

class MinPriorityQueue:
    def __init__(self):
        self.heap = []

    def push(self, x, y, H):
        # 将元素插入堆中，并按照 H 升序排序
        heapq.heappush(self.heap, (H, x, y))

    def pop(self):
        # 弹出最小 H 元素
        if not self.is_empty():
            H, x, y = heapq.heappop(self.heap)
            return (x, y, H)
        else:
            return None

    def peek(self):
        # 查看最小 H 元素（不弹出）
        if not self.is_empty():
            H, x, y = self.heap[0]
            return (x, y, H)
        return None

    def is_empty(self):
        # 判断队列是否为空
        return len(self.heap) == 0

    def size(self):
        # 返回队列的大小
        return len(self.heap)









def get_upstreamCells(i,j,fdir,rows,cols):
    """
    寻找8邻域内的上游栅格
    :param i:
    :param j:
    :param fdir:
    :return:
    """
    results = []
    for k in range(8):
        nextCell = (i+dmove[k][0],j+dmove[k][1])
        if not check_boundary(nextCell[0],nextCell[1],rows,cols):
            continue
        nextfdir = fdir[nextCell[0],nextCell[1]]
        if 2**((k+4)%8) == nextfdir:
            results.append(nextCell)

    return results

def check_boundary(x,y,rows,cols):
    if 0<=x<rows and 0<=y<cols:
        return True
    return False

def check_sink_boundary(cell,sink,sinkNodata):
    """
    To check whether a sink cell is at its boundary
    :param cell: (i,j)
    :param sink: [[],[]]
    :param sinkNodata:
    :return:
    """
    row,col = sink.shape

    for k in range(8):
        nextCell = (cell[0]+dmove[k][0],cell[1]+dmove[k][1])
        if not check_boundary(nextCell[0],nextCell[1],row,col):
            return 1

        if sink[nextCell[0],nextCell[1]] == sinkNodata:
            return 1
    return 0

def check_sink_outlet(cell,dir,sink,sinkNodata):
    """
    Check whether a cell is outlet: 边界栅格的下个点是nodata或者流出sink则是outlet，赋值1
    :param cell:
    :param dir:
    :param sink:
    :param sinkNodata:
    :return:
    """

    row,col =sink.shape
    nowDir = dir[cell[0],cell[1]]
    if nowDir in dmove_dic:
        nextCell = (cell[0]+dmove_dic[nowDir][0],cell[1]+dmove_dic[nowDir][1])
        if not check_boundary(nextCell[0],nextCell[1],row,col):
            return False
        if sink[nextCell[0],nextCell[1]] == sinkNodata:
            return True
    return False

def find_sink_outlet(sink,dem,sinkCells,sinkId):

    row,col = sink.shape
    maxHCell = sinkCells[0]
    for cell in sinkCells:
        if cell[3] == 0:
            continue
        if cell[2] > maxHCell[2]:
            maxHCell = cell

    # 给这种情况的出口点重新计算流向
    maxUpslope = -1
    returnDir = -1
    for k in range(8):
        nextCell = (maxHCell[0]+dmove[k][0],maxHCell[1]+dmove[k][1])
        if not check_boundary(nextCell[0],nextCell[1],row,col):
            continue
        if sink[nextCell[0],nextCell[1]] == sinkId:
            continue
        nextCellH = dem[nextCell[0], nextCell[1]]

        if k in [0, 2, 4, 6]:
            factor = 1
        else:
            factor = math.sqrt(2)

        nowUpslope = (dem[maxHCell[0],maxHCell[1]] - nextCellH) / factor

        if nowUpslope > maxUpslope:
            returnDir = 2 ** k
            maxUpslope = nowUpslope

    return maxHCell,returnDir

def d8(cell,dem,vis,sink):
    """
    sink的nodata必须是0
    :param cell:
    :param dem:
    :param vis:
    :param sink:
    :return:
    """
    row,col = dem.shape
    cellH = dem[cell[0],cell[1]]
    maxUpslope = -1
    returnDir = -1
    for k in range(8):
       nextCell = (cell[0]+dmove[k][0],cell[1]+dmove[k][1])

       if not check_boundary(nextCell[0],nextCell[1],row,col):
           continue
       # if sink[nextCell[0], nextCell[1]] == 0:
       #     # 不计算sink外的坡度，保证水流都是流向sink内的
       #     continue
       if vis[nextCell[0],nextCell[1]] == 0:
           # 未得到流向的栅格
           continue

       nextCellH = dem[nextCell[0],nextCell[1]]

       if nextCellH > cellH:
           # 只需要高程小于中心高程的邻域栅格
           continue

       if k in [0,2,4,6]:
           factor = 1
       else:
           factor = math.sqrt(2)

       nowUpslope = (cellH - nextCellH)/factor
       if nowUpslope > maxUpslope:
           returnDir = 2**k
           maxUpslope = nowUpslope

    return returnDir

def reverseD8(cell,dem,vis,sink):


    row,col = dem.shape
    cellH = dem[cell[0],cell[1]]
    minUpslope = math.inf
    returnDir = -1
    for k in range(8):
       nextCell = (cell[0]+dmove[k][0],cell[1]+dmove[k][1])

       if not check_boundary(nextCell[0],nextCell[1],row,col):
           continue
       # if sink[nextCell[0],nextCell[1]] == 0:
       #     continue
       if vis[nextCell[0],nextCell[1]] == 0:
           # 未得到流向的栅格
           continue

       nextCellH = dem[nextCell[0],nextCell[1]]

       if nextCellH < cellH:
           # 只需要高程大于中心高程的邻域栅格
           continue

       if k in [0,2,4,6]:
           factor = 1
       else:
           factor = math.sqrt(2)

       nowUpslope = (nextCellH - cellH)/factor
       if nowUpslope < minUpslope:
           returnDir = 2**k
           minUpslope = nowUpslope

    return returnDir

def Get_extent_out(DEM,sink,DEM_nodata,sink_cells,row,col,EXTENT):

    extent_cells=[]
    for cell in sink_cells:
        for k in range(8):
            next_cell=(cell[0]+dmove[k][0],cell[1]+dmove[k][1])
            if check_boundary(next_cell[0],next_cell[1],row,col):
                if sink[next_cell[0],next_cell[1]]!=sink[cell[0],cell[1]] and DEM[next_cell[0],next_cell[1]]!=DEM_nodata:
                    # 会存在无效值，后面使用时，要注意判断无效值
                    EXTENT[next_cell[0],next_cell[1]]=sink[cell[0],cell[1]]
                    extent_cells.append((next_cell[0],next_cell[1],DEM[next_cell[0],next_cell[1]]))

    return EXTENT,extent_cells
def Get_extent_out1(DEM,sink,DEM_nodata,sink_cells,row,col,sinknodata):

    extent_cells=[]
    for cell in sink_cells:
        for k in range(8):
            next_cell=(cell[0]+dmove[k][0],cell[1]+dmove[k][1])
            if not check_boundary(next_cell[0],next_cell[1],row,col):
                break
            if sink[next_cell[0],next_cell[1]] == sinknodata:
                extent_cells.append((next_cell[0], next_cell[1], DEM[next_cell[0], next_cell[1]]))
                break
            if sink[next_cell[0],next_cell[1]]!=sink[cell[0],cell[1]] and DEM[next_cell[0],next_cell[1]]!=DEM_nodata:
                # 会存在无效值，后面使用时，要注意判断无效值
                extent_cells.append((next_cell[0],next_cell[1],DEM[next_cell[0],next_cell[1]]))

    return extent_cells
def getOutlet(extentCells,sink,dem,dir,demNodata,sinkId):
    row,col = dem.shape
    outlet = [-1, -1, math.inf]
    for cell in extentCells:
        if cell[2] == demNodata:
            continue
        if cell[2] >= outlet[2]:
            continue
        cellDir = dir[cell[0],cell[1]]
        if cellDir not in dmove_dic:
            return cell
        nextCell = (cell[0]+dmove_dic[cellDir][0],cell[1]+dmove_dic[cellDir][1])
        if not check_boundary(nextCell[0],nextCell[1],row,col):
            continue
        if sink[nextCell[0],nextCell[1]] == sinkId:
            continue
        outlet = cell

    return outlet



#  处理全球数据分块使用的代码
def readExtent(shp_venu):
    """
    读取每个要素的四至
    :param shp_venu:
    :return:
    """
    Extents = {}    # {PFaf_id:[left,right,down,top]}
    shpFiles = os.listdir(shp_venu)
    for shpName in shpFiles:
        if shpName.split('.')[-1] != 'shp':
            continue

        shpPath = os.path.join(shp_venu,shpName)

        ds = ogr.Open(shpPath)
        layer = ds.GetLayer(0)

        feature = layer.GetFeature(0)

        left, right, down, up = layer.GetExtent()
        print(left, right, down, up)
        tempExtent = [math.floor(left)-3,math.ceil(right)+2,math.floor(down)-3,math.ceil(up)+2]

        Extents.setdefault(int(feature.GetField('Pfaf_id')),tempExtent)

    return Extents
def selectImages(venu, extent, outpath, continent):
    target = os.path.join(outpath, str(continent))
    if not os.path.exists(target):
        os.mkdir(target)
        os.chmod(target, 0o777)

    # extent = (左经，右经，下维度，上维度)
    listdirs = os.listdir(venu)
    files = []
    for dirs in listdirs:
        tempFiles = os.listdir(os.path.join(venu, dirs))

        for file in tempFiles:
            if file.split('.')[-1] == 'tif':
                files.append(os.path.join(venu, dirs, file))

    for file in files:
        basename = os.path.basename(file)
        info = basename.split('_')[0]

        lfactor = 1 if info[3] == 'E' else -1
        ll = float(info[4:]) * lfactor
        rl = ll + 1

        lafactor = 1 if info[0] == 'N' else -1
        lla = float(info[1:3]) * lafactor
        rla = lla + 1

        if extent[0] <= ll <= extent[1] and extent[2] <= lla <= extent[3]:
            copy(file, target)
    return target
def mosaic(input_dir, outputVenu,LEWoutVenu, name,nodata_value=-9999):
    """
    批量镶嵌栅格文件，将多个栅格合并成一个。

    Parameters:
    input_dir (str): 输入栅格文件所在的文件夹路径。
    output_path (str): 输出镶嵌后栅格文件的保存路径。
    nodata_value (float, optional): 设置无数据值。如果未指定，保持默认值。
    """

    output_path = os.path.join(outputVenu,name+'.tif')
    LEWPath = os.path.join(LEWoutVenu,name+'.tif')
    # 找到文件夹中的所有tif文件
    raster_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.tif')]

    # 如果没有找到任何tif文件，报错
    if not raster_files:
        print("没有找到任何tif文件。")
        return

    # 使用gdal.Warp镶嵌栅格
    vrt_options = gdal.BuildVRTOptions(srcNodata=nodata_value, VRTNodata=nodata_value)
    vrt_path = os.path.join(input_dir, 'mosaic.vrt')
    gdal.BuildVRT(vrt_path, raster_files, options=vrt_options)

    # 将VRT文件保存为实际的GeoTIFF格式
    gdal.Translate(output_path, vrt_path)

    compress_raster(output_path,LEWPath)
def compress_raster(input_path, output_path, compression="LZW"):
    """
    使用压缩算法压缩栅格数据。

    Parameters:
    input_path (str): 原始栅格数据路径
    output_path (str): 压缩后的栅格数据保存路径
    compression (str): 压缩算法（LZW, DEFLATE等）
    """
    with rasterio.open(input_path) as src:
        profile = src.profile
        profile.update(compress=compression)

        with rasterio.open(output_path, 'w', **profile) as dst:
            for i in range(1, src.count + 1):
                data = src.read(i)
                dst.write(data, i)
def preprocess1(shp_venu,FABDEM_filesVenu,outPath_region):
    region_extent = readExtent(shp_venu)
    po = Pool(25)
    for Pfaf_id in region_extent:
        extent = region_extent[Pfaf_id]
        # selectImages(FABDEM_filesVenu,extent,outPath_region,Pfaf_id)
        po.apply_async(selectImages,(FABDEM_filesVenu,extent,outPath_region,Pfaf_id,))

    po.close()
    po.join()
def preprocess2(input_dir, outputVenu,LEWoutVenu,nodata_value=-9999):
    names = [str(i) for i in range(100)]
    Names = []
    for name in names:
        if not os.path.exists(os.path.join(input_dir,name)):
            continue
        Names.append(name)
    print(Names)
    po = Pool(25)
    for name in Names:

        po.apply_async(mosaic,(os.path.join(input_dir,name),outputVenu,LEWoutVenu,name,nodata_value,))

    po.close()
    po.join()


def COPY(demVenu,outpath,extent):
    # extent = (左经，右经，下低维度，上维度)
    files = []
    minLon = extent[0]
    maxLon = extent[1]
    minLat = extent[2]
    maxLat = extent[3]
    lons = [i for i in range(minLon,maxLon)]
    lats = [i for i in range(minLat,maxLat)]

    for lon in lons:
        for lat in lats:
            tempFileName = get_dem_name(lon, lat)
            tempFilePath = os.path.join(demVenu, tempFileName)
            if not os.path.exists(tempFilePath):
                continue
            files.append(tempFilePath)

    for file in files:
        copy(file, outpath)

def cal_globalsink(demVenu,outVenu):
    def Cal_dir_sink_wbt(DEM_file, Dir_file,sink_file):
        wbt.d8_pointer(DEM_file, Dir_file, esri_pntr=True)
        wbt.sink(DEM_file, sink_file)
        print("{:s}计算完成".format(DEM_file))

    Names = os.listdir(demVenu)
    dirVenuPath = os.path.join(outVenu,'globalDir')
    sinkVenuPath = os.path.join(outVenu,'globalSink')

    po = Pool(30)
    for name in Names:
        demFile = os.path.join(demVenu,name)
        dirPath = os.path.join(dirVenuPath, 'dir_' + name)
        sinkPath = os.path.join(sinkVenuPath, 'sink_' + name)
        po.apply_async(Cal_dir_sink_wbt, (demFile, dirPath, sinkPath,))

    po.close()
    po.join()


# global_sink.py
def get_sink_name(sinkFile):
    """
    根据文件名获取左下角的经纬度
    :param sinkFile:
    :return:
    """
    celllonlat = 0.00027777778
    sinkFile = os.path.basename(sinkFile)
    info = sinkFile.split('_')[1]
    if len(info) == 7:
        # 原始命名方式
        lon = int(info[4:])
        lat = int(info[1:3])
        # print(lat,lon)

        if info[0] == 'N':
            lat = int(float(lat))
        if info[0] == 'S':
            lat = int(float(lat)) * -1
        if info[3] == 'E':
            lon = int(float(lon))
        if info[3] == 'W':
            lon = int(float(lon)) * -1
        return lon,lat,lon+1-celllonlat,lat+1-celllonlat
    else:
        minlon = info[6:9]
        minlat = info[1:3]
        maxlon = info[9:]
        maxlat = info[3:5]

        if info[0] == 'N':
            minlat = int(float(minlat))
            maxlat = int(float(maxlat))
        if info[0] == 'S':
            minlat = int(float(minlat)) * -1
            maxlat = int(float(maxlat)) * -1
        if info[5] == 'E':
            minlon = int(float(minlon))
            maxlon = int(float(maxlon))
        if info[5] == 'W':
            minlon = int(float(minlon)) * -1
            maxlon = int(float(maxlon)) * -1
        return minlon,minlat,maxlon+1-celllonlat,maxlat+1-celllonlat
def get_dem_name(lon,lat):
    if lat >= 0:
        latNS = 'N'
    else:
        latNS = 'S'
    if lon >= 0:
        lonEW = 'E'
    else:
        lonEW = 'W'
    Name = "{:s}{:02d}{:s}{:03d}_FABDEM_V1-0.tif".format(latNS,abs(lat),lonEW,abs(lon))
    # print(Name)
    return Name
def get_merge_name(extends):
    """
    给合并的dem,sink起名字。是外包矩形边上的左下角经纬度
    :param extends:
    :return:
    """
    minLon = min(extends[0],extends[2])
    minLat = min(extends[1],extends[3])

    maxLon = max(extends[0], extends[2])
    maxLat = max(extends[1], extends[3])

    if minLat>=0:
        latNS = 'N'
    else:
        latNS = 'S'
    if minLon>=0:
        lonEW = 'E'
    else:
        lonEW = 'W'
    Name = "{:s}{:02d}{:02d}{:s}{:03d}{:03d}_FABDEM_V1-0.tif".format(latNS, abs(minLat),abs(maxLat),lonEW, abs(minLon),abs(maxLon))
    return Name
def search_max_rec(extents):
    """
    按覆盖面积从小到大的顺序的list，用来寻找互不包含的外包矩形list
    :param extents:
    :return:
    """
    extentsCopy = extents.copy()
    result = []
    while extentsCopy:
        popE = extentsCopy.pop()
        result.append(popE)
        extents = extentsCopy.copy()
        for info in extents:
            if popE[0] <= info[0] <= popE[2] and popE[0] <= info[2] <= popE[2] and popE[1] <= info[1] <= popE[3] and popE[1] <= info[3] <= popE[3]:
                extentsCopy.remove(info)

    return result






