# -*- coding: utf-8 -*-
"""
@Time ： 2024/6/4 9:19
@Auth ：
@File ：sink.py
@IDE ：PyCharm
"""
import math
import os.path
import shutil

import matplotlib.pyplot as plt
import Raster
import numpy as np
from osgeo import gdal,ogr
#import WBT.whitebox_tools
# import whitebox_tools.whitebox-tools-master.whitebox_tools
import whitebox
from shapely.geometry import LineString, Point

import valid

wbt = whitebox.WhiteboxTools()
import geopandas as gpd
wbt.exe_path = '/datanode05/zhangbin/.conda/envs/zhangbin/lib/python3.8/site-packages/whitebox' # Linux必加，否则无法运行

dmove=[(0,1),(1,1),(1,0),(1,-1),(0,-1),(-1,-1),(-1,0),(-1,1)]
dmove_dic = {1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1), 16: (0, -1), 32: (-1, -1), 64: (-1, 0), 128: (-1, 1)}


import scipy.ndimage as ndimage


def fill_sinks(dem_file):

    dem = Raster.get_raster(dem_file)
    proj,geo,nodata = Raster.get_proj_geo_nodata(dem_file)
    # 获取高程数据的大小
    rows, cols = dem.shape

    # 用于存储填洼后的 DEM 数据
    filled_dem = dem.copy()

    # 初始化一个标记数组，记录已经填过的像素
    filled = np.zeros_like(dem, dtype=bool)

    # 填洼过程：逐像素检查
    for row in range(1, rows - 1):
        for col in range(1, cols - 1):
            if dem[row,col] == nodata:
                continue
            if not filled[row, col]:  # 如果该位置未填洼
                # 获取该位置周围的邻域（3x3窗口）
                window = dem[row - 1:row + 2, col - 1:col + 2]

                # 找到该邻域中的最大值（模拟填洼）
                max_value = np.max(window)

                # 如果当前点的值低于周围邻域的最大值，则填洼
                if dem[row, col] < max_value:
                    filled_dem[row, col] = max_value
                    filled[row, col] = True  # 标记为已填洼

    return filled_dem

def Cal_Filldepression(DEM_file,outDEMfile):
    """
    临时代码，用于洼地太大超出内存限制，选择用填洼后的DEM-原始DEM获取sink，同时限制填洼深度
    :param DEM_file: 
    :param outDEMfile: 
    :return: 
    """
    wbt.fill_depressions(DEM_file,outDEMfile,fix_flats=True,max_depth=11)

    
def Cal_filled(DEM_file,FilledDem_file):
    wbt.fill_depressions_planchon_and_darboux(DEM_file,FilledDem_file,fix_flats=True)


def Cal_dir(DEM_file,Dir_file):
    wbt.d8_pointer(DEM_file,Dir_file,esri_pntr=True)

def Cal_sink(DEM_file,sink_file):
    wbt.sink(DEM_file,sink_file)

def Cal_acc(Dir_file,acc_file):
    wbt.d8_flow_accumulation(Dir_file,acc_file,esri_pntr=True)



def calculate_fdir_depression(DEM_venu):
    """
    根据run_data下Asia的每个文件夹来计算fdir和sink
    :param DEM_venu:
    :return:
    """
    if len(DEM_venu.split(".")) != 1:
        return
    baseName = os.path.basename(DEM_venu)
    DEM_file = os.path.join(DEM_venu,baseName+'_burnDEM.tif')
    fdir_file = os.path.join(DEM_venu,baseName+'_burndir.tif')
    sink_file = os.path.join(DEM_venu,baseName+'_burnsink.tif')

    Cal_dir(DEM_file, fdir_file)
    if not os.path.exists(fdir_file):
        Cal_dir(DEM_file,fdir_file)
        pass
    else:
        # pass
        return
    if not os.path.exists(sink_file):
        Cal_sink(DEM_file,sink_file)

    # wbt.fill_depressions_planchon_and_darboux(os.path.join(DEM_venu,baseName+'_dem.tif'),os.path.join(DEM_venu,baseName+'_FillDEM.tif'))
    # A = Raster.get_raster(os.path.join(DEM_venu,baseName+'_FillDEM.tif'))
    # proj,geo,nodata = Raster.get_proj_geo_nodata(os.path.join(DEM_venu,baseName+'_FillDEM.tif'))
    #
    # Raster.save_raster(os.path.join(DEM_venu,baseName+'_FillDEM.tif'),A,proj,geo,gdal.GDT_Int16,nodata)
    #
    # wbt.d8_pointer(os.path.join(DEM_venu,baseName+'_FillDEM.tif'),fdir_file,esri_pntr=True)



def Cal_Stream(Acc_file,Stream_file,thresold):
    Acc=Raster.get_raster(Acc_file)
    proj,geo,nodata=Raster.get_proj_geo_nodata(Acc_file)
    row,col=Acc.shape
    Stream=np.zeros((row,col),dtype = np.int8)
    Stream[Acc>=thresold]=1

    Raster.save_raster(Stream_file,Stream,proj,geo,gdal.GDT_Byte,0)



def streamFeature(dir_file,stream_file,outfile):
    # wbt.raster_streams_to_vector(stream_file,dir_file,outfile,esri_pntr=True)
    # wbt.raster_to_vector_lines(stream_file,outfile)
    wbt.stream_link_identifier(dir_file,stream_file,outfile,esri_pntr=True)

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

def Check_extent(row,col,x,y):
    if 0<=x<row and 0<=y<col:
        return True
    return False

def bfs(id,i,j,arr,vis):
    """
    搜索块集合，并标记vis=1
    :param id: 当前搜索集合的id
    :param i:
    :param j:
    :param arr:
    :param vis:
    :return:
    """

    row,col=arr.shape
    pop_list=[(i,j)]
    collection_list=[(i,j)]
    vis[i,j]=1
    while pop_list:
        pop_cell=pop_list.pop()
        for k in range(8):
            now_cell_x=pop_cell[0]+dmove[k][0]
            now_cell_y=pop_cell[1]+dmove[k][1]
            if Check_extent(row,col,now_cell_x,now_cell_y):
                if arr[now_cell_x,now_cell_y]==id:
                    if vis[now_cell_x,now_cell_y]==0:
                        pop_list.insert(0,(now_cell_x,now_cell_y))
                        collection_list.insert(0,(now_cell_x,now_cell_y))
                        vis[now_cell_x,now_cell_y]=1

    return vis,collection_list,len(collection_list)

def Get_extent_out(DEM,sink,DEM_nodata,sink_cells,row,col,EXTENT):

    extent_cells=[]
    for cell in sink_cells:
        for k in range(8):
            next_cell=(cell[0]+dmove[k][0],cell[1]+dmove[k][1])
            if Check_extent(row,col,next_cell[0],next_cell[1]):
                if sink[next_cell[0],next_cell[1]]!=sink[cell[0],cell[1]] and DEM[next_cell[0],next_cell[1]]!=DEM_nodata:
                    # 会存在无效值，后面使用时，要注意判断无效值
                    EXTENT[next_cell[0],next_cell[1]]=sink[cell[0],cell[1]]
                    extent_cells.append((next_cell[0],next_cell[1],DEM[next_cell[0],next_cell[1]]))

    return EXTENT,extent_cells

def Get_OUT(Dir,DEM,sink_id,oridin_sink,sink_nodata,extent_cells,sink_cells,DEM_nodata,OUTLET,row,col,size=30):
    Volume = 0
    if len(extent_cells)==0:
        return OUTLET,Volume

    outlet=[-1,-1,math.inf]
    for cell in extent_cells:
        if cell[2]!=DEM_nodata and cell[2]<=outlet[2]:
            now_dir=Dir[cell[0],cell[1]]
            if now_dir in dmove_dic:
                next_cell=(cell[0]+dmove_dic[now_dir][0],cell[1]+dmove_dic[now_dir][1])
                if Check_extent(row,col,next_cell[0],next_cell[1]):
                    if oridin_sink[next_cell[0],next_cell[1]]==sink_id:
                        continue
            outlet = cell
            # if oridin_sink[cell[0],cell[1]]==sink_nodata:
            #     # print(oridin_sink[cell[0], cell[1]])
            #     outlet=cell
    # 计算最大蓄水容量，只计算低于出水口高程的像元体积
    Max_H=outlet[2]

    for cell in sink_cells:
        if DEM[cell[0],cell[1]]<=Max_H:
            Volume+=size*size*(Max_H-DEM[cell[0],cell[1]])
    OUTLET[outlet[0],outlet[1]]=sink_id
    return outlet,Volume,OUTLET

def Correct_sink(Dir,outlet,origin_sink,sink_id,Vis_correct,row,col):

    pop_cells=[(outlet[0],outlet[1])]

    while pop_cells:
        pop_cell=pop_cells.pop()
        for k in range(8):
            next_cell=(pop_cell[0]+dmove[k][0],pop_cell[1]+dmove[k][1])
            if Check_extent(row,col,next_cell[0],next_cell[1]):
                if origin_sink[next_cell[0],next_cell[1]]==sink_id :
                    if Vis_correct[next_cell[0],next_cell[1]]==0:
                        Dir[next_cell[0], next_cell[1]] = 2 ** ((k + 4) % 8)
                        Vis_correct[next_cell[0],next_cell[1]]=1
                        pop_cells.insert(0,next_cell)
    return Dir

def Correct_delete_sink(Dir,outlet,origin_sink,sink,sink_id,sink_nodata,Vis_correct,row,col):

    pop_cells=[(outlet[0],outlet[1])]

    while pop_cells:
        pop_cell=pop_cells.pop()
        for k in range(8):
            next_cell=(pop_cell[0]+dmove[k][0],pop_cell[1]+dmove[k][1])
            if Check_extent(row,col,next_cell[0],next_cell[1]):
                if origin_sink[next_cell[0],next_cell[1]]==sink_id :
                    if Vis_correct[next_cell[0],next_cell[1]]==0:
                        Dir[next_cell[0], next_cell[1]] = 2 ** ((k + 4) % 8)
                        Vis_correct[next_cell[0],next_cell[1]]=1
                        pop_cells.insert(0,next_cell)
                        sink[next_cell[0],next_cell[1]]=sink_nodata
    return Dir,sink

def Process_sink(sink_file,Dir_file,DEM_file,Out_venu,Volume_thre):
    sink=Raster.get_raster(sink_file)
    Dir=Raster.get_raster(Dir_file)
    DEM=Raster.get_raster(DEM_file)
    proj,geo,sink_nodata=Raster.get_proj_geo_nodata(sink_file)
    _,_,DEM_nodata=Raster.get_proj_geo_nodata(DEM_file)
    row,col=sink.shape
    Vis = np.zeros((row, col))
    EXTENT=np.zeros((row,col))
    EXTENT[:,:]=-9999
    OUTLET = np.zeros((row, col))
    OUTLET[:, :] = -9999
    Vis_correct=np.zeros((row,col))
    origin_sink=sink.copy()
    for i in range(row):
        for j in range(col):
            if sink[i,j]!=sink_nodata and Vis[i,j]==0:
                # 还需要考虑湖泊（F1添加湖泊指示；F2认为是大于1km2的为湖泊）
                sink_id=sink[i,j]
                Vis,sink_cells,sink_cells_len=bfs(sink[i,j],i,j,sink,Vis)  # 搜索同一片sink
                EXTENT,temp_extent_cells=Get_extent_out(DEM,sink,DEM_nodata,sink_cells,row,col,EXTENT)   # 搜索sink的外围边界，最后输出会有重叠，所以对每个sink单独处理
                temp_outlet,temp_Volume,OUTLET=Get_OUT(Dir,DEM,sink_id,origin_sink,sink_nodata,temp_extent_cells,sink_cells,DEM_nodata,OUTLET,row,col)
                print(sink_id,len(sink_cells)*28.33*28.33/1000/1000,temp_Volume)
                if temp_Volume>=Volume_thre:
                    # 修正流向
                    Dir = Correct_sink(Dir, temp_outlet, origin_sink, sink_id, Vis_correct, row, col)

                    pass
                else:
                    # 将sink填洼，并流出
                    Dir,sink = Correct_delete_sink(Dir, temp_outlet, origin_sink,sink, sink_id, sink_nodata,Vis_correct, row, col)

                    pass
    # Out=os.path.join(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\RESULT\20240818','New_Dir'+str(Volume_thre)+'.tif')
    Out = os.path.join(Out_venu,
                       'New_Dir' + str(Volume_thre) + '.tif')
    Raster.save_raster(Out,Dir,proj,geo,gdal.GDT_Byte,255)
    print(Volume_thre,len(np.unique(sink)))
    # out=os.path.join(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\RESULT\20240818','New_sink'+str(Volume_thre)+'.tif')
    out = os.path.join(Out_venu,
                       'New_sink' + str(Volume_thre) + '.tif')
    # print(out)
    Raster.save_raster(out,sink,proj,geo,gdal.GDT_Float32,sink_nodata)
    # Acc_file = os.path.join(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\RESULT\20240818','New_Acc'+str(Volume_thre)+'.tif')
    # Stream_file = os.path.join(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\RESULT\20240818','New_Stream'+str(Volume_thre)+'.tif')
    Acc_file = os.path.join(Out_venu,
                            'New_Acc' + str(Volume_thre) + '.tif')
    Stream_file = os.path.join(Out_venu,
                               'New_Stream' + str(Volume_thre) + '.tif')
    Cal_acc(Out, Acc_file)
    Cal_Stream(Acc_file, Stream_file, 800)
    return len(np.unique(sink))

def Judge_key_point_type(i,j,Dir,stream,stream_nodata,sink,sink_nodata,row,col):
    """
    判断河流与sink的交点类型：计算待判断点的下游，若为sink，并且判断上游是否为sink_nodata，若不是，则该点为河流汇入sink的点，编码为2；若为非sink，则为sink的出口点，标记为1.
    :param i:
    :param j:
    :param Dir:
    :param sink:
    :param sink_nodata:
    :param row:
    :param col:
    :return:
    """
    now_dir=Dir[i,j]
    if now_dir in dmove_dic:
        next_cell=(i+dmove_dic[now_dir][0],j+dmove_dic[now_dir][1])
        if Check_extent(row,col,next_cell[0],next_cell[1]):
            if sink[next_cell[0],next_cell[1]]!=sink_nodata:
                # 判断上游是否为sink_nodata，，若不是，则为河流汇入sink点
                up_cells=get_rever_D8(Dir,i,j)
                for cell in up_cells:
                    if sink[cell[0],cell[1]]==sink_nodata and stream[cell[0],cell[1]]!=stream_nodata:
                        return 2# 河流汇入sink点
                return 0

            else:
                return 1 # sink出流点
    else:
        return 0 # 无效点

def search_key_point(sink_file,stream_file,Dir_file):
    sink=Raster.get_raster(sink_file)
    stream=Raster.get_raster(stream_file)
    Dir=Raster.get_raster(Dir_file)

    proj,geo,sink_nodata=Raster.get_proj_geo_nodata(sink_file)
    _,_,stream_nodata=Raster.get_proj_geo_nodata(stream_file)
    row,col=sink.shape


    Result=np.zeros((row,col))
    for i in range(row):
        for j in range(col):
            if sink[i,j]!=sink_nodata and stream[i,j]!=stream_nodata:
                # 找到交点，判断是类型
                Result[i,j]=Judge_key_point_type(i,j,Dir,stream,stream_nodata,sink,sink_nodata,row,col)
                if Result[i,j]==0:
                    stream[i,j]=stream_nodata

    Raster.save_raster(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\RESULT\points.tif',Result,proj,geo,gdal.GDT_Byte,0)
    Raster.save_raster(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\RESULT\New_stream18000.tif',stream,proj,geo,gdal.GDT_Byte,stream_nodata)


def stream_class(stream_file,fdr_file):
    """
    给河流分级并转成矢量
    :param stream_file:
    :param fdr_file:
    :return:
    """
    baseDir = os.path.dirname(stream_file)
    #
    # stream_link = os.path.join(baseDir,'stream_link.tif')
    out_stream = os.path.join(baseDir,'stream_order.tif')
    out_shp = os.path.join(baseDir,'stream_clsss.shp')
    #
    # # wbt.stream_link_identifier(fdr_file,stream_file,stream_link,esri_pntr=True)
    #
    # # wbt.stream_link_class(fdr_file,stream_file,out_stream,esri_pntr=True)
    # wbt.strahler_stream_order(fdr_file,stream_link,out_stream,esri_pntr=True)
    #
    # streamFeature(fdr_file,out_stream,out_shp)

    import numpy as np

    # ----------------------------------------
    # 1. 定义 ESRI D8 方向编码（可按需修改）
    # ----------------------------------------
    D8_DIR = {
        1: (0, 1),  # East
        2: (1, 1),  # Southeast
        4: (1, 0),  # South
        8: (1, -1),  # Southwest
        16: (0, -1),  # West
        32: (-1, -1),  # Northwest
        64: (-1, 0),  # North
        128: (-1, 1)  # Northeast
    }

    # ----------------------------------------
    # 2. 构建所有河流像元的 upstream（上游列表）
    # ----------------------------------------
    def build_upstream(flowdir, river):
        rows, cols = flowdir.shape
        river_idx = np.argwhere(river > 0)

        # 初始化 upstream 字典，所有河流像元都要有 Key
        upstream = {(int(r), int(c)): [] for r, c in river_idx}

        # 遍历每个河流像元，查它的下游
        for r, c in river_idx:
            fd = int(flowdir[r, c])
            if fd not in D8_DIR:
                continue

            dr, dc = D8_DIR[fd]
            rr, cc = r + dr, c + dc

            # 如果下游也是河流，则下游的 upstream 中加入当前像元
            if 0 <= rr < rows and 0 <= cc < cols and river[rr, cc] > 0:
                upstream[(int(rr), int(cc))].append((int(r), int(c)))

        return upstream

    # ----------------------------------------
    # 3. 计算 Strahler 河流等级
    # ----------------------------------------
    def compute_strahler(flowdir, river,geo):
        upstream = build_upstream(flowdir, river)

        # 初始化所有河流像元的等级为 0
        order = {cell: 0 for cell in upstream.keys()}

        # 拓扑排序起点：所有无上游的源头
        queue = [cell for cell, ups in upstream.items() if len(ups) == 0]

        while queue:
            r, c = queue.pop()
            ups = upstream[(r, c)]

            # 1) 源头等级 = 1
            if len(ups) == 0:
                order[(r, c)] = 1
            else:
                up_orders = [order[u] for u in ups]
                m = max(up_orders)

                # Strahler 规则
                if up_orders.count(m) >= 2:
                    order[(r, c)] = m + 1
                else:
                    order[(r, c)] = m

            # 处理下游像元
            fd = int(flowdir[r, c])
            if fd in D8_DIR:
                dr, dc = D8_DIR[fd]
                rr, cc = r + dr, c + dc

                if (rr, cc) in upstream:
                    # 如果下游所有上游都已赋值，则加入队列
                    if all(order[u] > 0 for u in upstream[(rr, cc)]):
                        queue.append((rr, cc))


        # --------------------------- 栅格转矢量的代码 ------------------------------
        # 拓扑排序起点：所有无上游的源头
        queue = [cell for cell, ups in upstream.items() if len(ups) == 0]
        # 初始化所有河流像元的状态为 0
        vis = {cell: 0 for cell in upstream.keys()}
        lines = {}
        rid = 1
        stream_order = {}
        while queue:
            r, c = queue.pop()
            ups = upstream[(r, c)]

            # 1) 非源头跳过
            if len(ups) != 0:
                continue
            pop_cells = [(r,c)]
            while pop_cells:
                rr,cc = pop_cells.pop()
                lines.setdefault(rid, []).append(valid.pixel_to_lonlat(geo,rr,cc))
                stream_order.setdefault(rid,order[(rr,cc)])
                stream_order[rid] = max(stream_order[rid],order[(rr,cc)])
                vis[(rr,cc)] = 1
                fd = int(flowdir[rr, cc])

                if fd in D8_DIR:
                    dr, dc = D8_DIR[fd]
                    rrr, ccc = rr + dr, cc + dc
                    if (rrr,ccc) not in vis:
                        # 非河流，pass
                        continue
                    if vis[(rrr,ccc)] == 1:
                        # 已经查过，pass
                        continue
                    # 没查过，就继续追溯
                    if order[(rrr,ccc)] != order[rr,cc]:
                        # 新河流
                        rid += 1

                    pop_cells.append((rrr,ccc))
            rid += 1
            # break

        # 生成 LineString
        results = []
        for rid in lines:
            if len(lines[rid])>1:
                results.append({
                    "river_id": rid,
                    "geometry": LineString(lines[rid]),
                    "order":float(stream_order[rid])
                })



        # --------------------------- 栅格转矢量的代码 ------------------------------
        return order, gpd.GeoDataFrame(results, crs="EPSG:4326")

    # ----------------------------------------
    # 4. 将字典形式等级写回栅格
    # ----------------------------------------
    def order_to_raster(order, shape):
        out = np.zeros(shape, dtype=np.int16)

        for (r, c), val in order.items():
            out[r, c] = val

        return out

    # flowdir: numpy array，D8流向
    # river:   numpy array，河流掩膜（1=河流，0=非河流）

    flowdir = Raster.get_raster(fdr_file)


    river = Raster.get_raster(stream_file)
    proj,geo,nodata = Raster.get_proj_geo_nodata(stream_file)

    order_dict,gdf = compute_strahler(flowdir, river,geo)

    # gdf = raster_to_vector_rivers(flowdir, river_id, transform)
    gdf.to_file(out_shp)


    # order_raster = order_to_raster(order_dict, flowdir.shape)
    # Raster.save_raster(out_stream,order_raster,proj,geo,gdal.GDT_Int16,0)

    # streamFeature(fdr_file, out_stream, out_shp)
    # import stream2featutre
    # stream2featutre.stream_raster2feature(fdr_file,out_stream,out_shp)


def zip(venu):
    """
    压缩生成的region和new_acc，
    :return:
    """

    fdr_venu = os.path.join(venu,'region')
    acc_venu = os.path.join(venu,'new_acc')

    out_zip_fdr = os.path.join(venu,'flow_firection')
    out_zip_acc = os.path.join(venu,'flow_accumulation')
    # folder_path = r"D:/data/my_folder"  # 要压缩的文件夹
    # zip_path = r"D:/data/my_folder"  # 输出 zip（不要加 .zip）

    shutil.make_archive(
        base_name=out_zip_fdr,
        format="zip",
        root_dir=fdr_venu
    )

    shutil.make_archive(
        base_name=out_zip_acc,
        format="zip",
        root_dir=acc_venu
    )


if __name__=='__main__':

    # Get_extent_sink(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\sink.tif',r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\Dir.tif',
    #                 r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\DEM.tif',100)

    # Process_sink(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\sink.tif',
    #                 r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\Dir.tif',
    #                 r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\DEM.tif', 2000)
    # Process_sink(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\sink.tif',
    #                 r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\Dir.tif',
    #                 r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\DEM.tif', 18000)
    # Process_sink(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\sink.tif',
    #                 r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\Dir.tif',
    #                 r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\DEM.tif', 93000)
    # search_key_point(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\RESULT\New_sink18000.tif',
    #                  r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\RESULT\New_Stream18000.tif',
    #                  r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\RESULT\New_Dir18000.tif')

    # Procedure
    # DEM_file = '//datanode05/zhangbin/TBP_Stream/TBP/TBP_FABDEM.tif'
    # Dir_file = '/datanode05/zhangbin/TBP_Stream/DATA/20240818/TBP_Dir.tif'
    # sink_file = '/datanode05/zhangbin/TBP_Stream/DATA/20240818/TBP_sink.tif'
    DEM_file = r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\DEM.tif'
    Dir_file = r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\Dir.tif'
    sink_file = r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\sink1.tif'
    # 1、使用未填洼的DEM计算流向和sink
    # Cal_dir(DEM_file,Dir_file)
    # Cal_sink(DEM_file,sink_file)
    # # 2、处理洼地
    Process_sink(sink_file,
                    Dir_file,
                    DEM_file,'/datanode05/zhangbin/TBP_Stream/DATA/20240818', 93001)
    Process_sink(sink_file,
                 Dir_file,
                 DEM_file, r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604', 93001)



    pass