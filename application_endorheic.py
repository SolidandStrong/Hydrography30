# -*- coding: utf-8 -*-
"""
@Time ： 2025/4/2 14:36
@Auth ：
@File ：application_endorheic.py
@IDE ：PyCharm

处理内流区的nodata问题及内流区的相关应用函数, nodata = 0
"""
import csv
import os.path
import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import Geod

import Raster
from osgeo import gdal,ogr,osr
import split_cal_acc
import genral_functions
import sink
import valid
from genral_functions import *

def get_endorheic_final_point(dir_file,out_point_file,out_point_file_tif):
    """
    识别dir中的内流终点：8邻域内都不为空值的nodata像元
    :param dir_file:流向
    :param out_point_file: 矢量点文件
    :param out_point_file_tif: 栅格文件
    :return:
    """

    fdir = Raster.get_raster(dir_file)
    proj,geo,nodate = Raster.get_proj_geo_nodata(dir_file)
    print(nodate)

    row,col = fdir.shape

    # 2. 创建一个输出的 Shapefile 矢量文件
    # 3. 创建坐标系对象
    spatial_ref = osr.SpatialReference()
    spatial_ref.ImportFromWkt(proj)  # 使用栅格的投影信息
    shapefile_driver = ogr.GetDriverByName("ESRI Shapefile")
    out_shapefile = shapefile_driver.CreateDataSource(out_point_file)


    # 3. 创建一个点的图层
    layer = out_shapefile.CreateLayer("points", spatial_ref,geom_type=ogr.wkbPoint)


    result = np.zeros((row,col),dtype = np.int8)

    # for i in range(row):
    #     for j in range(col):
    #         if fdir[i,j] == nodate:
    #             continue
    #         flag = False
    #         if fdir[i,j] == 0:
    #             flag = True
    #         # mask = fdir[i-1:i+2,j-1:j+2]
    #         # if np.all(mask == nodate) == 1:
    #         #     continue
    #         # flag = True
    #         # for k in range(8):
    #         #     next_cell = (i + dmove[k][0], j + dmove[k][1])
    #         #     if not check_boundary(next_cell[0], next_cell[1], row, col):
    #         #         flag = False
    #         #         break
    #         #     if fdir[next_cell[0], next_cell[1]] == nodate:
    #         #         flag = False
    #         #         break
    #         if flag:
    #             x_coord = geo[0] + j * geo[1] + j * geo[2]
    #             y_coord = geo[3] + i * geo[4] + i * geo[5]
    #
    #             # 创建点几何对象
    #             point = ogr.Geometry(ogr.wkbPoint)
    #             point.AddPoint(x_coord, y_coord)
    #
    #             # 创建属性字段（可选）
    #             feature = ogr.Feature(layer.GetLayerDefn())
    #             feature.SetGeometry(point)
    #             layer.CreateFeature(feature)
    #
    #             result[i, j] = 1

    # 2. 找出所有位置 (行列坐标)
    positions = np.argwhere(fdir == 0)

    print("找到的位置数量:", len(positions))
    for i in positions:
        x_coord = geo[0] + i[1] * geo[1] + i[1] * geo[2]
        y_coord = geo[3] + i[0] * geo[4] + i[0] * geo[5]

        # 创建点几何对象
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(x_coord, y_coord)

        # 创建属性字段（可选）
        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetGeometry(point)
        layer.CreateFeature(feature)

        result[i[0], i[1]] = 1


    out_shapefile = None
    print("栅格转点矢量完成！")
    Raster.save_raster(out_point_file_tif,result,proj,geo,gdal.GDT_Byte,0)

def get_endorheic_final_point_csv(dir_file,out_csv_file):
    """
    识别dir中的内流终点：fdr = 0
    :param dir_file:流向
    :param out_point_file: 矢量点文件
    :param out_point_file_tif: 栅格文件
    :return:
    """

    fdir = Raster.get_raster(dir_file)
    proj, geo, nodate = Raster.get_proj_geo_nodata(dir_file)
    positions = np.argwhere(fdir == 0)
    result = []

    print("找到的位置数量:", len(positions))
    for i in positions:
        x_coord = geo[0] + i[1] * geo[1] + i[1] * geo[2]
        y_coord = geo[3] + i[0] * geo[4] + i[0] * geo[5]

        result.append([x_coord,y_coord])

    df = pd.DataFrame(result)
    df.columns = ['Lon','Lat']
    df.to_csv(out_csv_file)

def po_get_endorheic_final_point(dir_file,outvenu):
    if not os.path.exists(outvenu):
        os.mkdir(outvenu)
    region = os.path.join(outvenu,'region')
    split_cal_acc.split_raster(dir_file, region, block_height=10000)

    out_point_lonlat_venu = os.path.join(outvenu,'endorheic_point')
    if not os.path.exists(out_point_lonlat_venu):
        os.mkdir(out_point_lonlat_venu)

    paras = []
    for file in os.listdir(region):
        paras.append([os.path.join(region,file),os.path.join(out_point_lonlat_venu,file.split('.')[0]+'.csv')])

    Po = Pool(10)
    for para in paras:
        Po.apply_async(get_endorheic_final_point_csv,(para[0],para[1],))

    Po.close()
    Po.join()

    merge = []
    for file in os.listdir(out_point_lonlat_venu):

        df = pd.read_csv(os.path.join(out_point_lonlat_venu,file))
        merge.append(df)
    df = pd.concat(merge, ignore_index=True)
    df.to_csv(os.path.join(outvenu,'endorheic_point.csv'))


def extract_endorheic_watershed(fdr_file,endorheic_csv,out_venu):
    """
    根据endorheic_csv中的内流终点提取内流流域
    可能需要对点校正，检查像元8-邻域内，找到fdr=0的像元再提取上游
    :param fdr_file:
    :param endorheic_csv:
    :return:
    """
    if not os.path.exists(out_venu):
        os.mkdir(out_venu)

    fdr = Raster.get_raster(fdr_file)
    proj,geo,f_nodata = Raster.get_proj_geo_nodata(fdr_file)
    row_num,col_num = fdr.shape


    lon_lats = []
    # with open(endorheic_csv,'r') as f:
    #     # (118.5981249999991, 0.0002777777777777746, 0.0, 32.340694444444416, 0.0, -0.0002777777777777746)
    #     reader = csv.reader(f)
    #     n = 0
    #     for i in reader:
    #         if n == 0:
    #             n += 1
    #             continue
    #
    #         col = int((float(i[2]) - geo[0]) / geo[1])
    #         row = int((geo[3] - float(i[3])) / abs(geo[5]))
    #
    #         lon_lats.append([row, col,os.path.join(out_venu,"Endorheic_basin_"+str(i[0])+'.tif')])
    #
    #     f.close()

    # ------------------------------- 提取里海内流区的部分-----------------------------------------------
    lihai = (49.1957467,41.0520217)
    col = int((float(lihai[0]) - geo[0]) / geo[1])
    row = int((geo[3] - float(lihai[1])) / abs(geo[5]))
    lon_lats.append([row,col,"/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/Europe_ACC_final/Endoheic_application/里海/Lihai_Basin.tif"])

    # lihai = (34.016351,60.96088)
    # col = int((float(lihai[0]) - geo[0]) / geo[1])
    # row = int((geo[3] - float(lihai[1])) / abs(geo[5]))
    # lon_lats.append([row, col,
    #                  "/datanode05/zhangbin/FAB_hydrography/Final_results/Fdr/Europe/修正缺失的海/2/modified_dir_2.tif"])
    # ------------------------------- 提取里海内流区的部分-----------------------------------------------

    for lon_lat in lon_lats:

        row = lon_lat[0]
        col = lon_lat[1]

        try:
            # if fdr[row,col] != 0:
            #     for k in range(8):
            #         next_cell = (row+dmove[k][0],col+dmove[k][1])
            #         if fdr[next_cell[0],next_cell[1]] == 0:
            #             row = next_cell[0]
            #             col = next_cell[1]
            #             break

            vis = np.zeros((row_num, col_num), dtype=np.int8)
            max_row = 0
            min_row = row_num
            max_col = 0
            min_col = col_num
            # row1, col1 = lonlat_to_pixel(geo, lon_lat[0], lon_lat[1])
            row1 = row
            col1 = col
            # print(row1,col1)
            popCells = [(row1, col1)]
            basinCells = [(row1, col1)]

            # ------------------------------- 提取里海内流区的部分-----------------------------------------------
            while popCells:
                pop_cell = popCells.pop()
                for k in range(8):
                    next_cell = (pop_cell[0]+dmove[k][0],pop_cell[1]+dmove[k][1])
                    if vis[next_cell[0],next_cell[1]] == 1:
                        continue
                    flag = False
                    for kk in range(8):
                        nnext_cell = (next_cell[0]+dmove[kk][0],next_cell[1]+dmove[kk][1])
                        if fdr[nnext_cell[0],nnext_cell[1]] == f_nodata:
                            flag = True
                            break
                    if flag:

                        popCells.append(next_cell)
                        basinCells.append(next_cell)
                        vis[next_cell[0],next_cell[1]] = 1
            popCells = basinCells.copy()
            # print(popCells)
            # ------------------------------- 提取里海内流区的部分-----------------------------------------------

            # if vis[row1, col1] == 1:
            #     break
            # vis[row1, col1] = 1

            while popCells:
                popCell = popCells.pop()
                max_row = max(max_row, popCell[0])
                max_col = max(max_col, popCell[1])
                min_row = min(min_row, popCell[0])
                min_col = min(min_col, popCell[1])
                # if vis[popCell[0], popCell[1]] == 1:
                #     continue
                vis[popCell[0], popCell[1]] = 1

                upCells = valid.get_rever_D8(fdr, popCell[0], popCell[1])

                # basinCells += upCells
                # popCells += upCells

                for upCell in upCells:
                    if vis[upCell[0], upCell[1]] == 1:
                        continue
                    # basinCells.append(upCell)
                    popCells.append(upCell)


            # # 坐标转换至小坐标系
            new_row = max_row - min_row + 1
            new_col = max_col - min_col + 1
            basin = vis[min_row:max_row+1,min_col:max_col+1].copy()
            # basin = np.zeros((new_row, new_col), dtype=np.int8)
            # for basinCell in basinCells:
            #     basin[basinCell[0] - min_row, basinCell[1] - min_col] = 1
            #
            new_lon, new_lat = valid.pixel_to_lonlat(geo, min_row, min_col)
            new_geo = (new_lon, geo[1], geo[2], new_lat, geo[4], geo[5])
            #
            Raster.save_raster(lon_lat[2], basin, proj, new_geo, gdal.GDT_Byte, 0)
        except Exception as e:
            print(lon_lat, e)

def sbatch_raster_endorheic_venu(raster_venu,shp_venu):
    """
    批量将文件夹内的栅格矢量化
    :param raster_venu:
    :param shp_venu:
    :return:
    """
    if not os.path.exists(shp_venu):
        os.mkdir(shp_venu)

    paras = []
    for name in os.listdir(raster_venu):
        if name.split('.')[-1] != 'tif':
            continue
        paras.append([os.path.join(raster_venu,name),os.path.join(shp_venu,name.split('.')[0]+'.shp')])

    po = Pool(10)
    for para in paras:
        po.apply_async(valid.raster2Feature,(para[0],para[1],))

    po.close()
    po.join()

def statis_endorheic_information(csv_file,shp_venu):
    """
    统计内流流域的信息
    :param shp_venu:
    :return:
    """

    baseDir = os.path.dirname(shp_venu)

    paras = [['ID','Lon','Lat','Area']]
    with open(csv_file,'r') as f:

        reader = csv.reader(f)
        n = 0
        for i in reader:
            if n==0:
                n += 1
                continue
            vector_file = os.path.join(shp_venu, "Endorheic_basin_" + str(i[0]) + '.shp')
            temp = [i[0],float(i[2]),float(i[3])]
            if os.path.exists(vector_file):
                # 读取矢量数据
                gdf = gpd.read_file(vector_file)

                # 定义 WGS84 地理坐标系
                geod = Geod(ellps="WGS84")

                # 计算每个多边形的面积（单位：平方米）
                gdf['area_m2'] = gdf.geometry.apply(lambda geom: abs(geod.geometry_area_perimeter(geom)[0]))

                # 以平方千米为单位
                gdf['area_km2'] = gdf['area_m2'] / 1e6

                temp.append(gdf['area_km2'].sum())

                paras.append(temp)

    with open(os.path.join(baseDir,'Endorheic_basin_information.csv'),'w') as f:
        writer = csv.writer(f)
        writer.writerows(paras)
        f.close()


def merge_endorheic_basin_to_gpkg(venu):
    """
    把目录下的shp融合并输出到一个gpkg下，并计算面积
    :param venu:
    :return:
    """

    import geopandas as gpd
    import pandas as pd
    import os

    # 输入文件夹
    folder = venu

    # 输出 GPKG 文件
    out_gpkg = os.path.join(os.path.dirname(venu),"endorheic_basins.gpkg")
    layer_name = "endorheic"

    records = []

    for file in os.listdir(folder):
        if file.lower().endswith(".shp"):
            shp_path = os.path.join(folder, file)

            name = file.split('.')[0].split('_')[-1]#os.path.splitext(file)[0]  # 作为ID

            gdf = gpd.read_file(shp_path)

            gdf = gdf.set_crs("EPSG:4326")

            # 融合成一个单一 Geometry
            merged_geom = gdf.unary_union

            # 定义 WGS84 地理坐标系
            geod = Geod(ellps="WGS84")
            # 计算每个多边形的面积（单位：平方米）
            gdf['area_m2'] = gdf.geometry.apply(lambda geom: abs(geod.geometry_area_perimeter(geom)[0]))

            # 以平方千米为单位
            gdf['area_km2'] = gdf['area_m2'] / 1e6

            # temp.append(gdf['area_km2'].sum())

            # 存储记录
            records.append({"id": name, "geometry": merged_geom,"area":gdf['area_km2'].sum()})

    # 创建最终输出 GeoDataFrame
    out_gdf = gpd.GeoDataFrame(records, crs=gdf.crs)

    # 写入到同一个 layer
    out_gdf.to_file(out_gpkg, layer=layer_name, driver="GPKG")

    print("Done:", out_gpkg)


def get_river_mouth_point(dir_file,stream_file,out_point_file,out_point_file_tif):
    """
    识别dir中的河口：8邻域内都不为空值的nodata像元
    :param dir_file:
    :param out_point_file:
    :param out_point_file_tif:
    :return:
    """
    fdir = Raster.get_raster(dir_file)
    proj,geo,f_nodata = Raster.get_proj_geo_nodata(dir_file)

    stream = Raster.get_raster(stream_file)
    proj,geo,s_nodata = Raster.get_proj_geo_nodata(stream_file)

    row,col = fdir.shape

    # 3. 创建坐标系对象
    spatial_ref = osr.SpatialReference()
    spatial_ref.ImportFromWkt(proj)  # 使用栅格的投影信息
    shapefile_driver = ogr.GetDriverByName("ESRI Shapefile")
    out_shapefile = shapefile_driver.CreateDataSource(out_point_file)

    # 3. 创建一个点的图层
    layer = out_shapefile.CreateLayer("points", spatial_ref, geom_type=ogr.wkbPoint)
    for i in range(row):
        for j in range(col):
            if stream[i,j] == s_nodata:
                continue
            now_dir = fdir[i,j]
            if now_dir not in dmove_dic:
                continue

            next_Cell = (i+dmove_dic[now_dir][0],j+dmove_dic[now_dir][1])
            if not check_boundary(next_Cell[0],next_Cell[1],row,col):
                continue
            if fdir[next_Cell[0],next_Cell[1]] == f_nodata:
                x_coord = geo[0] + j * geo[1] + j * geo[2]
                y_coord = geo[3] + i * geo[4] + i * geo[5]

                # 创建点几何对象
                point = ogr.Geometry(ogr.wkbPoint)
                point.AddPoint(x_coord, y_coord)

                # 创建属性字段（可选）
                feature = ogr.Feature(layer.GetLayerDefn())
                feature.SetGeometry(point)
                layer.CreateFeature(feature)


def get_endorheic_from_final_point(fdir_file,final_point_file,out_file,outcsv):
    """
    根据内流终点和流向提取endorheic，并记录面积、id
    :param fdir_file:
    :param final_point_file:
    :return:
    """

    fdir = Raster.get_raster(fdir_file)
    proj,geo,nodata = Raster.get_proj_geo_nodata(fdir_file)
    point = Raster.get_raster(final_point_file)
    proj,geo,p_nodata = Raster.get_proj_geo_nodata(final_point_file)

    row,col = point.shape
    area = []
    endorheic = np.zeros((row,col))
    vis = np.zeros((row,col),dtype = np.int8)
    endorheic_id = 1
    for i in range(row):
        for j in range(col):
            if point[i,j] == p_nodata:
                continue
            if vis[i,j] == 1:
                continue
            # 回溯上游
            vis[i,j] = 1
            popCells = [(i,j)]
            temp_area = 0
            while popCells:
                popCell = popCells.pop()
                temp_area += 1
                up_cells = sink.get_rever_D8(fdir,popCell[0],popCell[1])

                for up_cell in up_cells:
                    endorheic[up_cell[0],up_cell[1]] = endorheic_id
                    if vis[up_cell[0],up_cell[1]] == 0:
                        popCells.append(up_cell)

                        vis[up_cell[0], up_cell[1]] = 1
            area.append([endorheic_id,temp_area])
            endorheic_id += 1

    Raster.save_raster(out_file, endorheic, proj, geo, gdal.GDT_Float32, 0)
    with open(outcsv,'w',newline='') as f:
        writer = csv.writer(f)
        writer.writerows(area)
        f.close()



def get_SA():
    get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/Final_FDIR.tif",
                              "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/endorheic/SA_endorheic_point.shp",
                              "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/endorheic/SA_endorheic_point.tif")
    get_endorheic_from_final_point("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/Final_FDIR.tif",
                                   "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/endorheic/SA_endorheic_point.tif",
                                   "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/endorheic/SA_endorheic.tif",
                                   "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/endorheic/SA_endorheic_point.csv")



if __name__ == '__main__':
    get_river_mouth_point("/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Peral_river/ZJ_fdir.tif",
                          "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Peral_river/ZJ_stream.tif",
                          "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Peral_river/ZJ_mouth.shp",
                          "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Peral_river/ZJ_stream.tif")






    pass






