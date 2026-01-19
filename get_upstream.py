# -*- coding: utf-8 -*-
"""
@Time ： 2024/11/20 18:44
@Auth ：
@File ：get_upstream.py
@IDE ：PyCharm
"""
import os
import pickle
import sqlite3
import networkx as nx
from osgeo import ogr
from multiprocessing import Pool
import csv
# import dbfread

def find_upstream_basins(outlet_basin, dictPath, nxModel=False):
    """
    查找上游流域id并返回[]
    :param outlet_basin: 带查询流域id
    :param dictPath: .dic
    :param nxModel:
    :return:
    """

    if nxModel is True:
        G = nx.read_gpickle(dictPath)
        subGraph = nx.bfs_tree(G, outlet_basin, reverse=True)
        result = subGraph.nodes()

    else:
        # 加载拓扑字典
        with open(dictPath, "rb") as fs:
            topoDict = pickle.load(fs)
        fs.close()
        # 初始化返回结果
        upDict = topoDict.upDict

        result = []
        query_queue = [outlet_basin]
        num = 1
        # 查询上游
        while num > 0:
            temp = query_queue.pop()
            result.append(temp)
            num -= 1
            ups = upDict[temp]
            up_num = len(ups)
            if up_num > 0:
                num += up_num
                query_queue.extend(ups)

    return result

def create_output_shp_topoCat(basinList, filePath, outPath):
    """
    输出流域范围
    :param basinList: [流域id]
    :param filePath: 每个大洲level12的shp
    :param outPath: 输出路径
    :return:
    """

    inDs = ogr.Open(filePath)
    inLayer = inDs.GetLayer()
    srs = inLayer.GetSpatialRef()
    featureDefn = inLayer.GetLayerDefn()

    shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    outDs = shpDriver.CreateDataSource(outPath)
    outLayer = outDs.CreateLayer("data", srs=srs, geom_type=ogr.wkbMultiPolygon)
    # Add field values from input Layer
    for i in range(0, featureDefn.GetFieldCount()):
        fieldDefn = featureDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)

    id_filter = "Pfaf_ID IN ({})".format(", ".join(map(str, basinList)))
    # Apply the filter to the layer
    inLayer.SetAttributeFilter(id_filter)

    for feature in inLayer:
        outLayer.CreateFeature(feature)
    # print(Total_area)
    outLayer.SyncToDisk()
    outDs.Destroy()
    # inDs.ReleaseResultSet(inLayer)
    inDs.Destroy()
    #
    # return Area,Total_area

def check(filePath):
    inDs = ogr.Open(filePath)
    inLayer = inDs.GetLayer()
    srs = inLayer.GetSpatialRef()
    featureDefn = inLayer.GetLayerDefn()

    field_count = featureDefn.GetFieldCount()
    print(f"  Number of fields: {field_count}")
    num = inLayer.GetFeatureCount()
    print(num)
    # 输出字段信息
    for j in range(field_count):
        field_defn = featureDefn.GetFieldDefn(j)
        field_name = field_defn.GetName()
        field_type = field_defn.GetTypeName()
        print(f"    Field {j + 1}: {field_name} ({field_type})")

    # shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    # outDs = shpDriver.CreateDataSource(outPath)
    # outLayer = outDs.CreateLayer("data", srs=srs, geom_type=ogr.wkbMultiPolygon)
    # Add field values from input Layer
    # for i in range(0, featureDefn.GetFieldCount()):
    #     fieldDefn = featureDefn.GetFieldDefn(i)
    #     print(fieldDefn.GetName())
        # outLayer.CreateField(fieldDefn)

def gpkg2shp(input_gpkg,output_shp):
    from osgeo import ogr

    # 设置输入和输出文件路径
    # input_gpkg = "your_file.gpkg"  # 替换为你的GeoPackage文件路径
    # output_shp = "output_shapefile.shp"  # 替换为你希望保存的Shapefile路径

    # 打开 GeoPackage 文件
    driver = ogr.GetDriverByName('GPKG')
    data_source = ogr.Open(input_gpkg)

    # 检查文件是否成功打开
    if data_source is None:
        print("Failed to open the input GeoPackage file.")
    else:
        # 获取图层数量
        layer_count = data_source.GetLayerCount()
        for i in range(layer_count):
            layer = data_source.GetLayerByIndex(i)
            print(i,layer.GetName())
            # 复制字段
            layer_defn = layer.GetLayerDefn()
            for i in range(layer_defn.GetFieldCount()):
                field_defn = layer_defn.GetFieldDefn(i)
                print(field_defn.GetName())
        print(f"Total layers: {layer_count}")

        # 获取第一个图层（或根据需要选择不同的图层）
        layer = data_source.GetLayerByIndex(2)

        # 获取 Shapefile 驱动
        shp_driver = ogr.GetDriverByName('ESRI Shapefile')

        # 创建 Shapefile 输出文件
        if shp_driver is None:
            print("Shapefile driver not available.")
        else:
            # 创建输出 Shapefile 数据源
            output_data_source = shp_driver.CreateDataSource(output_shp)

            # 创建图层（复制输入图层的结构和投影）
            output_layer = output_data_source.CreateLayer(layer.GetName(), geom_type=layer.GetGeomType())

            # 复制字段
            layer_defn = layer.GetLayerDefn()
            for i in range(layer_defn.GetFieldCount()):
                field_defn = layer_defn.GetFieldDefn(i)
                output_layer.CreateField(field_defn)

            # 复制特征（要素）
            for feature in layer:
                output_layer.CreateFeature(feature)

            print(f"GeoPackage to Shapefile conversion successful. Output saved to: {output_shp}")
