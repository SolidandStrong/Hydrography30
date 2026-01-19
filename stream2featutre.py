# -*- coding: utf-8 -*-
"""
@Time ： 2025/11/22 14:17
@Auth ：
@File ：stream2featutre.py
@IDE ：PyCharm
"""
import numpy as np
from shapely.geometry import LineString, Point
import geopandas as gpd
from osgeo import gdal,ogr,osr

# D8 流向编码 → 行列偏移
import Raster

d8_dict = {
    1:  (0, 1),   # E
    2:  (1, 1),   # SE
    4:  (1, 0),   # S
    8:  (1, -1),  # SW
    16: (0, -1),  # W
    32: (-1, -1), # NW
    64: (-1, 0),  # N
    128:(-1, 1),  # NE
}

def raster_to_vector_rivers(flowdir, river_id, transform):
    """
    flowdir: 2D numpy array, D8 directions
    river_id: 2D numpy array, river segment id
    transform: rasterio Affine，用于坐标转换
    """
    results = []

    unique_ids = np.unique(river_id)
    unique_ids = unique_ids[unique_ids > 0]  # 去除 0

    for rid in unique_ids:
        coords = np.argwhere(river_id == rid)

        # 将像元放入 dict 方便查找
        cells = {(r, c): True for r, c in coords}

        # 找到上游点（没有别的河流流向它）
        indegree = {tuple(x): 0 for x in coords}
        for r, c in coords:
            d = flowdir[r, c]
            if d not in d8_dict:
                continue
            dr, dc = d8_dict[d]
            nr, nc = r + dr, c + dc
            if (nr, nc) in indegree:
                indegree[(nr, nc)] += 1

        # 上游起点（indegree=0）
        sources = [key for key, val in indegree.items() if val == 0]
        if len(sources) == 0:
            continue  # 避免环

        start = sources[0]
        line = []

        # 追踪下游
        curr = start
        while True:
            r, c = curr
            # 转成坐标（中心点）
            x, y = transform * (c + 0.5, r + 0.5)
            line.append((x, y))
            if (r, c) not in cells:
                break  # 出河段
            d = flowdir[r, c]
            if d not in d8_dict:
                break

            dr, dc = d8_dict[d]
            nr, nc = r + dr, c + dc
            # if (nr, nc) not in cells:
            #     break  # 出河段
            curr = (nr, nc)

        # 生成 LineString
        if len(line) > 1:
            results.append({
                "river_id": rid,
                "geometry": LineString(line)
            })

    return gpd.GeoDataFrame(results, crs="EPSG:4326")

def stream_raster2feature(fdr_file,stream_file,out_file):
    import rasterio

    with rasterio.open(fdr_file) as ds:
        fdr = ds.read(1)
        transform = ds.transform  # 必为 Affine
        crs = ds.crs


    stream = Raster.get_raster(stream_file)

    gdf = raster_to_vector_rivers(fdr,stream,transform)

    # gdf = raster_to_vector_rivers(flowdir, river_id, transform)
    gdf.to_file(out_file)

# ---------------------
# 示例调用
# ---------------------
# 假设你已经读了 rasterio 打开的数据
# flowdir = ds.read(1)
# river_id = ds2.read(1)
# transform = ds.transform

# gdf = raster_to_vector_rivers(flowdir, river_id, transform)
# gdf.to_file("rivers.gpkg")

if __name__ == '__main__':

    ds = ogr.Open(r'F:\提取指定点上游的参数\Discription\Nandita\canada_river\test\3.shp')

    layer = ds.GetLayer(0)

    feature = layer.GetFeatureCount()

    for fea in layer:
        print(fea)
        print("The spatial geometry is",fea.geometry())

