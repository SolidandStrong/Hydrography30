# -*- coding: utf-8 -*-
"""
@Time ： 2025/11/19 11:19
@Auth ：
@File ：stream_valid.py
@IDE ：PyCharm
"""
import csv


from shapely.ops import nearest_points
from pyproj import Geod
import geopandas as gpd
from shapely.strtree import STRtree
from shapely.ops import nearest_points
from pyproj import Geod


def valid_stream(valid_point,valid_stream,out_point,out_csv):
    """
    计算SWROD点到valid河流的空间距离
    :param valid_point:
    :param valid_stream:
    :param out_point:
    :param out_csv:
    :return:
    """

    # # ------------------------------
    # # 1. 读取数据
    # # ------------------------------
    # points = gpd.read_file(valid_point)  # 验证点
    # rivers = gpd.read_file(valid_stream)  # 河流 shp (line)
    #
    # # 确保坐标是 EPSG:4326
    # points = points.to_crs(4326)
    # rivers = rivers.to_crs(4326)
    #
    # # 建立一个 Geod 对象用于测距
    # geod = Geod(ellps="WGS84")
    #
    # # ------------------------------
    # # 2. 找每个点最近的河流点，并计算距离
    # # ------------------------------
    # dx_list = []
    # dy_list = []
    # d_list = []
    #
    # results = [['dx_m','dy_m','ds_m']]
    # for pt in points.geometry:
    #     # 找最近河段
    #     nearest_geom = rivers.geometry.distance(pt).sort_values().index[0]
    #     river_line = rivers.geometry.loc[nearest_geom]
    #
    #     # 获取 pt 到河流最近点（投影点）
    #     nearest_pt = nearest_points(pt, river_line)[1]
    #
    #     # 提取经纬度
    #     lon1, lat1 = pt.x, pt.y
    #     lon2, lat2 = nearest_pt.x, nearest_pt.y
    #
    #     # ---- 计算 dx (东西方向距离) ----
    #     # 只改经度，经向距离 = constant latitude
    #     dx = geod.line_length([lon1, lon2], [lat1, lat1])
    #
    #     # ---- 计算 dy (南北方向距离) ----
    #     # 只改纬度，纬向距离 = constant longitude
    #     dy = geod.line_length([lon1, lon1], [lat1, lat2])
    #
    #     # ---- 计算最短球面距离 ----
    #     _, _, d = geod.inv(lon1, lat1, lon2, lat2)
    #
    #     dx_list.append(dx)
    #     dy_list.append(dy)
    #     d_list.append(d)
    #     results.append([dx,dy,d])
    #
    #
    # # 保存到新的列
    # points["dx_m"] = dx_list  # 经向距离（m）
    # points["dy_m"] = dy_list  # 纬向距离（m）
    # points["dist_m"] = d_list  # 最短距离（m）
    #
    # # ------------------------------
    # # 3. 输出结果
    # # ------------------------------
    # points.to_file(out_point)
    #
    # with open(out_csv,'w') as f:
    #     writer = csv.writer(f)
    #     writer.writerows(results)
    #     f.close()

    # 读取
    points = gpd.read_file(valid_point).to_crs(4326)

    rivers = gpd.read_file(valid_stream)  # 河流 shp (line)
    rivers = rivers.set_crs("EPSG:4326")
    # rivers = gpd.read_file(valid_stream).to_crs(4326)

    geod = Geod(ellps="WGS84")

    # 构建 STRtree（非常快）
    # tree = STRtree(rivers.geometry.values)
    # geoms = rivers.geometry.values  # list of river lines

    rivers_geom = rivers.geometry.values
    tree = STRtree(rivers_geom)

    # 存储
    dx_list = []
    dy_list = []
    dist_list = []
    results = [['lon','lat','dx_m','dy_m','ds_m']]
    for pt in points.geometry:
        # 找最近河段
        # idx = tree.nearest(pt)  # 返回索引
        # nearest_line = rivers.geometry.iloc[idx]  # 用索引取几何
        # # nearest_line = tree.nearest(pt)
        # nearest_pt = nearest_points(pt, nearest_line)[1]
        idx = tree.nearest(pt)  # 最近河段的索引
        nearest_line = rivers_geom[idx]  # 取真正的几何对象
        nearest_pt = nearest_points(pt, nearest_line)[1]

        lon1, lat1 = pt.x, pt.y
        lon2, lat2 = nearest_pt.x, nearest_pt.y

        # 经向偏差
        dx = geod.line_length([lon1, lon2], [lat1, lat1])
        # 纬向偏差
        dy = geod.line_length([lon1, lon1], [lat1, lat2])
        # 最短距离
        _, _, d = geod.inv(lon1, lat1, lon2, lat2)

        dx_list.append(dx)
        dy_list.append(dy)
        dist_list.append(d)
        results.append([lon1,lat1,dx, dy, d])

    points["dx_m"] = dx_list
    points["dy_m"] = dy_list
    points["dist_m"] = dist_list

    points.to_file(out_point)

    with open(out_csv,'w') as f:
        writer = csv.writer(f)
        writer.writerows(results)
        f.close()

