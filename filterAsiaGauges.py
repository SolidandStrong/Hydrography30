import os
import glob
from osgeo import ogr, osr
import numpy as np


"""
    提取位于亚洲区域的GRDC水文站流域
    取所有与亚洲区域外包矩形相交的流域

"""


def main(fpath,shp,outtxt):

    # Asia_ds = ogr.Open(r"D:\data\raw_data\HydroBASINS\hybas_as_lev01-12_v1c\hybas_as_lev01_v1c.shp")
    Asia_ds = ogr.Open(shp)
    Asia_layer = Asia_ds.GetLayer(0)
    feature = Asia_layer.GetFeature(0)
    Asia_geom = feature.GetGeometryRef()
    Asia_envelope = Asia_geom.GetEnvelope()
    Asia_ds.Destroy()

    # 筛选亚洲区域的流域
    os.chdir(fpath)
    pattern = "grdc_basins_smoothed_md_no_???????.shp"
    files = glob.glob(pattern)

    filter_stations = []
    # 筛选位于HydroBASINS亚洲区域的流域
    # 且流域面积大于100km2
    for fn in files:
        ds = ogr.Open(fn)
        layer = ds.GetLayer(0)
        inLayerDefn = layer.GetLayerDefn()

        if inLayerDefn.GetFieldCount() > 0:
            feature = layer.GetFeature(0)
            # 提取面积
            area = feature.GetField("AREA")
            # 提取编号
            grdc_no = feature.GetField("GRDC_NO")

            # 判断面积是否大于100km2
            if area > 100.0:
                lon = feature.GetField("LONG_NEW")
                lat = feature.GetField("LAT_NEW")
                if Asia_envelope[0] < lon < Asia_envelope[1] and Asia_envelope[2] < lat < Asia_envelope[3]:
                    filter_stations.append(grdc_no)

        ds.Destroy()

    print(filter_stations)
    print(len(filter_stations))

    with open(outtxt, "w") as fs:
        for grdc_no in filter_stations:
            fs.write("%d\n" % grdc_no)
    fs.close()

def sbatch_main():
    GRDC_folder = "/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/statbas_shp_zip/"
    continents = ['Asia','NorthAmerica','SouthAmerica','Africa','Europe','Siberia','Greenland','Arctic','Australia']
    shp_venu = "/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/hyb_shp/"
    filter_stations = "/datanode05/zhangbin/FAB_hydrography/initial_data/GRDC/filter_stations/"
    for continent in continents:
        try:
            main(GRDC_folder,os.path.join(shp_venu,continent+'.shp'),os.path.join(filter_stations,continent+'.txt'))
        except Exception as e:
            print(continent,e)
            pass


def filter_station():

    path = r"D:\map\MasterThesis\Validation\vector\filterStations.shp"
    ds = ogr.Open(path)
    layer = ds.GetLayer()

    filter_stations = []
    for feature in layer:
        grdc_no = feature.GetField("GRDC_NO")
        filter_stations.append(grdc_no)

    with open(r"D:\论文\毕业论文\作图分析代码\数据集评估\filter_stations.txt", "w") as fs:
        for grdc_no in filter_stations:
            fs.write("%d\n" % grdc_no)
    fs.close()


def merge_stations(fpath, spath):
    """
    将所有的站点范围合并到一个图层
    polygon和point各一个
    :return:
    """

    # 读取站点列表
    stations = np.loadtxt(spath, dtype=np.int32)
    stations = stations.ravel()

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)


    # 创建两个数据集，polygon和point
    driver = ogr.GetDriverByName("ESRI Shapefile")
    polygonDs = driver.CreateDataSource("polygon.shp")
    polygonLayer = polygonDs.CreateLayer("data", srs=srs, geom_type=ogr.wkbPolygon)
    pointDs = driver.CreateDataSource("point.shp")
    pointLayer = pointDs.CreateLayer("data", srs=srs, geom_type=ogr.wkbPoint)

    # 获取属性字段定义
    station = stations[0]
    fn = os.path.join(fpath, "grdc_basins_smoothed_md_no_%d.shp" % station)
    refDs = ogr.Open(fn)
    inLayer = refDs.GetLayer(0)
    # Add input Layer Fields to the output Layer if it is the one we want
    inLayerDefn = inLayer.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        polygonLayer.CreateField(fieldDefn)
        pointLayer.CreateField(fieldDefn)
    refDs.Destroy()

    polygonLayerDefn = polygonLayer.GetLayerDefn()
    pointLayerDefn = pointLayer.GetLayerDefn()


    for station in stations:
        fn = os.path.join(fpath, "grdc_basins_smoothed_md_no_%d.shp" % station)
        inDs = ogr.Open(fn)
        inLayer = inDs.GetLayer(0)
        inLayerDefn = inLayer.GetLayerDefn()
        feature = inLayer.GetFeature(0)

        lon = feature.GetField("LONG_NEW")
        lat = feature.GetField("LAT_NEW")

        polygonFeature = ogr.Feature(polygonLayerDefn)
        pointFeature = ogr.Feature(pointLayerDefn)

        for i in range(0, inLayerDefn.GetFieldCount()):
            polygonFeature.SetField(i, feature.GetField(i))
            pointFeature.SetField(i, feature.GetField(i))

        geom = feature.GetGeometryRef()
        polygonFeature.SetGeometry(geom)
        polygonLayer.CreateFeature(polygonFeature)
        pgeom = ogr.Geometry(ogr.wkbPoint)
        pgeom.AddPoint(lon, lat)
        pointFeature.SetGeometry(pgeom)
        pointLayer.CreateFeature(pointFeature)

        inDs.Destroy()
        inLayer = None

    pointDs.Destroy()
    polygonDs.Destroy()



if __name__ == "__main__":
    folder = r'F:\全球数据集\图表\素材\画图data\GRDC\statbas_shp_zip'#r"D:\Downloads\statbas_shp_zip"
    # main(folder)
    # filter_station()
    # merge_stations(folder, "filter_stations.txt")
