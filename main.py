# -*- coding: utf-8 -*-
"""
@Time ： 2024/9/20 19:39
@Auth ：
@File ：main.py.py
@IDE ：PyCharm
"""
import os
import shutil
import numpy as np

import postprocess
import special_process
import Merge
import calculate_sink_area_volume
import db
import process_sink
from global_sink import *
from genral_functions import *
from Run_all import *
def A(venuV,outVenu):
    """
    批量处理文件夹中的DEM:优先搜索边界上的0值加入队列，使用图搜索将8邻域内的0值栅格赋值为nodata
    :param venu:
    :param outVenu:
    :return:
    """
    # if not os.path.exists(outVenu):
    #     os.mkdir(outVenu)

    venus = os.listdir(venuV)
    for venu1 in venus:
        venu = os.path.join(venuV,venu1)
        files = os.listdir(venu)
        # po = Pool(25)
        for file in files:
            if file.split('.')[-1] != 'tif':
                continue

            demFile = os.path.join(venu,file)
            outFile = os.path.join(outVenu,file)
            # process_0(demFile,outFile)
            # po.apply_async(process_0,(demFile,outFile,))
        # po.close()
        # po.join()

# ---------------------------------------------------1.< 废弃代码 --------------------------------------------------------------- #
def main(db_name,dirVenu,sinkVenu,mergeDemPath,demVenu):


    dems = os.listdir(demVenu)
    visDEM = []
    for dem in dems:
        if dem.split('.')[-1] != 'tif':
            continue
        demfile = os.path.join(demVenu,dem)
        if demfile in visDEM:
            continue

        while True:
            outdirfile,outsinkfile = cal_dir_sink(demfile,dirVenu,sinkVenu)
            print(outsinkfile)
            get_sink(outsinkfile,db_name)
            inputRasters,extend = cal_neighbor(db_name,demVenu)
            if len(inputRasters) == 0:
                #
                break
            visDEM += inputRasters
            flag,demfile = Merge_cal_sir_sink(extend,inputRasters,mergeDemPath)


        # outdir,outsink = cal_dir_sink(outdemfile)
        # get_sink(outsink,db_name)
def judge_sink_type_main(demvenu,dirvenu,sinkvenu):

    demfiles = os.listdir(demvenu)

    for demfile in demfiles:
        if demfile.split('.')[-1] != 'tif':
            continue
        demPath = os.path.join(demvenu,demfile)
        sinkfile = 'sink_'+demfile
        sinkPath = os.path.join(sinkvenu,sinkfile)
        dirfile = 'dir_'+demfile
        dirPath = os.path.join(dirvenu,dirfile)
        calculate_sink_area_volume.calculate_sink_volume(demPath,dirPath,sinkPath)
def process_main(db_name,dirvenu,sinkvenu,outDir):
    """
    :param db_name:sinkdb.db
    :param dirvenu:
    :param sinkvenu:
    :param outDir:
    :return:
    """
    # 修正主体运行函数：根据每个数据库筛选dem文件，对每个文件进行修正，最后进行合并
    infos = db.query_sink(db_name)
    deleteInfos = db.query_sinkdelete(db_name)

    sinkFilesDic = {}
    sinkDeleteDic = {}
    for info in infos:
        sinkFilesDic.setdefault(info[-1],{}).setdefault(info[7],[]).append(info)
    for info in deleteInfos:
        sinkDeleteDic.setdefault(info[-1],[]).append(info[0])
    for demPath in sinkFilesDic:
        print(demPath,'...')
        sinkInfos = sinkFilesDic[demPath]
        demfile = os.path.basename(demPath)
        sinkPath = os.path.join(sinkvenu,'sink_'+demfile)
        dirPath = os.path.join(dirvenu,'dir_'+demfile)
        if demfile in sinkFilesDic:
            sinkDelete = sinkDeleteDic[demfile]
        else:
            sinkDelete = []
        process_sink.Allmain(demPath,dirPath,sinkPath,sinkInfos,sinkDelete,outDir)
def sbatch_main(extent,demDirVenu,OutDirVenu):
    """
    并行计算
    :param extent:[L,R,D,U]，有正负
    :param demDirVenu:
    :param OutDirVenu:
    :return:
    """

    left = extent[0]
    right = extent[1]
    down = extent[2]
    up = extent[3]

    Lon = list(range(left,right,10))
    Lat = list(range(down,up,10))
    print(Lon,Lat)
    demDirName = []  # 构建输入路径
    for i in range(len(Lon)):
        sLon1 = "{:03d}".format(abs(Lon[i]))
        sLon2 = "{:03d}".format(abs(Lon[i]+10))
        if Lon[i] >= 0:
            dLon1 = 'E'
        else:
            dLon1 = 'W'
        if Lon[i]+10 >= 0:
            dLon2 = 'E'
        else:
            dLon2 = 'W'
        for j in range(len(Lat)):
            sLat1 = "{:02d}".format(abs(Lat[j]))
            sLat2 = "{:02d}".format(abs(Lat[j]+10))
            if Lat[j] >= 0:
                dLat1 = 'N'
            else:
                dLat1 = 'S'
            if Lat[j]+10 >= 0:
                dLat2 = 'N'
            else:
                dLat2 = 'S'
            demDirPath = dLat1+sLat1+dLon1+sLon1+'-'+dLat2+sLat2+dLon2+sLon2+'_FABDEM_V1-0'
            temp = os.path.join(demDirVenu,demDirPath)
            temp_input = []
            if os.path.exists(temp):
                # demDirName.append(temp)
                outPath = os.path.join(OutDirVenu,demDirPath)
                if not os.path.exists(outPath):
                    os.mkdir(outPath)
                    os.chmod(outPath,0o777)
                dirPath = os.path.join(outPath,'originDir')
                sinkPath = os.path.join(outPath,'sink')
                modifiedDir = os.path.join(outPath,'modifiedDir')
                if not os.path.exists(dirPath):
                    os.mkdir(dirPath)
                    os.chmod(dirPath, 0o777)

                if not os.path.exists(sinkPath):
                    os.mkdir(sinkPath)
                    os.chmod(sinkPath, 0o777)
                if not os.path.exists(modifiedDir):
                    os.mkdir(modifiedDir)
                    os.chmod(modifiedDir, 0o777)
                temp_input += [temp,dirPath,sinkPath,os.path.join(outPath,dLat1+sLat1+dLon1+sLon1+'-'+dLat2+sLat2+dLon2+sLon2+'.db'),modifiedDir]
                demDirName.append(temp_input)

    # print(*demDirName)
    # 批量运行main计算洼地
    # po = Pool(25)
    # for info in demDirName:
    #     po.apply_async(main, (info[3],info[1],info[2],info[0],info[0],))
    # po.close()
    # po.join()
    #
    # # 批量运行cal计算洼地类型
    # po = Pool(25)
    # for info in demDirName:
    #     po.apply_async(judge_sink_type_main, (info[0], info[1], info[2],))
    # po.close()
    # po.join()

    # # 批量运行process_main计算洼地类型
    # po = Pool(25)
    # for info in demDirName:
    #     po.apply_async(process_main, (info[0], info[1], info[2],))
    # po.close()
    # po.join()


def mosica_main(ModifyVenu,output):

    files = os.listdir(ModifyVenu)
    mosicaFiles = [os.path.join(ModifyVenu,file) for file in files]
    mosaic1(mosicaFiles,output,0)
def caldir(demPath,outDir,outSink):


        wbt.d8_pointer(demPath,outDir,esri_pntr=True)
        wbt.sink(demPath,outSink)

def merge_America(extent,demDirVenu,OutMergeDEM):
    """
        并行计算
        :param extent:[L,R,D,U]，有正负
        :param demDirVenu:
        :param OutDirVenu:
        :return:
        """

    left = extent[0]
    right = extent[1]
    down = extent[2]
    up = extent[3]

    Lon = list(range(left, right, 10))
    Lat = list(range(down, up, 10))
    print(Lon, Lat)
    demName = []  # 构建输入路径
    for i in range(len(Lon)):
        sLon1 = "{:03d}".format(abs(Lon[i]))
        sLon2 = "{:03d}".format(abs(Lon[i] + 10))
        if Lon[i] >= 0:
            dLon1 = 'E'
        else:
            dLon1 = 'W'
        if Lon[i] + 10 >= 0:
            dLon2 = 'E'
        else:
            dLon2 = 'W'
        for j in range(len(Lat)):
            sLat1 = "{:02d}".format(abs(Lat[j]))
            sLat2 = "{:02d}".format(abs(Lat[j] + 10))
            if Lat[j] >= 0:
                dLat1 = 'N'
            else:
                dLat1 = 'S'
            if Lat[j] + 10 >= 0:
                dLat2 = 'N'
            else:
                dLat2 = 'S'
            demDirPath = dLat1 + sLat1 + dLon1 + sLon1 + '-' + dLat2 + sLat2 + dLon2 + sLon2 + '_FABDEM_V1-0'


            temp = os.path.join(demDirVenu, demDirPath)
            temp_input = []
            if os.path.exists(temp):
                dems = os.listdir(temp)
                for dem in dems:
                    if dem.split('.')[-1] != 'tif':
                        continue
                    demName.append(os.path.join(temp,dem))

    mosaic1(demName,OutMergeDEM)
    # print(demName)

def merge_SouthAmerica(extent,demDirVenu,process_demVenu,OutMergeDEM):
    """
    并行计算
    :param extent:[L,R,D,U]，有正负
    :param demDirVenu:
    :param OutDirVenu:
    :return:
    """

    left = extent[0]
    right = extent[1]
    down = extent[2]
    up = extent[3]

    Lon = list(range(left, right, 10))
    Lat = list(range(down, up, 10))
    print(Lon, Lat)
    demName = []  # 构建输入路径
    for i in range(len(Lon)):
        sLon1 = "{:03d}".format(abs(Lon[i]))
        sLon2 = "{:03d}".format(abs(Lon[i] + 10))
        if Lon[i] >= 0:
            dLon1 = 'E'
        else:
            dLon1 = 'W'
        if Lon[i] + 10 >= 0:
            dLon2 = 'E'
        else:
            dLon2 = 'W'
        for j in range(len(Lat)):
            sLat1 = "{:02d}".format(abs(Lat[j]))
            sLat2 = "{:02d}".format(abs(Lat[j] + 10))
            if Lat[j] >= 0:
                dLat1 = 'N'
            else:
                dLat1 = 'S'
            if Lat[j] + 10 >= 0:
                dLat2 = 'N'
            else:
                dLat2 = 'S'
            demDirPath = dLat1 + sLat1 + dLon1 + sLon1 + '-' + dLat2 + sLat2 + dLon2 + sLon2 + '_FABDEM_V1-0'

            temp = os.path.join(demDirVenu, demDirPath)
            temp_input = []
            if os.path.exists(temp):
                dems = os.listdir(temp)
                for dem in dems:
                    if dem.split('.')[-1] != 'tif':
                        continue
                    demName.append(os.path.join(process_demVenu, dem))

    mosaic1(demName, OutMergeDEM)
    # print(demName)

# --------------------------------------------------- 废弃代码 > 1.--------------------------------------------------------------- #

def copy_DEM_flat(venuV, outVenu):
    """
    将该目录下所有DEM复制到另一个文件夹
    :param venu:
    :param outVenu:
    :return:
    """
    if not os.path.exists(outVenu):
        os.mkdir(outVenu)

    venus = os.listdir(venuV)
    for venu1 in venus:
        venu = os.path.join(venuV, venu1)
        files = os.listdir(venu)

        for file in files:
            if file.split('.')[-1] != 'tif':
                continue
            demFile = os.path.join(venu, file)
            shutil.copy(demFile, outVenu)  # 将文件复制到目标文件夹，文件名保持不变
def mergeDEM(Venu,DEMPath):
    result = []
    files = os.listdir(Venu)
    for file in files:
        if file.split('.')[-1] != 'tif':
            continue
        result.append(os.path.join(Venu, file))

    mosaic1(result, DEMPath)

# 水体burning相关代码
def merge_GSWO(Venu,outGSWOPath):
    result = []
    files = os.listdir(Venu)
    for file in files:
        if file.split('.')[-1] != 'tif':
            continue
        result.append(os.path.join(Venu, file))

    mosaic1(result, outGSWOPath,255)
def merge_OSM(Venu,outOSMPath):

    result = []
    files = os.listdir(Venu)
    for file in files:
        if file.split('.')[-1] != 'tif':
            continue
        result.append(os.path.join(Venu,file))

    mosaic1(result,outOSMPath,-9)
def discretize_OSM(OSMfile,outOSMfile):
    """
    将OSM原本的1-5离散化为连续型:
    2\3 = 30%
    5  = 10%
    :param OSMfile:
    :param outOSMfile:
    :return:
    """
    OSM = Raster.get_raster(OSMfile)
    proj,geo,nodata = Raster.get_proj_geo_nodata(OSMfile)

    row,col = OSM.shape
    OSM[OSM == 4] = 10
    OSM[OSM == 1] = 0  # 将海洋的高程降低多点，保证沿海水流入海里
    OSM[OSM == 2] = 30
    OSM[OSM == 3] = 30
    OSM[OSM == 5] = 10
    OSM[OSM == nodata] = 0

    Raster.save_raster(outOSMfile,OSM,proj,geo,gdal.GDT_Byte,0)

def resample(input_raster,output_raster):
    """
    重采样
    :param input_raster:
    :param output_raster:
    :return:
    """
    from osgeo import gdal

    # 使用 GDAL Warp 调整分辨率
    gdal.Warp(
        output_raster,
        input_raster,
        xRes=0.00027777778,  # 设置目标分辨率（像素大小）
        yRes=0.00027777778,
        resampleAlg="bilinear",  # 选择重采样方法（这里是双线性插值）
        creationOptions=['COMPRESS=LZW', "TILED=True","BIGTIFF=IF_SAFER"]
    )
    print("Resampling complete. Output saved to:", output_raster)
def rerange_GSWO(GSWOfile,outGSWOfile):
    """
    将GSWO的区间从0-100变换为0-70，参考Yamazaki
    :param GSWOfile:
    :param outGSWOfile:
    :return:
    """
    GSWO = Raster.get_raster(GSWOfile)
    proj,geo,nodata = Raster.get_proj_geo_nodata(GSWOfile)

    arr = np.array(GSWO)
    arr[arr > 100] = 0
    result = (arr * 0.7).astype(np.uint8)

    Raster.save_raster(outGSWOfile,result,proj,geo,gdal.GDT_Byte,255)
    # print(result.max(),result.min())
def Clip_OSM_GSWO(DEMVenu):

    def clip_raster_raster(input_raster,output_raster,base_raster):
        """
        用栅格裁剪栅格并对齐
        :param input_raster:
        :param output_raster:
        :param base_raster:
        :return:
        """

        # 打开基准栅格
        base_ds = gdal.Open(base_raster)
        geo_transform = base_ds.GetGeoTransform()
        proj = base_ds.GetProjection()
        x_res = geo_transform[1]
        y_res = -geo_transform[5]

        # 使用 gdal.Warp 对齐
        gdal.Warp(
            output_raster,
            input_raster,
            format="GTiff",
            xRes=x_res,
            yRes=y_res,
            targetAlignedPixels=True,  # 保持像素对齐
            outputBounds=(geo_transform[0], geo_transform[3] + (base_ds.RasterYSize-1) * geo_transform[5],
                          geo_transform[0] + (base_ds.RasterXSize-1) * geo_transform[1], geo_transform[3]),
            dstSRS=proj,
            creationOptions=['COMPRESS=LZW', "TILED=True", "BIGTIFF=IF_SAFER"]
        )

        print(f"栅格裁剪和对齐完成，输出文件：{output_raster}")

    GSWOFILE = "/datanode05/zhangbin/FAB_hydrography/initial_data/GSWO/GSWO_global.tif"
    OSMFILE = "/datanode05/zhangbin/FAB_hydrography/initial_data/OSM/OSM_global.tif"
    DEMdirs = os.listdir(DEMVenu)
    for DEMdir in DEMdirs:
        DEMdirPath = os.path.join(DEMVenu,DEMdir)

        basefile = os.path.join(DEMdirPath,DEMdir+"_dem.tif")
        outGSWOfile = os.path.join(DEMdirPath, DEMdir + "_GSWO.tif")
        outOSMfile = os.path.join(DEMdirPath, DEMdir + "_OSM.tif")
        print(basefile)
        clip_raster_raster(GSWOFILE,outGSWOfile,basefile)
        clip_raster_raster(OSMFILE, outOSMfile, basefile)

        gswo = Raster.get_raster(outGSWOfile)
        s_row,g_col = gswo.shape
        print("GSWO = {:d}x{:d}".format(s_row,g_col))
        osm = Raster.get_raster(outOSMfile)
        o_row,o_col = osm.shape
        print("OSM = {:d}x{:d}".format(o_row, o_col))
        dem = Raster.get_raster(basefile)
        d_row,d_col = dem.shape
        print("dem = {:d}x{:d}".format(d_row, d_col))

    pass
def Occurence(GSWOfile,OSMfile,outOccurencefile):


    GSWO = Raster.get_raster(GSWOfile)
    proj,geo,G_nodata = Raster.get_proj_geo_nodata(GSWOfile)

    OSM = Raster.get_raster(OSMfile)
    _,_,O_nodata = Raster.get_proj_geo_nodata(OSMfile)

    GSWO_row,GSWO_col = GSWO.shape
    OSM_row,OSM_col = OSM.shape

    if not (GSWO_row&OSM_row and GSWO_col&OSM_col) :
        print(GSWO_row,GSWO_col,OSM_row,OSM_col,"栅格未对齐:")
        print(GSWOfile,OSMfile)
        return

    result = np.zeros((GSWO_row,GSWO_col),dtype = np.float32)

    result[:,:] = 0

    # 仅对有效值区域进行加法运算（或其他操作）
    # result1 = np.multiply(GSWO[mask], 0.7, casting="unsafe")
    # result = np.multiply(GSWO, 0.7, casting="unsafe")+OSM
    result = GSWO + OSM
    result[GSWO == G_nodata] = 255
    result[OSM == O_nodata] = 255
    Raster.save_raster(outOccurencefile,result,proj,geo,gdal.GDT_Byte,255)
def prepare_burning_data(OSMfile,outOSMfile,GSWOfile,outGSWOfile,outOccurencefile):
    discretize_OSM(OSMfile,outOSMfile)
    rerange_GSWO(GSWOfile,outGSWOfile)
    Occurence(outGSWOfile,outOSMfile,outOccurencefile)
def sbatch_occ(DEMVenu):
    DEMdirs = os.listdir(DEMVenu)
    for DEMdir in DEMdirs:
        DEMdirPath = os.path.join(DEMVenu, DEMdir)

        outOccfile = os.path.join(DEMdirPath, DEMdir + "_Occ.tif")
        GSWOfile = os.path.join(DEMdirPath, DEMdir + "_GSWO.tif")
        OSMfile = os.path.join(DEMdirPath, DEMdir + "_OSM.tif")
        outGSWOfile = os.path.join(DEMdirPath, DEMdir + "_GSWO_.tif")
        outOSMfile = os.path.join(DEMdirPath, DEMdir + "_OSM_.tif")

        # Occurence(GSWOfile,OSMfile,outOccfile)
        try:
            prepare_burning_data(OSMfile,outOSMfile,GSWOfile,outGSWOfile,outOccfile)
        except:
            print(DEMdir,'is fail')
    pass

def PLD_lake_Occ(Occ_file,PLD_file,out_occ_file):
    """
    将湖泊数据耦合到Occ图层中，保证干旱区的湖泊连续性
    :param Occ_file:
    :param PLD_file:
    :param out_occ_file:
    :return:
    """
    Occ = Raster.get_raster(Occ_file)
    proj,geo,o_nodata = Raster.get_proj_geo_nodata(Occ_file)
    row1, col1 = Occ.shape

    PLD = Raster.get_raster(PLD_file)
    proj,geo,p_nodata = Raster.get_proj_geo_nodata(PLD_file)
    row2, col2 = PLD.shape

    if not (row1&row2 and col2&col1) :
        print(row1,col1,row2,col2,"栅格未对齐:")
        print(Occ_file,PLD_file)
        return

    result = np.zeros((row1, col1), dtype=np.float32)

    result[:, :] = 0

    # 仅对有效值区域进行加法运算（或其他操作）
    # result1 = np.multiply(GSWO[mask], 0.7, casting="unsafe")
    # result = np.multiply(GSWO, 0.7, casting="unsafe")+OSM
    result = Occ.copy()
    mask = PLD != p_nodata
    result[mask] = 100
    Raster.save_raster(out_occ_file, result, proj, geo, gdal.GDT_Byte, 255)


def get_dam_info(damShp, geo):
    def lon2index(lon, lat, geo):

        col = int((lon - geo[0]) / geo[1])
        row = int((lat - geo[3]) / geo[5])

        return [row, col]

    # 设置文件路径
    shapefile_path = damShp
    # 打开 Shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(shapefile_path, 0)  # 0 表示只读模式

    daminfos = []  # 存放[ lon, lat]
    if data_source is None:
        print("无法打开文件")
    else:
        print("成功打开文件")
        # 获取图层
        layer = data_source.GetLayer()
        # 打印图层信息
        print(f"图层中水库数量: {layer.GetFeatureCount()}")

        # 遍历要素
        for feature in layer:
            temp_info = []
            lon_DAM = feature.GetField("LONG_DAM")
            lat_DAM = feature.GetField("LAT_DAM")
            lon_Riv = feature.GetField("LONG_RIV")
            lat_Riv = feature.GetField("LAT_RIV")
            # print(geo,lon2index(lon_Riv,lat_Riv,geo))
            if (lon_DAM != 0) and (lat_DAM != 0):
                # print(lat_DAM, lon_DAM)
                daminfos.append(lon2index(lon_DAM, lat_DAM, geo))
            else:
                daminfos.append(lon2index(lon_Riv, lat_Riv, geo))

    # 关闭数据源
    data_source = None
    return daminfos
def burning(ocean_file,demfile, Occfile, outBurningDEMfile, damShp=None):
    from scipy.ndimage import label

    # DEM = Raster.get_raster(demfile)
    # proj, geo, d_nodata = Raster.get_proj_geo_nodata(demfile)
    # x,y = DEM.shape
    # DEM1 = DEM[:, :y-1].copy()
    # Raster.save_raster(demfile, DEM1, proj, geo, gdal.GDT_Float32, d_nodata)

    DEM = Raster.get_raster(demfile)
    proj, geo, d_nodata = Raster.get_proj_geo_nodata(demfile)
    dem_row, dem_col = DEM.shape


    # Occ = Raster.get_raster(Occfile)
    # _, _, O_nodata = Raster.get_proj_geo_nodata(Occfile)
    # Y = Occ[:,:dem_col].copy()
    # Raster.save_raster(Occfile,Y,proj,geo,gdal.GDT_Byte,O_nodata)
    Occ = Raster.get_raster(Occfile)
    _, _, O_nodata = Raster.get_proj_geo_nodata(Occfile)

    Occ_row, Occ_col = Occ.shape
    # Ocean = Raster.get_raster(ocean_file)
    # proj, geo, ocean_nodata = Raster.get_proj_geo_nodata(ocean_file)
    # X = Ocean[:, :dem_col].copy()
    # Raster.save_raster(ocean_file, X, proj, geo, gdal.GDT_Byte, ocean_nodata)

    Ocean = Raster.get_raster(ocean_file)
    Ocean_mask = (Ocean == 1)
    # print(Occ.max(), Occ.min())


    print(dem_row, dem_col, Occ_row, Occ_col)

    if (dem_row & Occ_row and dem_col & Occ_col) == 0:
        print(dem_row, dem_col, Occ_row, Occ_col, "栅格未对齐")
        return

    mask = (DEM != d_nodata) & (Occ != O_nodata) & (Occ != 0)
    DEM[mask] -= np.multiply(Occ[mask],0.06)-3  # Occ burning
    # print(DEM.max(),DEM.min())

    DEM[Ocean_mask] = d_nodata
    # 处理水坝
    # 使用连通性分析
    structure = np.ones((3, 3))  # 8邻域连通性
    if damShp != None:
        damInfos = get_dam_info(damShp, geo)  # dam信息
        for info in damInfos:

            dammask = Occ[info[0] - 20:info[0] + 21, info[1] - 20:info[1] + 21].copy()
            demmask = DEM[info[0] - 20: info[0] + 21, info[1] - 20: info[1] + 21].copy()

            # 寻找
            dam_water_mask = dammask > 0


            labeled_array, num_features = label(dam_water_mask, structure=structure)
            pointValue = labeled_array[20, 21]  # 有点dam不准确，可能没落在水体上，此方法处理
            if pointValue == 0:
                tempDem = np.zeros_like(demmask)
                tempDem[:, :] = 9999
                tempDem[dam_water_mask] = demmask[dam_water_mask]
                minValue = np.min(tempDem)
                demmask[dam_water_mask] = minValue
                DEM[info[0] - 20: info[0] + 21, info[1] - 20: info[1] + 21] = demmask
            else:

                decrease_dem_mask = (labeled_array == pointValue)
                tempDem = np.zeros_like(demmask)
                tempDem[:, :] = 9999
                tempDem[decrease_dem_mask] = demmask[decrease_dem_mask]
                minValue = np.min(tempDem)
                demmask[decrease_dem_mask] = minValue

                DEM[info[0] - 20: info[0] + 21, info[1] - 20: info[1] + 21] = demmask
            print(demmask.min(), demmask.max(), minValue)


    Raster.save_raster(outBurningDEMfile,DEM,proj,geo,gdal.GDT_Float32,d_nodata)
def sbatch_burning(DEMVenu,GDWVenu):
    DEMdirs = os.listdir(DEMVenu)
    for DEMdir in DEMdirs:
        DEMdirPath = os.path.join(DEMVenu, DEMdir)
        DEMfile = os.path.join(DEMdirPath, DEMdir + "_dem.tif")
        Oceanfile = os.path.join(DEMdirPath, DEMdir + "_OSM.tif")
        Occfile = os.path.join(DEMdirPath, DEMdir + "_Occ.tif")
        damfile = os.path.join(GDWVenu,DEMdir + ".shp")
        outburnDEMfile = os.path.join(DEMdirPath, DEMdir + "_burnDEM.tif")
        # burning(Oceanfile, DEMfile, Occfile, outburnDEMfile, damfile)
        if os.path.exists(outburnDEMfile):
            continue
        try:
            burning(Oceanfile,DEMfile,Occfile,outburnDEMfile,damfile)
        except:
            print(DEMdir,'burning is fail')
    pass

def clip_next_level(clip_dir,now_dir):

   def Po_clip_DEM(DEM_file,mask_file,out_dem_file,out_dem_png):
        # 设置压缩选项
        warp_options = gdal.WarpOptions(
            cutlineDSName=mask_file,
            cropToCutline=True,
            dstNodata=-9999,
            creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"]  # 使用LZW压缩并启用瓦片存储
        )
        # 打开矢量文件并读取图层
        gdal.Warp(
            out_dem_file,
            DEM_file,
            options=warp_options
        )

        read_and_save_tif(out_dem_file, out_dem_png)

   def Po_clip_Occ(DEM_file, mask_file, out_dem_file, out_dem_png):
       # 设置压缩选项
       warp_options = gdal.WarpOptions(
           cutlineDSName=mask_file,
           cropToCutline=True,
           dstNodata=255,
           creationOptions=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"]  # 使用LZW压缩并启用瓦片存储
       )
       # 打开矢量文件并读取图层
       gdal.Warp(
           out_dem_file,
           DEM_file,
           options=warp_options
       )

       read_and_save_tif(out_dem_file, out_dem_png)


   files = os.listdir(clip_dir)
   baseName = os.path.basename(now_dir)
   next_venu = os.path.dirname(now_dir)
   # if not os.path.exists(next_venu):
   #     os.mkdir(next_venu)
   DEM_file = os.path.join(now_dir, baseName + "_burnDEM.tif")
   Occ_file = os.path.join(now_dir, baseName + '_Occ.tif')
   for file in files:
       if file.split('.')[-1] != 'shp':
           continue
       mask_file = os.path.join(clip_dir,file)
       shp_base_name = file.split('.')[0]
       save_venu = os.path.join(next_venu,shp_base_name)
       if not os.path.exists(save_venu):
           os.mkdir(save_venu)
       outDEM = os.path.join(save_venu,shp_base_name + "_burnDEM.tif")
       outDEM_png = os.path.join(save_venu, shp_base_name + "_burnDEM.png")
       outOcc = os.path.join(save_venu, shp_base_name + "_Occ.tif")
       outOcc_png = os.path.join(save_venu, shp_base_name + "_Occ.png")

       Po_clip_DEM(DEM_file,mask_file,outDEM,outDEM_png)
       Po_clip_Occ(Occ_file,mask_file,outOcc,outOcc_png)








def check_rowcol(demfile,GSWOfile,OSMfile):
    DEM = Raster.get_raster(demfile)
    GSWO = Raster.get_raster(GSWOfile)
    OSM = Raster.get_raster(OSMfile)

    d_row,d_col = DEM.shape
    G_row,G_col = GSWO.shape
    O_row,O_col = OSM.shape

    print(d_row,d_col)
    print(G_row,G_col)
    print(O_row,O_col)

def read_and_save_tif(tif_file, output_file, scale_factor=1.0, region=None, cmap='viridis'):
    """
    读取大尺寸的 .tif 文件并保存为图形文件 (.png 或 .jpg)。

    参数:
        tif_file (str): 输入的 .tif 文件路径。
        output_file (str): 输出保存文件路径 (支持 .png 和 .jpg)。
        scale_factor (float): 缩放因子，<1 表示缩小；>1 表示放大。
        region (tuple): 截取区域 (x_start, y_start, width, height)。为 None 时读取整个图像。
        cmap (str): 用于渲染图像的颜色映射表 (默认 'viridis')。
    """
    import rasterio
    from rasterio.plot import show
    import matplotlib.pyplot as plt
    import numpy as np
    from rasterio.enums import Resampling
    from rasterio.windows import Window
    with rasterio.open(tif_file) as src:
        # 如果指定了区域，则读取相应窗口
        if region:
            x_start, y_start, width, height = region
            window = Window(x_start, y_start, width, height)
            data = src.read(window=window,
                            out_shape=(src.count, int(height * scale_factor), int(width * scale_factor)),
                            resampling=Resampling.bilinear)
            transform = src.window_transform(window)
        else:
            # 设定最大显示尺寸 (例如：最大宽度或高度不超过 1024)
            max_size = 1024
            scale_factor = min(max_size / src.width, max_size / src.height)
            # 整体缩放读取
            data = src.read(
                out_shape=(
                    src.count,
                    int(src.height * scale_factor),
                    int(src.width * scale_factor)
                ),
                resampling=Resampling.bilinear
            )
            transform = src.transform

        # 绘制并保存图片
        fig, ax = plt.subplots(figsize=(8, 8))
        if data.shape[0] == 1:
            # 单波段显示
            show(data[0], transform=transform, ax=ax, cmap=cmap)
        elif data.shape[0] >= 3:
            # 多波段显示 (RGB)
            img = np.dstack([data[0], data[1], data[2]])
            ax.imshow(img)
        ax.set_title('TIF 图像展示')
        plt.axis('off')

        # 保存图片
        plt.savefig(output_file, bbox_inches='tight', pad_inches=0.1, dpi=300)
        plt.close()

        print(f"✅ 图像已成功保存为：{output_file}")


if __name__ == '__main__':


    # -------------------------1 data preprocess-----------------------------------
    preprocess()
    clip_code()
    prepare_occurance()
    run_burning()
    # -------------------------1 data preprocess-----------------------------------

    # -------------------------2 calculate depression CDD-----------------------------------
    run_clip_next_level()
    run_fdir_drpression()
    # -------------------------2 calculate depression CDD-----------------------------------


    # -------------------------3 calculate imformation of CDD-----------------------------------
    run_calculate_sink_infos()
    # -------------------------3 calculate imformation of CDD-----------------------------------

    # -------------------------4 modify flow direction-----------------------------------
    run_process_endorheic()
    run_process_exorheic()
    run_clip_outDir()
    run_merge()
    # -------------------------4 modify flow direction-----------------------------------

    # -------------------------5 calculate flow accumulation,endorheic points and rivers-----------------------------------
    run_acc_calculation()
    run_extract_stream()
    get_stream_order()
    run_check_endorheic_point()
    # -------------------------5 calculate flow accumulation,endorheic points and rivers-----------------------------------

    # -------------------------6 Validation of streams and basins-----------------------------------
    GRDC_valid()
    run_extract_stream_valid_basin()
    run_endorheic_basin_postprocess()
    run_stream_valid()
    run_extract_region_stream()
    # -------------------------6 Validation of streams and basins-----------------------------------







    pass