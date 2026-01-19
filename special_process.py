# -*- coding: utf-8 -*-
"""

处理那些由于内存过大而导致提取sink失败的流域
@Time ： 2025/3/1 12:48
@Auth ：
@File ：special_process.py
@IDE ：PyCharm
"""
import csv
import heap_PF_D8
import math
import os
import shutil

import rasterio

import split_cal_acc
from genral_functions import *
import application_endorheic
import calculate_sink_area_volume
import Find
from multiprocessing import Pool
import sink
import Merge
import split
import valid
import Raster
import process_sink
import postprocess
import endorheic_process
from osgeo import gdal,ogr
import numpy as np
def Find_null_sink(venu):

    dirNames = os.listdir(venu)
    for dirName in dirNames:
        dirPath = os.path.join(venu,dirName)
        sink_file = os.path.join(dirPath,dirName+'_burnsink_.tif')
        if not os.path.exists(sink_file):
            print(dirName)

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
def preprocess(buffervenu,dem_file,occ_file,venu):
    """
    根据buffer文件夹来生成根目录，包含裁剪后的dem、occ
    :param buffervenu: 包含buffer的文件夹
    :param dem_file:
    :param occ_file:
    :param venu: 州文件夹.不加\
    :return:
    """

    baseName = os.path.basename(venu)
    buffer_Names = os.listdir(buffervenu)

    params = []
    for buffer_Name in buffer_Names:
        buffer_Name_s = buffer_Name.split('.')
        if buffer_Name_s[1] != 'shp':
            continue
        temp_venu = os.path.join(venu,baseName + '_' + buffer_Name_s[0])
        if not os.path.exists(temp_venu):
            os.mkdir(temp_venu)
            os.chmod(temp_venu,0o777)
        outdem = os.path.join(temp_venu,baseName + '_' + buffer_Name_s[0] + '_burnDEM_.tif')
        outocc = os.path.join(temp_venu, baseName + '_' + buffer_Name_s[0] + '_Occ_.tif')
        temp_paras = [os.path.join(buffervenu,buffer_Name),dem_file,outdem,occ_file,outocc]
        params.append(temp_paras)
        Find.clip(dem_file,os.path.join(buffervenu,buffer_Name),outdem)
        Find.clip(occ_file, os.path.join(buffervenu, buffer_Name), outocc)

def cal_dir_sink(buffervenu,venu):
    """
    根据buffer文件夹来计算根目录下的sink dir
    :param buffervenu:
    :return:
    """
    baseName = os.path.basename(venu)
    buffer_Names = os.listdir(buffervenu)
    for buffer_Name in buffer_Names:
        buffer_Name_s = buffer_Name.split('.')
        if buffer_Name_s[1] != 'shp':
            continue
        temp_venu = os.path.join(venu,baseName + '_' + buffer_Name_s[0])

        outdem = os.path.join(temp_venu,baseName + '_' + buffer_Name_s[0] + '_burnDEM_.tif')

        try:
            outsink = os.path.join(temp_venu, baseName + '_' + buffer_Name_s[0] + '_burnsink_.tif')
            if os.path.exists(outsink):
                continue
            sink.Cal_sink(outdem,outsink)
        except Exception as e:
            print(outdem,e)
        outdir = os.path.join(temp_venu, baseName + '_' + buffer_Name_s[0] + '_burndir_.tif')
        if os.path.exists(outdir):
            continue
        sink.Cal_dir(outdem, outdir)

def PFD8(burn_dem_file,outdirfile):
    """
    按照内流区使用PFD8算法计算流向
    :param burn_dem_file:
    :param outdirfile: 输出的MOdifiedDir.tif
    :return:
    """
    dem = Raster.get_raster(burn_dem_file)
    proj,geo,d_nodata = Raster.get_proj_geo_nodata(burn_dem_file)


    # 找到最小值
    minCell = [0,0]
    minValue = 1000
    row,col = dem.shape
    for i in range(row):
        for j in range(col):
            if dem[i,j] == d_nodata:
                continue
            if dem[i,j] <= minValue:
                minValue = dem[i,j]
                minCell = [i,j]
    mask_dir = np.zeros((row,col))
    mask_dir[minCell[0],minCell[1]] = 0

    A = endorheic_process.temp_queuePriorityFlood(dem,d_nodata,dem,minCell,mask_dir)

    Raster.save_raster(outdirfile,A,proj,geo,gdal.GDT_Byte,0)

def PFD8_exorheic(dem_file,sink_file,fdir_file,outfdir_path,out_lon,out_lat,new_fdr):

    """
    postprocess中对对流栅格的上游实施breach策略：
    1）搜索边界上的高程最低的栅格
    2）计算初始流向：非掩膜区域外的坡度最低栅格
    3)构建数据参数：masksink、maskdem、maskdir
    4）保存修正后的流向
    :param dem_file:
    :param sink_file:
    :param fdir_file:
    :return:
    """
    dem = Raster.get_raster(dem_file)
    proj,geo,d_nodata = Raster.get_proj_geo_nodata(dem_file)

    sink = Raster.get_raster(sink_file)
    proj,geo,s_nodata = Raster.get_proj_geo_nodata(sink_file)

    fdir = Raster.get_raster(fdir_file)
    proj,geo,f_nodata = Raster.get_proj_geo_nodata(fdir_file)

    row,col = sink.shape
    # minHCell = (0,0)
    # minH = math.inf
    # for i in range(row):
    #     for j in range(col):
    #         if sink[i,j] == s_nodata:
    #             continue
    #         # 判断是否为边界栅格
    #         H = dem[i,j]
    #         flag = False
    #         for k in range(8):
    #             nextCell = (i+dmove[k][0],j+dmove[k][1])
    #             if not check_boundary(nextCell[0],nextCell[1],row,col):
    #                 continue
    #             if sink[nextCell[0],nextCell[1]] == s_nodata:
    #                 # 则(i,j) 是边界栅格
    #                 flag = True
    #                 break
    #         if flag:
    #             # 比较当前栅格，取最低高程栅格
    #             if H <= minH:
    #                 minH = H
    #                 minHCell = (i,j)
    # # 计算初始流向
    # initial_Dir = f_nodata
    # min_hillslope = - math.inf
    # for k in range(8):
    #     nextCell = (minHCell[0]+dmove[k][0],minHCell[1]+dmove[k][1])
    #     if not check_boundary(nextCell[0], nextCell[1], row, col):
    #         continue
    #     if sink[nextCell[0],nextCell[1]] != s_nodata:
    #         continue
    #     if k in [1,3,5,7]:
    #         kk = math.sqrt(2)
    #     else:
    #         kk = 1
    #     hillslope = (minH - dem[nextCell[0],nextCell[1]]) / kk
    #     if hillslope > min_hillslope:
    #         min_hillslope = hillslope
    #         initial_Dir = int(2**k)
    # fdir[minHCell[0], minHCell[1]] = initial_Dir
    # ----------------------- 手动计算行列号 ----------------------
    lon = out_lon#49.4902437#-58.9836659#-58.1573169#-59.8978350#49.4904803#-135.6474676#-115.0928235##-59.9152881##
    lat = out_lat#53.4589140#-3.2430144#-32.796551#-3.1296820#53.4617709#68.7559798#61.3897567##-3.1375043#
    col1 = int((lon - geo[0]) / geo[1])
    row1 = int((lat - geo[3]) / geo[5])
    minHCell = (row1,col1)
    print(row1, col1)
    print(fdir[minHCell[0],minHCell[1]])
    fdir[minHCell[0], minHCell[1]] = new_fdr#2#1#2#4#1
    # ----------------------- 手动计算行列号 ----------------------


    # 构建数据参数：masksink、maskdem、maskdir
    sink[sink != 1] = 0  # 写成0 是因为好用：否则后面查最小值的时候会错误
    maskDEM = dem.copy()
    maskDEM[sink != 1] = d_nodata
    maskDir = fdir.copy()
    maskDir[sink != 1] = f_nodata

    # 计算淹没内的流向
    # A = process_sink.priorityFlood_special(sink, 0, maskDEM, minHCell, maskDir)  # Python
    A = heap_PF_D8.optimized_flow_repair(maskDEM,sink,0,minHCell)   # 堆栈优化后的版本
    # print(A)
    A[sink != 1] = f_nodata
    A[minHCell[0], minHCell[1]] = new_fdr  # 2#1#2#4#1

    # A = C_PriorityFlood(maskSink, 0, maskDem, outCell[0] - minX, outCell[1] - minY,row,col)  # C语言
    # print(A)

    Raster.save_raster(outfdir_path,A,proj,geo,gdal.GDT_Byte,f_nodata)
    # # os.chmod(outfdir_path,0o777)

def modify_conuter_flow(dem_file,sink_file,fdir_file,outfdir_file):
    """
        postprocess中对对流栅格实施breach策略：
        1）找到对流栅格所有的上游并进行唯一编码1-00
        2）让对流栅格的流向流出：搜索两个对流栅格的邻域，如果有一个对流栅格存在不同的编码则流出，另一个不变；如果两个都没有，则更新高程较低的那个对流栅格为内流终点。
        3）保存修正后的流向
        :param dem_file:
        :param sink_file:
        :param fdir_file:
        :return:
        """
    # 正常填洼

    # if os.path.exists(outfdir_file):
    #     return
    # print(outfdir_file,' is Filled!')
    proj, geo, f_nodata = Raster.get_proj_geo_nodata(fdir_file)

    import whitebox
    wbt = whitebox.WhiteboxTools()
    venu = os.path.dirname(dem_file)

    sink = Raster.get_raster(sink_file)
    proj, geo, s_nodata = Raster.get_proj_geo_nodata(sink_file)
    # new_dem_file = os.path.join(venu,'new_dem.tif')
    # dem = Raster.get_raster(dem_file)
    # proj,geo,d_nodata = Raster.get_proj_geo_nodata(dem_file)
    # dem[sink == s_nodata] = d_nodata
    # Raster.save_raster(new_dem_file,dem,proj,geo,gdal.GDT_Float32,d_nodata)

    temp_outFillfile = os.path.join(venu,'filled_dem.tif')
    temp_outtdirfile = os.path.join(venu,'temp_fdir1.tif')
    wbt.fill_depressions_planchon_and_darboux(dem_file, temp_outFillfile, fix_flats=True)
    wbt.d8_pointer(temp_outFillfile, temp_outtdirfile, esri_pntr=True)


    fdir = Raster.get_raster(temp_outtdirfile)

    fdir[sink == s_nodata] = f_nodata

    Raster.save_raster(outfdir_file,fdir,proj,geo,gdal.GDT_Byte,f_nodata)

def LZW(input_file,output_file):

    input1 = Raster.get_raster(input_file)
    proj,geo,nodata = Raster.get_proj_geo_nodata(input_file)

    Raster.save_raster(output_file,input1,proj,geo,gdal.GDT_Byte,nodata)

def sbatch_clip_ModifiedDir(DEM_file,venu,maskvenu):
    """
    原来是裁剪流向的代码，注释掉的部分

    现在是裁剪initial DEM的
    :param venu:
    :param maskvenu:
    :return:
    """
    # dirNames = os.listdir(venu)
    # for dirName in dirNames:
    #     if len(dirName.split('.'))>1:
    #         continue
    #
    #     dirPath = os.path.join(venu,dirName)
    #     mask_file = os.path.join(maskvenu,dirName+'.shp')
    #     DIRdile = os.path.join(dirPath,dirName+'_ModifiedDir_.tif')
    #     Merge.clip_modifiedDir(mask_file,DIRdile,dirPath)


    if not os.path.exists(venu):
        os.mkdir(venu)

    def Po_clip(DEM_file,mask_file,out_dem_file,out_dem_png):
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

    # DEM_file = "/data/zhangbin/TBP_Stream_Data/origin_raster/continent_origin_dem/Asia.tif" #"/data/zhangbin/TBP_Stream_Data/origin_raster/FABDEM_global.tif"
    # GSWO_file = "/data/zhangbin/TBP_Stream_Data/origin_raster/GSWO_global.tif"
    # OSM_file = "/data/zhangbin/TBP_Stream_Data/origin_raster/OSM_global.tif"

    paras = []
    for dirName in os.listdir(maskvenu):
        if dirName.split('.')[-1] != 'shp':
            continue

        outpath = os.path.join(venu,dirName.split('.')[0])
        if not os.path.exists(outpath):
            os.mkdir(outpath)
        mask_file = os.path.join(maskvenu,dirName)
        out_dem_file = os.path.join(outpath,dirName.split('.')[0]+'_dem.tif')
        out_dem_png = os.path.join(outpath, dirName.split('.')[0] + '_dem.png')
        paras.append([DEM_file,mask_file,out_dem_file,out_dem_png])
        # print([DEM_file,mask_file,out_dem_file,out_dem_png])

    # Po = Pool(10)
    for para in paras:
        Po_clip(para[0], para[1], para[2], para[3])
        # Po.apply_async(Po_clip,(para[0],para[1],para[2],para[3],))

    # Po.close()
    # Po.join()




def s_shp_si(shpFile,out_shp):
    # out_shp = shpFile[:-4] + '_ssi' + '.shp'
    # out_shp =shpFile
    # print(out_shp)
    # 打开数据
    ds = ogr.Open(shpFile, 0)
    if ds is None:
        print("打开文件 %s 失败！" % shpFile)
        return
    print("打开文件%s成功！" % shpFile)
    # 获取该数据源中的图层个数，一般shp数据图层只有一个，如果是mdb、dxf等图层就会有多个
    m_layer_count = ds.GetLayerCount()
    m_layer = ds.GetLayerByIndex(0)
    if m_layer is None:
        print("获取第%d个图层失败！\n", 0)
        return

    # 创建输出文件
    driver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(out_shp):
        driver.DeleteDataSource(out_shp)
    outds = driver.CreateDataSource(out_shp)
    outlayer = outds.CreateLayer(out_shp[:-4], m_layer.GetSpatialRef(),geom_type=ogr.wkbPolygon)
    # 获取输出层的要素定义
    outLayerDefn = outlayer.GetLayerDefn()
    # 对图层进行初始化，如果对图层进行了过滤操作，执行这句后，之前的过滤全部清空
    m_layer.ResetReading()
    # 获取投影
    prosrs = m_layer.GetSpatialRef()
    # 添加字段
    inLayerDefn = m_layer.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outlayer.CreateField(fieldDefn)

    # loop through the input features
    m_feature = m_layer.GetNextFeature()
    while m_feature:
        # print(m_feature)
        o_geometry = m_feature.GetGeometryRef()
        # 关键，合并几何
        o_geometry = o_geometry.Union(o_geometry)
        outfeature = ogr.Feature(outLayerDefn)
        outfeature.SetGeometry(o_geometry)
        # 遍历每个要素的字段，并设置字段属性
        for i in range(0, outLayerDefn.GetFieldCount()):
            # print(outLayerDefn.GetFieldDefn(i).GetNameRef())
            outfeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), m_feature.GetField(i))
        outlayer.CreateFeature(outfeature)
        # dereference the features and get the next input feature
        outfeature = None
        m_feature = m_layer.GetNextFeature()


    outds.Destroy()
def Check_self_intersect(Basin_shp_venu):
    """
    检查矢量文件是否存在空间自相交的问题，并解决
    :param Basin_shp_venu:
    :return:
    """

    Basin_list=os.listdir(Basin_shp_venu)
    shp_list=[file for file in Basin_list if file.endswith(".shp")]
    del Basin_list

    temp=os.path.join(os.path.dirname(Basin_shp_venu),"temp_1")
    if not os.path.exists(temp):
        os.mkdir(temp)
    for file in shp_list:
        save_path=os.path.join(temp,file)
        s_shp_si(os.path.join(Basin_shp_venu,file),save_path)

    shutil.rmtree(Basin_shp_venu)
    os.rename(temp,Basin_shp_venu)

def BB1(venu,outfile):
    """
    生产update文件，将所有的内流流域改为外流
    :param venu:
    :return:
    """
    if not os.path.exists(venu):
        return
    result = [['record','newtype']]
    dirNames = os.listdir(venu)
    for dirName in dirNames:
        if dirName.split('.')[-1] != 'tif':
            continue
        result.append([dirName.split('_')[0],3])

    # outfile = r'F:\青藏高原水体数据集\DATA\visiualCheck\Greendland\Greenland_9020000010_9030013780_update_.csv'
    with open(outfile,'w',newline='') as f:
        writer = csv.writer(f)
        writer.writerows(result)
        f.close()
def sbatch_BB1(venu):
    dirNames = os.listdir(venu)
    for dirName in dirNames:
        if len(dirName.split('.')) > 1:
            continue
        if len(dirName.split('_')) == 0:
            continue
        maskpath = os.path.join(venu,dirName,'mask')
        outfile = os.path.join(venu,dirName,dirName+'_update_.csv')
        if os.path.join(maskpath):
                BB1(maskpath,outfile)


def spatial_DEM_fdr(ref_path,src_path,out_path):
    """
    根据fdr裁剪DEM，保持像元对齐
    :param ref_path:
    :param src_path:
    :param out_path:
    :return:
    """
    import rasterio
    from rasterio.warp import reproject, Resampling

    # 输入文件
    # ref_path = "reference.tif"  # 参考栅格
    # src_path = "to_align.tif"  # 待对齐栅格
    # out_path = "aligned.tif"  # 输出文件

    # 打开参考与待对齐栅格
    # 打开参考和源文件
    with rasterio.open(ref_path) as ref, rasterio.open(src_path) as src:
        kwargs = src.meta.copy()
        kwargs.update({
            'height': ref.height,
            'width': ref.width,
            'transform': ref.transform,
            'compress': 'lzw',  # 启用压缩，可改为 'deflate'
            'tiled': True  # 建议加快读取速度
        })

        # 对齐并保存
        with rasterio.open(out_path, 'w', **kwargs) as dst:
            reproject(
                source=rasterio.band(src, 1),
                destination=rasterio.band(dst, 1),
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=ref.transform,
                dst_crs=ref.crs,
                resampling=Resampling.nearest
            )






if __name__ == '__main__':


    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030034520")
    #NA
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\NA')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica","/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/NA")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica")
    # Merge.merge_continent2("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica")
    # Find.draw_tif("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_dir.tif","/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_dir.png")
    # mask_extent_sink("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020000010")
    # Merge.clip_modifiedDir('/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020000010/010.shp',"/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020000010/SouthAmerica_6020000010_ModifiedDir_.tif",
    #                        "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020000010")
    Merge.clip_modifiedDir('/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020006540/540.shp',
                           "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020006540/SouthAmerica_6020006540_ModifiedDir_.tif",
                           "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020006540")
    # Merge.clip_modifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020008320/320.shp",
    #                        "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020008320/SouthAmerica_6020008320_ModifiedDir_.tif",
    #                        "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020008320")
    # Merge.clip_modifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020014330/330.shp",
    #                        "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020014330/SouthAmerica_6020014330_ModifiedDir_.tif",
    #                        "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020014330")
    # Merge.clip_modifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020017370/370.shp",
    #                        "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020017370/SouthAmerica_6020017370_ModifiedDir_.tif",
    #                        "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020017370")
    # Merge.clip_modifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020021870/870.shp",
    #                        "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020021870/SouthAmerica_6020021870_ModifiedDir_.tif",
    #                        "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020021870")
    # Merge.clip_modifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020029280/280.shp",
    #                        "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020029280/SouthAmerica_6020029280_ModifiedDir_.tif",
    #                        "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020029280")
    # valid.get_basin_main('/datanode05/zhangbin/NWEI/data/GRDC',"/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_fdir.tif","/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/valid/")
    # valid.sbatch_cal_area("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/valid",
    #                       '/datanode05/zhangbin/NWEI/data/GRDC')



    # Merge.get_sink_mask("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_6020000010")
    # valid.get_basin_main('/datanode05/zhangbin/NWEI/data/GRDC',"/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_fdir.tif",'/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/valid')
    # valid.sbatch_cal_area("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/valid",'/datanode05/zhangbin/NWEI/data/GRDC')
    # LZW("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030034520_1/NorthAmerica_7030034520_ModifiedDir_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030034520_1/NorthAmerica_7030034520_ModifiedDir_1.tif")
    #
    # sink.Cal_acc("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030034520_1/NorthAmerica_7030034520_ModifiedDir_1.tif",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030034520_1/NorthAmerica_7030034520_ModifiedDir_acc.tif")
    #
    # sink.Cal_Stream("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030034520_1/NorthAmerica_7030034520_ModifiedDir_acc.tif",
    #                 "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030034520_1/NorthAmerica_7030034520_ModifiedDir_10000.tif",10000)

    # Asia
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Asia/230","/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4030050230/Asia_4030050230_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4030050230/Asia_4030050230_Occ_.tif","/datanode05/zhangbin/TBP_Stream/DATA/Asia")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Asia/230","/datanode05/zhangbin/TBP_Stream/DATA/Asia")
    PFD8("/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4030050230/Asia_4030050230_burnDEM_.tif",
        "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4030050230/Asia_4030050230_ModifiedDir_.tif")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\AS')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Asia",
    #                         "/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/AS")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Asia")
    # valid.get_basin_main('/datanode05/zhangbin/NWEI/data/GRDC',
    #                      "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_fdir.tif",
    #                      "/datanode05/zhangbin/TBP_Stream/DATA/Asia/valid")
    # valid.sbatch_cal_area("/datanode05/zhangbin/TBP_Stream/DATA/Asia/valid",
    #                       '/datanode05/zhangbin/NWEI/data/GRDC')

    # NorthAmerica
    # PFD8("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030034520_1/NorthAmerica_7030034520_burnDEM_.tif",
    #    "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030034520_1/NorthAmerica_7030034520_ModifiedDir_.tif")
    # PFD8("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030022240_1/NorthAmerica_7030022240_burnDEM_.tif",
    #      "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030022240_1/NorthAmerica_7030022240_burnDEM_ModifiedDir_.tif")
    # LZW("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030034520_1/NorthAmerica_7030034520_ModifiedDir_1.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030034520_1/NorthAmerica_7030034520_ModifiedDir_.tif")

    # Find_null_sink("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/NorthAmerica/6750",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020046750/NorthAmerica_7020046750_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020046750/NorthAmerica_7020046750_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/NorthAmerica/6750",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")

    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/NorthAmerica/1430",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020021430/NorthAmerica_7020021430_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020021430/NorthAmerica_7020021430_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/NorthAmerica/1430",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")
    #
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/NorthAmerica/4600",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020024600/NorthAmerica_7020024600_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020024600/NorthAmerica_7020024600_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/NorthAmerica/4600",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/NorthAmerica/010",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020000010/NorthAmerica_7020000010_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7020000010/NorthAmerica_7020000010_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/NorthAmerica/010",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/NorthAmerica/1430_1",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030021430/NorthAmerica_7030021430_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030021430/NorthAmerica_7030021430_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/NorthAmerica/1430_1",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/NorthAmerica/4520",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030034520/NorthAmerica_7030034520_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030034520/NorthAmerica_7030034520_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/NorthAmerica/4520",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/NorthAmerica/2240",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030022240/NorthAmerica_7030022240_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030022240/NorthAmerica_7030022240_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/NorthAmerica/2240",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")

    # Africa
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Africa/0010",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020000010/Africa_1020000010_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020000010/Africa_1020000010_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Africa/0010",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Africa")

    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Africa/1530",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020011530/Africa_1020011530_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020011530/Africa_1020011530_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Africa/1530",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Africa")
    #
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Africa/4170",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170/Africa_1020034170_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170/Africa_1020034170_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Africa/4170",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Africa")
    #
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Africa/8110",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110/Africa_1020018110_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110/Africa_1020018110_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Africa/8110",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Africa")

    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Africa/1660",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1030011660/Africa_1030011660_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1030011660/Africa_1030011660_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Africa/1660",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Africa")

    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Africa/9940",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1041259940/Africa_1041259940_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1041259940/Africa_1041259940_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Africa/9940",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Africa")

    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Africa/3040",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1040843040/Africa_1040843040_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1040843040/Africa_1040843040_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Africa/3040",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Africa")
    # Find_null_sink("/datanode05/zhangbin/TBP_Stream/DATA/Africa")
    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020000010_sink")
    process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020000010_sink")
    postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020000010_sink")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\Af\0010')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020000010_sink","/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/Af/0010")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020000010_sink")

    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020011530_1/Africa_1030011660_sink")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020011530_1/Africa_1030011660_sink")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020011530_1/Africa_1030011660_sink")
    # # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\Af\1660')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020011530_1/Africa_1030011660_sink",
    #                         "/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/Af/1660")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020011530_1/Africa_1030011660_sink")

    calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170_1/Africa_1040843040_sink")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170_1/Africa_1040843040_sink")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170_1/Africa_1040843040_sink")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\Af\3040')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170_1/Africa_1040843040_sink",
    # "/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/Af/3040")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170_1/Africa_1040843040_sink")

    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170_1")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170_1")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170_1")
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170_1","/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/Af/4170")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170_1")


    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110_1")
    process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110_1")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110_1")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\Af\8110')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110_1",
    #                         "/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/Af/8110")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110_1")



    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110_1/Africa_1041259940_sink")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110_1/Africa_1041259940_sink")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110_1/Africa_1041259940_sink")
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110_1/Africa_1041259940_sink","/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/Af/9940")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110_1/Africa_1041259940_sink")



    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020011530_1")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020011530_1")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020011530_1")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\Af\1530')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020011530_1",
    #                         "/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/Af/1530")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020011530_1")

    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\Africa\7430')
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Africa/7430_1",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430.1/Africa_1020027430_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430.1/Africa_1020027430_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430")

    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Africa/7430_1","/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430")
    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430")
    # calculate_sink_area_volume.calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_1030040310")
    # calculate_sink_area_volume.calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_1030040300")
    # calculate_sink_area_volume.calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_1030040280")
    #
    # calculate_sink_area_volume.calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_1030040250")
    # calculate_sink_area_volume.calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_1030040220")
    # calculate_sink_area_volume.calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_1030034610")
    #
    # calculate_sink_area_volume.calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_1030031860")
    # calculate_sink_area_volume.calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_1030029810")
    # calculate_sink_area_volume.calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_1030027430")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430")


    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170/Africa_1040843040/Africa_1051159450")
    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170/Africa_1040843040/Africa_1051089880")
    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020011530/Africa_1030011660/Africa_1041472390")
    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110/Africa_1041259940/Africa_1051265610")

    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa")





    # Europe
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/1190",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190/Europe_2020071190_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190/Europe_2020071190_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/1190",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Europe")
    #
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/5840",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840/Europe_2020065840_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840/Europe_2020065840_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/5840",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Europe")
    #
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/4230",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020024230/Europe_2020024230_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020024230/Europe_2020024230_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/4230",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Europe")

    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/2030065840",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030065840/Europe_2030065840_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030065840/Europe_2030065840_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/2030065840",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Europe")

    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/2030066850",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030066850/Europe_2030066850_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030066850/Europe_2030066850_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/2030066850",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Europe")


    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/2050066490",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2050066490/Europe_2050066490_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2050066490/Europe_2050066490_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2050066490")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/2050066490","/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2050066490")
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/2040067740",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2040067740/Europe_2040067740_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2040067740/Europe_2040067740_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2040067740")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/2040067740","/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2040067740")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\Europe\2040067740')
    def AAA(dem_file,outmask,outdem):
        dem = Raster.get_raster(dem_file)
        proj,geo,nodata = Raster.get_proj_geo_nodata(dem_file)

        row,col = dem.shape
        mask = np.zeros((row,col))
        for i in range(row):
            for j in range(col):
                if dem[i,j] != nodata:
                    mask[i,j] = 1
    #
        Raster.save_raster(outmask,mask,proj,geo,gdal.GDT_Byte,0)
        # Raster.save_raster(outdem,dem,proj,geo,gdal.GDT_Float32,nodata)
    # AAA("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2050066500/Europe_2050066500_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2050066500/Europe_2050066500_mask_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2050066500/Europe_2050066500_dem_.tif")
    # sink.Cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2050066500/Europe_2050066500_dem_.tif","/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2050066500/Europe_2050066500_dir_.tif")
    # endorheic_process.process_endorheic("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2050066500")

    # AAA("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2040067280/Europe_2040067280_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2040067280/Europe_2040067280_mask_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2040067280/Europe_2040067280_dem_.tif")
    # sink.Cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2040067280/Europe_2040067280_dem_.tif","/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2040067280/Europe_2040067280_dir_.tif")
    # endorheic_process.process_endorheic("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2040067280")

    # AAA("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2040085880/Europe_2040085880_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2040085880/Europe_2040085880_mask_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2040085880/Europe_2040085880_dem_.tif")
    # sink.Cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2040085880/Europe_2040085880_dem_.tif",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2040085880/Europe_2040085880_dir_.tif")
    # endorheic_process.process_endorheic("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2040085880")

    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/2030068690",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030068690_2/Europe_2030068690_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030068690_2/Europe_2030068690_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030068690_2")
    # AAA("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020003440/Europe_2020003440_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020003440/Europe_2020003440_mask_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030068690_2/Europe_2030068690_2_2040085690/Europe_2030068690_2_2040085690_dem_.tif")
    # sink.Cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030068690_2/Europe_2030068690_2_2040085690/Europe_2030068690_2_2040085690_dem_.tif",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030068690_2/Europe_2030068690_2_2040085690/Europe_2030068690_2_2040085690_dir_.tif")
    # endorheic_process.process_endorheic("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030068690_2/Europe_2030068690_2_2040085690")
    # import calculate_sink_area_volume
    # try:
    #     calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2050066490")
    # except Exception as e:
    #     # 捕获并输出错误信息
    #     print(f"An error occurred: {e}")
    import process_sink
    # try:
    #     process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2050066490")
    # except Exception as e:
    #     # 捕获并输出错误信息
    #     print(f"An error occurred: {e}")

    # try:
    #     postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2050066490")
    # except Exception as e:
    #     # 捕获并输出错误信息
    #     print(f"An error occurred: {e}")

    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\EU\1190')
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190")
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190","/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/EU/1190")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190")


    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030066850_sink/Europe_2040067740_visual")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030066850_sink/Europe_2040067740_visual")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\EU\2040067740')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030066850_sink/Europe_2040067740_visual",
    #                         "/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/EU/2040067740")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030066850_sink/Europe_2040067740_visual")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030066850_sink/Europe_2040067740_visual")

    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\EU\2050066490')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030065840_sink/Europe_2050066490",
    #                         "/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/EU/2050066490")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030065840_sink/Europe_2050066490")
    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030065840_sink")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030065840_sink")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030065840_sink")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\EU\2030065840')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030065840_sink","/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/EU/2030065840")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030065840_sink")

    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/2040067570","/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030066850_sink/Europe_2040067570/Europe_2040067570_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030066850_sink/Europe_2040067570/Europe_2040067570_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030066850_sink/Europe_2040067570")
    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030066850_sink/Europe_2040067570")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030066850_sink/Europe_2040067570")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030066850_sink/Europe_2040067570")
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030066850_sink/Europe_2040067570","/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/EU/2040067570")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030066850_sink/Europe_2040067570")

    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030066850_sink")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030066850_v")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030066850_v")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\EU\2030066850')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030066850_v",
    #                         "/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/EU/2030066850")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink/Europe_2030066850_v")

    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020024230_sink")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020024230_sink")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\Europe\4230')
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020024230_sink")
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020024230_sink","/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/EU/4230")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020024230_sink")

    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\Europe\5840')
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840_sink")
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840","/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/EU/5840")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840")

    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030068690_2")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030068690_2")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\EU\2030068690')
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030068690_2")
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030068690_2","/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/EU/2030068690")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030068690_2")

    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2")

    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\EU\level2')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe","/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/EU/level2")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Europe")


    # valid.get_basin_main('/datanode05/zhangbin/NWEI/data/GRDC',
    #                      "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_fdir.tif",
    #                      "/datanode05/zhangbin/TBP_Stream/DATA/Europe/valid/")

    # valid.sbatch_cal_area("/datanode05/zhangbin/TBP_Stream/DATA/Europe/valid",
    #                       '/datanode05/zhangbin/NWEI/data/GRDC')

    def BB(venu):
        # files = []
        # dirNames = os.listdir(venu)
        # for dirName in dirNames:
        #     if dirName.split('.')[-1] != 'tif':
        #         continue
        #     files.append(os.path.join(venu,dirName))
        files = ["/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110_1/Africa_1041259940_sink/Africa_1051259200/Africa_1051259200_burnDEM_.tif",
                 "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110_1/Africa_1041259940_sink/Africa_1051265610/Africa_1051265610_burnDEM_.tif",
                 "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110_1/Africa_1041259940_sink/Africa_1051265730/Africa_1051265730_burnDEM_.tif"]
        gdal.Warp("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110_1/Africa_1041259940_sink/Africa_1051265610_burnDEM_.tif", files, options=[
            "-co", "COMPRESS=LZW",  # 选择LZW压缩
            "-co", "BIGTIFF=YES"])  # 允许大文件])
        print("合并完成，输出文件：")
    # BB("/datanode05/zhangbin/data/SRTM30DEM")
    # sink.Cal_sink("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030068690_2/fill.europe/Europe_fill.tif",
    #               "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030068690_2/fill.europe/europeFill_burnsink.tif")
    # sink.Cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030068690_2/fill.europe/Europe_fill.tif",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030068690_2/fill.europe/europeFill_burndir.tif")
    # import calculate_sink_area_volume
    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2030068690_2")





    # Siberia:未处理
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Siberia/9230",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Siberia/9230",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320")
    #
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Siberia/8670",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670/Siberia_3020008670_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670/Siberia_3020008670_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Siberia/8670",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670")
    #
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Siberia/3020003790",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020003790/Siberia_3020003790_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020003790/Siberia_3020003790_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020003790")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Siberia/3020003790",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020003790")
    # valid.sbatch_cal_area("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/valid",
    #                       '/datanode05/zhangbin/NWEI/data/GRDC')


    # try:
    #     calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670")
    # except Exception as e:
    #     # 捕获并输出错误信息
    #     print(f"An error occurred: {e}")
    # try:
    #     process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670")
    # except Exception as e:
    #     # 捕获并输出错误信息
    #     print(f"An error occurred: {e}")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\Si\3020008670')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670","/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/Si/3020008670")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670")
    # try:
    #     calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320")
    # except Exception as e:
    #     # 捕获并输出错误信息
    #     print(f"An error occurred: {e}")
    # try:
    #     process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320")
    # except Exception as e:
    #     # 捕获并输出错误信息
    #     print(f"An error occurred: {e}")
    # # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\Si\Siberia_3020009320')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320",
    #                         "/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/Si/Siberia_3020009320")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320")

    # try:
    #     calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020003790")
    # except Exception as e:
    #     # 捕获并输出错误信息
    #     print(f"An error occurred: {e}")
    # try:
    #     process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020003790")
    # except Exception as e:
    #     # 捕获并输出错误信息
    #     print(f"An error occurred: {e}")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\Si\Siberia_3020003790')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020003790",
    #                         "/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/Si/Siberia_3020003790")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020003790")


    # try:
    #     postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320")
    # except Exception as e:
    #     # 捕获并输出错误信息
    #     print(f"An error occurred: {e}")

    # try:
    #     postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670")
    # except Exception as e:
    #     # 捕获并输出错误信息
    #     print(f"An error occurred: {e}")

    # try:
    #     postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020003790")
    # except Exception as e:
    #     # 捕获并输出错误信息
    #     print(f"An error occurred: {e}")

    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\Si\level2')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia",
    #                         "/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/Si/level2")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia")
    # valid.get_basin_main('/datanode05/zhangbin/NWEI/data/GRDC',"/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_fdir.tif",
    #                      "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/valid")



    # Greenland
    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Greenland")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\GreenLand\0010')
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Greenland/0010",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Greenland/Greenland_9020000010.1/Greenland_9020000010_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Greenland/Greenland_9020000010.1/Greenland_9020000010_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Greenland/Greenland_9020000010")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Greenland/0010","/datanode05/zhangbin/TBP_Stream/DATA/Greenland/Greenland_9020000010")
    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Greenland/Greenland_9020000010")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Greenland/Greenland_9020000010")
    # sbatch_BB1("/datanode05/zhangbin/TBP_Stream/DATA/Greenland/Greenland_9020000010")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Greenland/Greenland_9020000010")






    # AU

    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\Au\9720')
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Australia/9720/",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_5020049720.1/Australia_5020049720_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_5020049720.1/Australia_5020049720_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_5020049720")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Australia/9720","/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_5020049720")
    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_5020049720")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_5020049720")
    # calculate_sink_area_volume.calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_5020049720/Australia_5020049720_5030079570")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_5020049720")


    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_5020055870/Australia_5020055870_dem_.tif",
    #           "/datanode05/zhangbin/data/5020055870.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_5020055870/Australia_5020055870_burnDEM_.tif")
    # sink.Cal_sink("/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_5020055870/Australia_5020055870_burnDEM_.tif",
    #               "/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_5020055870/Australia_5020055870_burnsink_.tif")
    # sink.Cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_5020055870/Australia_5020055870_burnDEM_.tif",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_5020055870/Australia_5020055870_burndir_.tif")
    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Australia")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Australia")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Australia")

    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_5020049720",
    #                         "/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/AU/9720")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\Au\Australia')
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_5020049720")
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Australia",
    #                         "/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/AU/Australia")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Australia")

    # valid.get_basin_main(r'',"/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_fdir.tif","/datanode05/zhangbin/TBP_Stream/DATA/Australia/valid")
    # valid.sbatch_cal_area("/datanode05/zhangbin/TBP_Stream/DATA/Australia/valid")




    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_fdir.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/valid/mssi/mssi.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/valid/mssi/mssi_fdir.tif")
    # sink.Cal_acc("/datanode05/zhangbin/TBP_Stream/DATA/valid/mssi/mssi_fdir.tif",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/valid/mssi/mssi_acc.tif")
    # sink.Cal_Stream("/datanode05/zhangbin/TBP_Stream/DATA/valid/mssi/mssi_acc.tif",
    #                 "/datanode05/zhangbin/TBP_Stream/DATA/valid/mssi/mssi_stream_5km2.tif",5555)
    # sink.streamFeature("/datanode05/zhangbin/TBP_Stream/DATA/valid/mssi/mssi_fdir.tif",
    #                    "/datanode05/zhangbin/TBP_Stream/DATA/valid/mssi/mssi_stream_5km2.tif",
    #                    "/datanode05/zhangbin/TBP_Stream/DATA/valid/mssi/mssi_stream_5km2.shp")


    def Check(venu):
        result = []
        dirNames = os.listdir(venu)
        popPaths = [os.path.join(venu,dirName) for dirName in dirNames]
        while popPaths:
            popPath = popPaths.pop()

            dirName = os.path.basename(popPath)
            if len(dirName.split('.')) > 1:
                continue
            if len(dirName.split('_')) == 0:
                continue
            if os.path.exists(os.path.join(popPath,dirName+'_ModifiedDir_.tif')):
                read_and_save_tif(os.path.join(popPath,dirName+'_ModifiedDir_.tif'),os.path.join(popPath,dirName+'_ModifiedDir_.png'))
                os.chmod(os.path.join(popPath,dirName+'_ModifiedDir_.png'),0o777)

            sink_file = os.path.join(popPath,dirName+'_burnsink_.tif')
            fdir_file = os.path.join(popPath,dirName+'_ModifiedDir_.tif')
            if (not os.path.exists(sink_file)) and (os.path.exists(fdir_file)):
                result.append([popPath])
            for temp_dirName in os.listdir(popPath):
                if len(temp_dirName.split('.')) > 1:
                    continue
                if len(temp_dirName.split('_')) == 0:
                    continue
                popPaths.append(os.path.join(popPath,temp_dirName))

        with open("/datanode05/zhangbin/TBP_Stream/DATA/CHECK.csv",'w',newline='') as f:
            writer = csv.writer(f)
            writer.writerows(result)
            f.close()
        os.chmod("/datanode05/zhangbin/TBP_Stream/DATA/CHECK.csv",0o777)




    # def temp_A(venu):
    #     pathNames = os.listdir(venu)
    #     for pathName in pathNames:
    #         if len(pathName.split('.')) != 1:
    #             continue
    #         file = os.path.join(venu,pathName,pathName+'_burnDEM_.tif')
    #         outfile = os.path.join(venu,pathName,pathName+'_burnDEM_.png')
    #         read_and_save_tif(file,outfile)
    #
    # temp_A("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430")

    # Check("/datanode05/zhangbin/TBP_Stream/DATA/Africa")
    # Check("/datanode05/zhangbin/TBP_Stream/DATA/Europe")


    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_fdir.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_fdir.png")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030022240/NorthAmerica_7030022240_burnDEM_.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030022240/NorthAmerica_7030022240_burnDEM_.png")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7040569640/NorthAmerica_7040569640_burnDEM_.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7040569640/NorthAmerica_7040569640_burnDEM_.png")
    # Merge.merge_DEM("/datanode05/zhangbin/FABDEM_preprocess",[-120,58,-87,36],"/datanode05/zhangbin/TBP_Stream/DATA/regionDem/NA_bu.tif")



    # BB1(r'F:\青藏高原水体数据集\DATA\visiualCheck\Greendland\maskGreenland_9020000010_9030013780')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Greenland/Greenland_9020000010","/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/Gr/0010")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Greenland/Greenland_9020000010")

    # Arctic
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\Arctic\8900')
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Arctic/8900",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_2/Arctic_8020008900_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_2/Arctic_8020008900_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Arctic/8900","/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1")
    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1")
    # sbatch_BB1("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1")

    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1/Arctic_8020008900_1_8040307300")
    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1/Arctic_8020008900_1_8040264670")
    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1/Arctic_8020008900_1_8040252610")
    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1/Arctic_8020008900_1_8040183070")


    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020020760/Arctic_8020020760_burnDEM_.tif", "/datanode05/zhangbin/TBP_Stream/DATA/mask/Ar/8020020760.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020020760/Arctic_8020020760_burnDEM_1.tif")
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020020760/Arctic_8020020760_Occ_.tif", "/datanode05/zhangbin/TBP_Stream/DATA/mask/Ar/8020020760.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020020760/Arctic_8020020760_Occ_1.tif")
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020020760/Arctic_8020020760_burnsink_.tif","/datanode05/zhangbin/TBP_Stream/DATA/mask/Ar/8020020760.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020020760/Arctic_8020020760_burnsink_1.tif")
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020020760/Arctic_8020020760_burndir_2.tif","/datanode05/zhangbin/TBP_Stream/DATA/mask/Ar/8020020760.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020020760/Arctic_8020020760_burndir_.tif")
    # calculate_sink_area_volume.calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020020760")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Arctic")
    # BB1("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020020760","/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020020760/Arctic_8020020760_update_.csv")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Arctic")

    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\Arctic\0010')
    # preprocess("/datanode05/zhangbin/TBP_Stream/DATA/mask/Arctic/0010",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020000010_1/Arctic_8020000010_burnDEM_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020000010_1/Arctic_8020000010_Occ_.tif",
    #            "/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020000010_1")
    # cal_dir_sink("/datanode05/zhangbin/TBP_Stream/DATA/mask/Arctic/0010","/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020000010_1")
    # calculate_sink_area_volume.sbatch_calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020000010_1")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020000010_1")
    # sbatch_BB1("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020000010_1")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020000010_1")

    # sbatch_BB1("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1")

    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1/Arctic_8020008900_1_8040264670/Arctic_8020008900_1_8040264670_ModifiedDir_new_.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1/Arctic_8020008900_1_8040264670/Arctic_8020008900_1_8040264670_ModifiedDir_new_.png")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1/Arctic_8020008900_1_8040307300/Arctic_8020008900_1_8040307300_ModifiedDir_new_.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1/Arctic_8020008900_1_8040307300/Arctic_8020008900_1_8040307300_ModifiedDir_new_.png")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1/Arctic_8020008900_1_8040183070/Arctic_8020008900_1_8040183070_ModifiedDir_new_.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1/Arctic_8020008900_1_8040183070/Arctic_8020008900_1_8040183070_ModifiedDir_new_.png")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1/Arctic_8020008900_1_8040252610/Arctic_8020008900_1_8040252610_ModifiedDir_new_.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1/Arctic_8020008900_1_8040252610/Arctic_8020008900_1_8040252610_ModifiedDir_new_.png")

    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020000010_1",
    #                         "/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/Ar/0010")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\Ar\0010')
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\Ar\8900')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1",
    #                         "/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/Ar/8900")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020000010_1")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\Ar\Arctic')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Arctic","/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/Ar/Arctic")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Arctic")

    # valid.get_basin_main(r'',"/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_fdir.tif","/datanode05/zhangbin/TBP_Stream/DATA/Arctic/valid")


    # Africa补

    # Merge.merge_DEM("/datanode05/zhangbin/FABDEM_preprocess", [-12.8, 34.1, 12.6, 14.9],
    #                 "/datanode05/zhangbin/TBP_Stream/DATA/regionDem/Af3_0310.tif")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/Af3_0310.tif","/datanode05/zhangbin/TBP_Stream/DATA/regionDem/Af3_0310.png")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_burnDEM_.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_burnDEM_.png")
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/Af3_0310.tif", "/datanode05/zhangbin/TBP_Stream/DATA/mask/Africa/FIll/3_0310_buffer.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_1030040310/Africa_1020027430_1030040310_burnDEM_.tif")
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_Occ_.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/mask/Africa/FIll/3_0310_buffer.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_1030040310/Africa_1020027430_1030040310_Occ_.tif")
    # read_and_save_tif(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_1030040310/Africa_1020027430_1030040310_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_1030040310/Africa_1020027430_1030040310_burnDEM_.png")
    # sink.Cal_sink("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_1030040310/Africa_1020027430_1030040310_burnDEM_.tif",
    #               "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_1030040310/Africa_1020027430_1030040310_burnsink_.tif")
    # sink.Cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_1030040310/Africa_1020027430_1030040310_burnDEM_.tif",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_1030040310/Africa_1020027430_1030040310_burndir_.tif")
    # calculate_sink_area_volume.calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_1030040310")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\Africa\7430')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430",
    #                         "/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/Af/7430")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\Af\Africa')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa","/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/Af/Africa")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Africa")


    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_1030040310/Africa_1020027430_1030040310_Occ_.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430/Africa_1020027430_1030040310/Africa_1020027430_1030040310_Occ_.png")
    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110/Africa_1041259940/Africa_1051265610")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110/Africa_1041259940/Africa_1051265610/Africa_1051265610_ModifiedDir_new_.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110/Africa_1041259940/Africa_1051265610/Africa_1051265610_ModifiedDir_new_.png")
    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020011530/Africa_1030011660/Africa_1041472390")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020011530/Africa_1030011660/Africa_1041472390/Africa_1041472390_ModifiedDir_new_.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020011530/Africa_1030011660/Africa_1041472390/Africa_1041472390_ModifiedDir_new_.png")
    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170/Africa_1040843040/Africa_1051159450")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170/Africa_1040843040/Africa_1051089880/Africa_1051089880_ModifiedDir_new_.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170/Africa_1040843040/Africa_1051089880/Africa_1051089880_ModifiedDir_new_.png")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170/Africa_1040843040/Africa_1051159450/Africa_1051159450_ModifiedDir_new_.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170/Africa_1040843040/Africa_1051159450/Africa_1051159450_ModifiedDir_new_.png")
    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170/Africa_1040843040/Africa_1051089880")

    # valid.get_basin_main(r'',"/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_fdir.tif","/datanode05/zhangbin/TBP_Stream/DATA/Africa/valid")
    # valid.sbatch_cal_area("/datanode05/zhangbin/TBP_Stream/DATA/Africa/valid")


    def A(venuV, outVenu):
        """
        批量处理文件夹中的DEM:优先搜索边界上的0值加入队列，使用图搜索将8邻域内的0值栅格赋值为nodata
        :param venu:
        :param outVenu:
        :return:
        """
        # if not os.path.exists(outVenu):
        #     os.mkdir(outVenu)
        def C(demFile,outFile):
            shutil.copyfile(demFile,outFile)
        venus = os.listdir(venuV)
        for venu1 in venus:
            venu = os.path.join(venuV, venu1)
            files = os.listdir(venu)
            for file in files:
                if file.split('.')[-1] != 'tif':
                    continue
                demFile = os.path.join(venu, file)
                outFile = os.path.join(outVenu, file)
                C(demFile, outFile)

    # SI 补
    # A("/datanode05/zhangbin/FABDEM","/datanode05/zhangbin/FABDEM_flat")
    # Merge.merge_DEM("/datanode05/zhangbin/FABDEM_flat",[143,72,180,59],"/datanode05/zhangbin/TBP_Stream/DATA/regionDem/SI_143_72_180_59.tif")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/SI_143_72_180_59.tif","/datanode05/zhangbin/TBP_Stream/DATA/regionDem/SI_143_72_180_59.png")

    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/Siberia.tif","/datanode05/zhangbin/TBP_Stream/DATA/regionDem/Siberia.png")
    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020003790/Siberia_3020003790_3040483880")

    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670/Siberia_3020008670_Occ_.tif","/datanode05/zhangbin/TBP_Stream/DATA/mask/Siberia/FILL/9300.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670/Siberia_3020008670_3040299300/Siberia_3020008670_3040299300_Occ_.tif")
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/Siberia.tif","/datanode05/zhangbin/TBP_Stream/DATA/mask/Siberia/FILL/9300.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670/Siberia_3020008670_3040299300/Siberia_3020008670_3040299300_burnDEM_.tif")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670/Siberia_3020008670_3040299300/Siberia_3020008670_3040299300_burnDEM_.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670/Siberia_3020008670_3040299300/Siberia_3020008670_3040299300_burnDEM_.png")
    # sink.Cal_sink("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670/Siberia_3020008670_3040299300/Siberia_3020008670_3040299300_burnDEM_.tif",
    #               "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670/Siberia_3020008670_3040299300/Siberia_3020008670_3040299300_burnsink_.tif")
    # sink.Cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670/Siberia_3020008670_3040299300/Siberia_3020008670_3040299300_burnDEM_.tif",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670/Siberia_3020008670_3040299300/Siberia_3020008670_3040299300_burndir_.tif")
    # calculate_sink_area_volume.calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670/Siberia_3020008670_3040299300")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670")
    # BB1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670/Siberia_3020008670_3040299300/mask",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670/Siberia_3020008670_3040299300/Siberia_3020008670_3040299300_update_.csv")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670")

    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670/Siberia_3020008670_Occ_.tif","/datanode05/zhangbin/TBP_Stream/DATA/mask/Siberia/FILL/9300.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670/Siberia_3020008670_3040299300/Siberia_3020008670_3040299300_Occ_.tif")
    #
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/SI_143_72_180_59.tif","/datanode05/zhangbin/TBP_Stream/DATA/mask/Siberia/FILL/0910.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040010910/Siberia_3020009320_3040010910_burnDEM_.tif")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040010910/Siberia_3020009320_3040010910_burnDEM_.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040010910/Siberia_3020009320_3040010910_burnDEM_.png")
    # sink.Cal_sink("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040010910/Siberia_3020009320_3040010910_burnDEM_.tif",
    #               "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040010910/Siberia_3020009320_3040010910_burnsink_.tif")
    # sink.Cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040010910/Siberia_3020009320_3040010910_burnDEM_.tif",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040010910/Siberia_3020009320_3040010910_burndir_.tif")
    # calculate_sink_area_volume.calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040010910")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320")
    # BB1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040010910/mask",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040010910/Siberia_3020009320_3040010910_update_.csv")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320")


    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/SI_143_72_180_59.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/mask/Siberia/FILL/0920.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040010920/Siberia_3020009320_3040010920_burnDEM_.tif")
    # read_and_save_tif(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040010920/Siberia_3020009320_3040010920_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040010920/Siberia_3020009320_3040010920_burnDEM_.png")
    # sink.Cal_sink(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040010920/Siberia_3020009320_3040010920_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040010920/Siberia_3020009320_3040010920_burnsink_.tif")
    # sink.Cal_dir(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040010920/Siberia_3020009320_3040010920_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040010920/Siberia_3020009320_3040010920_burndir_.tif")
    # calculate_sink_area_volume.calculate_sink_volume_byDir(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040010920")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320")
    # BB1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040010920/mask",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040010920/Siberia_3020009320_3040010920_update.csv")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320")
    #
    #
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/SI_143_72_180_59.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/mask/Siberia/FILL/7570.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040017570/Siberia_3020009320_3040017570_burnDEM_.tif")
    # read_and_save_tif(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040017570/Siberia_3020009320_3040017570_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040017570/Siberia_3020009320_3040017570_burnDEM_.png")
    # sink.Cal_sink(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040017570/Siberia_3020009320_3040017570_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040017570/Siberia_3020009320_3040017570_burnsink_.tif")
    # sink.Cal_dir(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040017570/Siberia_3020009320_3040017570_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040017570/Siberia_3020009320_3040017570_burndir_.tif")
    # calculate_sink_area_volume.calculate_sink_volume_byDir(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040017570")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320")
    # BB1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040017570/mask",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040017570/Siberia_3020009320_3040017570_update.csv")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320")

    #
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/SI_143_72_180_59.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/mask/Siberia/FILL/2270.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040012270/Siberia_3020009320_3040012270_burnDEM_.tif")
    # read_and_save_tif(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040012270/Siberia_3020009320_3040012270_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040012270/Siberia_3020009320_3040012270_burnDEM_.png")
    # sink.Cal_sink(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040012270/Siberia_3020009320_3040012270_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040012270/Siberia_3020009320_3040012270_burnsink_.tif")
    # sink.Cal_dir(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040012270/Siberia_3020009320_3040012270_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040012270/Siberia_3020009320_3040012270_burndir_.tif")
    # calculate_sink_area_volume.calculate_sink_volume_byDir(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040012270")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670")
    # BB1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040012270/mask",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040012270/Siberia_3020009320_3040012270_update_.csv")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670")
    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040012270")
    #
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/SI_143_72_180_59.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/mask/Siberia/FILL/1770.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040011770/Siberia_3020009320_3040011770_burnDEM_.tif")
    # read_and_save_tif(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040011770/Siberia_3020009320_3040011770_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040011770/Siberia_3020009320_3040011770_burnDEM_.png")
    # sink.Cal_sink(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040011770/Siberia_3020009320_3040011770_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040011770/Siberia_3020009320_3040011770_burnsink_.tif")
    # sink.Cal_dir(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040011770/Siberia_3020009320_3040011770_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040011770/Siberia_3020009320_3040011770_burndir_.tif")
    # calculate_sink_area_volume.calculate_sink_volume_byDir(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040011770")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320")
    # BB1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040011770/mask",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320/Siberia_3020009320_3040011770/Siberia_3020009320_3040011770_update_.csv")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320")
    # Check_self_intersect(r'F:\青藏高原水体数据集\DATA\region_buffer\mask\Si\9230')
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320",
    #                         "/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/Si/9230/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia")


    # NA补
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/NA_bu.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/mask/NorthAmerica/FILL/9640.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7040569640/NorthAmerica_7040569640_burnDEM_.tif")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7040569640/NorthAmerica_7040569640_burnDEM_.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7040569640/NorthAmerica_7040569640_burnDEM_.png")
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7040569640.1/NorthAmerica_7040569640_Occ_.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/mask/NorthAmerica/FILL/9640.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7040569640/NorthAmerica_7040569640_Occ_.tif")
    # sink.Cal_sink("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7040569640/NorthAmerica_7040569640_burnDEM_.tif",
    #               "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7040569640/NorthAmerica_7040569640_burnsink_.tif")
    # sink.Cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7040569640/NorthAmerica_7040569640_burnDEM_.tif",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7040569640/NorthAmerica_7040569640_burndir_.tif")
    # calculate_sink_area_volume.calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7040569640")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")



    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/NA_bu.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/mask/NorthAmerica/FILL/2240.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030022240/NorthAmerica_7030022240_burnDEM_.tif")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030022240/NorthAmerica_7030022240_burnDEM_.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030022240/NorthAmerica_7030022240_burnDEM_.png")
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030022240.1/NorthAmerica_7030022240_Occ_.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/mask/NorthAmerica/FILL/2240.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030022240/NorthAmerica_7030022240_Occ_.tif")
    # sink.Cal_sink("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030022240/NorthAmerica_7030022240_burnDEM_.tif",
    #               "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030022240/NorthAmerica_7030022240_burnsink_.tif")
    # sink.Cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030022240/NorthAmerica_7030022240_burnDEM_.tif",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030022240/NorthAmerica_7030022240_burndir_.tif")
    # calculate_sink_area_volume.calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030022240")
    # process_sink.priority_D8_cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_7030022240")

    # EU补
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/Europe.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/regionDem/Europe.png")

    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/Europe.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/FILL/1190.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190/Europe_2030071190/Europe_2030071190_burnDEM_.tif")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190/Europe_2030071190/Europe_2030071190_burnDEM_.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190/Europe_2030071190/Europe_2030071190_burnDEM_.png")
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190/Europe_2020071190_Occ_.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/FILL/1190.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190/Europe_2030071190/Europe_2030071190_Occ_.tif")
    # sink.Cal_sink("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190/Europe_2030071190/Europe_2030071190_burnDEM_.tif",
    #               "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190/Europe_2030071190/Europe_2030071190_burnsink_.tif")
    # sink.Cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190/Europe_2030071190/Europe_2030071190_burnDEM_.tif",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190/Europe_2030071190/Europe_2030071190_burndir_.tif")
    # calculate_sink_area_volume.calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190/Europe_2030071190")
    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190")
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190",
    #                         "/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/EU/1190")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190")

    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/Europe.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/FILL/5840.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840/Europe_2020065840_burnDEM_.tif")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840/Europe_2020065840_burnDEM_.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840/Europe_2020065840_burnDEM_.png")
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840.1/Europe_2020065840_Occ_.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/FILL/5840.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840/Europe_2020065840_Occ_.tif")
    # sink.Cal_sink("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840/Europe_2020065840_burnDEM_.tif",
    #               "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840/Europe_2020065840_burnsink_.tif")
    # sink.Cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840/Europe_2020065840_burnDEM_.tif",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840/Europe_2020065840_burndir_.tif")
    # calculate_sink_area_volume.calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840")
    process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Europe")
    postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Europe")

    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840/Europe_2020065840_burnDEM_.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/5840/2030068680.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840/Europe_2030068680/Europe_2030068680_burnDEM_.tif")
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840/Europe_2020065840_Occ_.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/mask/Europe/5840/2030068680.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840/Europe_2030068680/Europe_2030068680_Occ_.tif")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840/Europe_2030068680/Europe_2030068680_burnDEM_.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840/Europe_2030068680/Europe_2030068680_burnDEM_.png")
    # sink.Cal_sink("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840/Europe_2030068680/Europe_2030068680_burnDEM_.tif",
    #               "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840/Europe_2030068680/Europe_2030068680_burnsink_.tif")
    # sink.Cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840/Europe_2030068680/Europe_2030068680_burnDEM_.tif",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020065840/Europe_2030068680/Europe_2030068680_burndir_.tif")

    process_sink.priority_D8_cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4030050230")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4030050230/Asia_4030050230_ModifiedDir_new_.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4030050230/Asia_4030050230_ModifiedDir_new_.png")

    def process_DEM(A,B):
        FAB = Raster.get_raster(A)
        SRTM = Raster.get_raster(B)

        proj,geo,F_nodata = Raster.get_proj_geo_nodata(A)
        proj,geo,SRTM_nodata = Raster.get_proj_geo_nodata(B)

        diffs = []
        for i in range(3600):
            for j in range(3600):
                if FAB[i,j] != F_nodata and SRTM[i,j] != SRTM_nodata:
                    diff = SRTM[i,j] - FAB[i,j]
                    diffs.append(diff)

        import matplotlib.pyplot as plt
        plt.scatter(np.arange(len(diffs)),diffs)
        plt.show()
    # process_DEM(r'F:\FABDEM\Global_FABDEM\N30E040-N40E050_FABDEM_V1-0\N34E042_FABDEM_V1-0.tif',
    #             r'F:\FABDEM\EU补SRTM\n34_e042_1arc_v3.tif')


    #######################################################  应用代码 #####################################################
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_fdir.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_fdir.png")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_fdir.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_fdir.png")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_fdir.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_fdir.png")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_fdir.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_fdir.png")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Greenland/Greenland_fdir.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Greenland/Greenland_fdir.png")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_fdir.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_fdir.png")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_fdir.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_fdir.png")


    # 识别内流区域
    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_fdir.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_endor.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_endor.tif")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_endor.tif","/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_endor.png")
    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_fdir.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_endor.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_endor.tif")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_endor.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_endor.png")
    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_fdir.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_endor.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_endor.tif")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_endor.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_endor.png")
    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_fdir.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_endor.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_endor.tif")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_endor.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_endor.png")
    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_fdir.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_endor.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_endor.tif")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_endor.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_endor.png")
    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_fdir.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_endor.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_endor.tif")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_endor.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_endor.png")
    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_fdir.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_endor.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_endor.tif")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_endor.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_endor.png")
    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_fdir.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_endor.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_endor.tif")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_endor.tif","/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_endor.png")


    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/BurnDATA/OSM/GlobalGSWO.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Peral_river/ZJ.shp","/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Peral_river/ZJ_OSM.tif")
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_fdir.tif","/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Peral_river/ZJ.shp","/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Peral_river/ZJ_fdir.tif")
    # sink.Cal_acc("/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Peral_river/ZJ_fdir.tif",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Peral_river/ZJ_acc.tif")

    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_fdir.tif","/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Adelaide_river/Adf.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Adelaide_river/Adelaide_fdir.tif")
    # sink.Cal_acc("/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Adelaide_river/Adelaide_fdir.tif",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Adelaide_river/Adelaide_acc.tif")

    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_fdir.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/TBP_river/TBP.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/TBP_river/TBP_fdir.tif")
    # sink.Cal_acc("/datanode05/zhangbin/TBP_Stream/DATA/valid/river/TBP_river/TBP_fdir.tif",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/TBP_river/TBP_acc.tif")

    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_fdir.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/NDvina/NDvina.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/NDvina/NDvina_fdir.tif")
    # sink.Cal_acc("/datanode05/zhangbin/TBP_Stream/DATA/valid/river/NDvina/NDvina_fdir.tif",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/NDvina/NDvina_acc.tif")

    def get_stream_order_2shp():
        import whitebox
        wbt = whitebox.WhiteboxTools()
        wbt.exe_path = '/datanode05/zhangbin/.conda/envs/zhangbin/lib/python3.8/site-packages/whitebox'  # Linux必加，否则无法运行
        # 5. 设定阈值提取河网（例如流量 > 100）
        wbt.extract_streams("/datanode05/zhangbin/TBP_Stream/DATA/valid/river/NDvina/NDvina_acc.tif",
                            "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/NDvina/NDvina_stream.tif", threshold=5555)

        # 6. 计算河网分级（斯特拉勒）
        # wbt.strahler_stream_order("/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Adelaide_river/Adelaide_fdir.tif",
        #                           "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Adelaide_river/Adelaide_stream.tif",
        #                           "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Adelaide_river/Adelaide_order_stream.tif")


        # 7. 栅格转矢量（输出为 Shapefile）
        wbt.raster_streams_to_vector("/datanode05/zhangbin/TBP_Stream/DATA/valid/river/NDvina/NDvina_stream.tif",
                                     "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/NDvina/NDvina_fdir.tif",
                                     "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/NDvina/NDvina_stream.shp")

        print("河网提取完成，矢量已保存为 river_network.shp")

    # get_stream_order_2shp()
    # application_endorheic.get_river_mouth_point("/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Peral_river/ZJ_fdir.tif",
    #                       "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Peral_river/ZJ_stream.tif",
    #                       "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Peral_river/ZJ_mouth.shp",
    #                       "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Peral_river/ZJ_stream.tif")
    # application_endorheic.get_river_mouth_point("/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Adelaide_river/Adelaide_fdir.tif",
    #                                             "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Adelaide_river/Adelaide_stream.tif",
    #                                             "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Adelaide_river/Adelaide_mouth.shp",
    #                                             "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/Adelaide_river/Adelaide_mouth.tif")
    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/valid/river/TBP_river/TBP_fdir.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/TBP_river/TBP_endorheic_point.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/TBP_river/TBP_endorheic_point.tif")
    # application_endorheic.get_river_mouth_point("/datanode05/zhangbin/TBP_Stream/DATA/valid/river/NDvina/NDvina_fdir.tif",
    #                                             "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/NDvina/NDvina_stream.tif",
    #                                             "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/NDvina/NDvina_mouth.shp",
    #                                             "/datanode05/zhangbin/TBP_Stream/DATA/valid/river/NDvina/NDvina_mouth.tif")

    # Merge.merge_sinks("/datanode05/zhangbin/TBP_Stream/DATA/Australia")
    # Merge.merge_sinks("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")
    # Merge.merge_sinks("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica")
    # Merge.merge_sinks("/datanode05/zhangbin/TBP_Stream/DATA/Asia")
    # Merge.merge_sinks("/datanode05/zhangbin/TBP_Stream/DATA/Europe")
    # Merge.merge_sinks("/datanode05/zhangbin/TBP_Stream/DATA/Siberia")
    # Merge.merge_sinks("/datanode05/zhangbin/TBP_Stream/DATA/Africa")


    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_depressions.tif",
    #                   "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_depressions.png")





    # 重新合并每个大洲，每个文件夹做0.05度缓冲区进行合并
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Asia/","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/New_As/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Asia")
    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_fdir1.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_endor1.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_endor1.tif")
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190/","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/EU/new_1190/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190/")
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020024230/","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/EU/new_4230/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020024230/")
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Europe/","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/EU/new_level2/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Europe")
    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_fdir1.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_endor1.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_endor1.tif")

    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/NA/new_NA")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")
    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_fdir1.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_endor1.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_endor1.tif")

    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110/Africa_1041259940","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/Af/new_9940/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110/Africa_1041259940")
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/Af/new_8110/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020018110")

    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020011530/Africa_1030011660","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/Af/new_1660/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020011530/Africa_1030011660")
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020011530","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/Af/new_1530/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020011530")

    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170/Africa_1040843040","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/Af/new_3040/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170/Africa_1040843040")

    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020000010","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/Af/new_0010/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020000010")

    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/Af/new_4170/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020034170")

    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/Af/new_7430/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_1020027430")

    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Africa","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/Af/new_Africa/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Africa")
    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_fdir2.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_endor2.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_endor2.tif")
    # application_endorheic.get_endorheic_from_final_point("/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_fdir1.tif",
    #                                                      "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_endor1.tif",
    #                                                      "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_endor_area.tif")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_endor_area.tif","/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_endor_area.png")
    # endorheic_process.Check1("/datanode05/zhangbin/TBP_Stream/DATA/Africa")
    # endorheic_process.sbatch_buchong_endorheic("/datanode05/zhangbin/TBP_Stream/DATA/Africa/buchong.1") #>15:11444   0-15:
    # endorheic_process.process_endorheic("/datanode05/zhangbin/TBP_Stream/DATA/Africa/buchong.1/41613")
    # Merge.merge_fill("/datanode05/zhangbin/TBP_Stream/DATA/Africa/buchong.1")


    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/new_SA/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica")
    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_fdir1.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_endor1.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_endor1.tif")


    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_5020049720","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/AU/new_9720/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_5020049720")
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Australia","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/AU/new_Australia/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Australia")
    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_fdir1.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_endor1.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_endor1.tif")

    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/AR/new_8900/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020008900_1")
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_8020000010_1","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/AR/new_0010/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Arctic")
    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_fdir1.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_endor1.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Arctic/Arctic_endor1.tif")

    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/SI/new_9230/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020009320")
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/SI/new_8670/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020008670")
    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020003790/","/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/SI/new_3790/")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_3020003790")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia")
    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_fdir1.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_endor1.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_endor1.tif")




    #  最后补充
    # endorheic_process.Check1("/datanode05/zhangbin/TBP_Stream/DATA/Afica")
    # endorheic_process.sbatch_buchong_endorheic("/datanode05/zhangbin/TBP_Stream/DATA/Africa/buchong.1") #>15:11444   0-15:
    # endorheic_process.process_endorheic("/datanode05/zhangbin/TBP_Stream/DATA/Africa/buchong.1/41613")
    # Merge.merge_fill("/datanode05/zhangbin/TBP_Stream/DATA/Africa/buchong.1")

    # endorheic_process.Check1("/datanode05/zhangbin/TBP_Stream/DATA/Asia")
    # endorheic_process.Check1("/datanode05/zhangbin/TBP_Stream/DATA/Europe")
    # endorheic_process.Check1("/datanode05/zhangbin/TBP_Stream/DATA/Australia")
    # endorheic_process.Check1("/datanode05/zhangbin/TBP_Stream/DATA/Arctic")
    # # endorheic_process.Check1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia")
    # endorheic_process.Check1("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")
    # endorheic_process.Check1("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica")
    # endorheic_process.Check1("/datanode05/zhangbin/TBP_Stream/DATA/Greenland")


    def change_permissions(root_dir, file_mode=0o777, dir_mode=0o777):
        for root, dirs, files in os.walk(root_dir):
            # 修改目录权限
            for d in dirs:
                dir_path = os.path.join(root, d)
                os.chmod(dir_path, dir_mode)
                print(f"Changed directory: {dir_path} to {oct(dir_mode)}")

            # 修改文件权限
            for f in files:
                file_path = os.path.join(root, f)
                os.chmod(file_path, file_mode)
                print(f"Changed file: {file_path} to {oct(file_mode)}")


    # 示例：更改 /path/to/your/directory 下所有文件为 644，目录为 755
    # change_permissions("/datanode05/zhangbin/TBP_Stream/DATA/Asia/buchong.1")
    # change_permissions("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/buchong.1")
    # change_permissions("/datanode05/zhangbin/TBP_Stream/DATA/Greenland/buchong.1")
    # change_permissions("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/buchong.1")
    # change_permissions("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/buchong.1")
    # change_permissions("/datanode05/zhangbin/TBP_Stream/DATA/Australia/buchong.1")
    # change_permissions("/datanode05/zhangbin/TBP_Stream/DATA/Europe/buchong.1")

    # endorheic_process.sbatch_buchong_endorheic("/datanode05/zhangbin/TBP_Stream/DATA/Asia/buchong.1")
    # endorheic_process.sbatch_buchong_endorheic("/datanode05/zhangbin/TBP_Stream/DATA/Greenland/buchong.1")
    # endorheic_process.sbatch_buchong_endorheic("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/buchong.1")
    # endorheic_process.sbatch_buchong_endorheic("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/buchong.1")
    # endorheic_process.sbatch_buchong_endorheic("/datanode05/zhangbin/TBP_Stream/DATA/Arctic/buchong.1")
    # endorheic_process.sbatch_buchong_endorheic("/datanode05/zhangbin/TBP_Stream/DATA/Australia/buchong.1")
    # endorheic_process.sbatch_buchong_endorheic("/datanode05/zhangbin/TBP_Stream/DATA/Europe/buchong.1")




    # Merge.merge_fill("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")
    # application_endorheic.get_endorheic_final_point(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_fdir2.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_endor2.shp",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_endor2.tif")
    # Merge.merge_fill("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica")
    # application_endorheic.get_endorheic_final_point(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_fdir2.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_endor2.shp",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/SouthAmerica_endor2.tif")
    # Merge.merge_fill("/datanode05/zhangbin/TBP_Stream/DATA/Australia")
    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_fdir2.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_endor2.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Australia/Australia_endor2.tif")

    # Merge.merge_fill("/datanode05/zhangbin/TBP_Stream/DATA/Asia")
    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_fdir2.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_endor2.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_endor2.tif")

    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_fdir2.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_endorheic2.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_endorheic2.tif")




    # Merge.merge_DEM("/datanode05/zhangbin/FABDEM_flat", [140, 74, 179, 55],
    #                 "/datanode05/zhangbin/TBP_Stream/DATA/regionDem/SI_140_74_179_55.tif")

    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/SI_140_74_179_55.tif","/datanode05/zhangbin/TBP_Stream/DATA/regionDem/SI_140_74_179_55.png")



    # Siberia_7570
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/SI_140_74_179_55.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/regionDem/Si_mask/3040017570.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040017570/Siberia_3040017570_burnDEM_.tif")

    # sink.Cal_sink("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040017570/Siberia_3040017570_burnDEM_.tif",
    #               "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040017570/Siberia_3040017570_burnsink_.tif")
    # sink.Cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040017570/Siberia_3040017570_burnDEM_.tif",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040017570/Siberia_3040017570_burndir_.tif")
    # calculate_sink_area_volume.calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040017570")
    #

    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia")
    # BB1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040017570/mask",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040017570/Siberia_3040017570_update_.csv")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia")

    # sbatch_clip_ModifiedDir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040017570",
    #                         "/datanode05/zhangbin/TBP_Stream/DATA/clip_Modifid_mask/EU/1190")
    # Merge.merge_continent1("/datanode05/zhangbin/TBP_Stream/DATA/Europe/Europe_2020071190")
    # application_endorheic.get_endorheic_final_point("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040017570/Siberia_3040017570_ModifiedDir_.tif",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040017570/Siberia_3040017570_endorheic_.shp",
    #                                                 "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040017570/Siberia_3040017570_endorheic_.tif")


    # 0910
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/SI_140_74_179_55.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/regionDem/Si_mask/3040010910.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010910/Siberia_3040010910_burnDEM_.tif")
    #
    # sink.Cal_sink("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010910/Siberia_3040010910_burnDEM_.tif",
    #               "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010910/Siberia_3040010910_burnsink_.tif")
    # sink.Cal_dir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010910/Siberia_3040010910_burnDEM_.tif",
    #              "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010910/Siberia_3040010910_burndir_.tif")
    # calculate_sink_area_volume.calculate_sink_volume_byDir("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010910")

    # # 0920
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/SI_140_74_179_55.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/regionDem/Si_mask/3040010920.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010920/Siberia_3040010920_burnDEM_.tif")
    #
    # sink.Cal_sink(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010920/Siberia_3040010920_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010920/Siberia_3040010920_burnsink_.tif")
    # sink.Cal_dir(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010920/Siberia_3040010920_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010920/Siberia_3040010920_burndir_.tif")
    # calculate_sink_area_volume.calculate_sink_volume_byDir(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010920")

    # 1770
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/regionDem/SI_140_74_179_55.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/regionDem/Si_mask/3040011770.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040011770/Siberia_3040011770_burnDEM_.tif")
    #
    # sink.Cal_sink(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040011770/Siberia_3040011770_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040011770/Siberia_3040011770_burnsink_.tif")
    # sink.Cal_dir(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040011770/Siberia_3040011770_burnDEM_.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040011770/Siberia_3040011770_burndir_.tif")
    # calculate_sink_area_volume.calculate_sink_volume_byDir(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040011770")


    # process_sink.sbatch_main1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia")
    # BB1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010920/mask",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010920/Siberia_3040010920_update_.csv")
    # BB1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010910/mask",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010910/Siberia_3040010910_update_.csv")
    # BB1("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040011770/mask",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040011770/Siberia_3040011770_update_.csv")
    # postprocess.sbatch_VisualCheck("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia")

    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040017570/Siberia_3040017570_ModifiedDir_.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/SI/new_9230/Siberia_3020009320_3040017570.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040017570/Siberia_3040017570_FABDIR1_.tif")
    #
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010920/Siberia_3040010920_ModifiedDir_.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/SI/new_9230/Siberia_3020009320_3040010920.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010920/Siberia_3040010920_FABDIR1_.tif")
    #
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010910/Siberia_3040010910_ModifiedDir_.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/SI/new_9230/Siberia_3020009320_3040010910.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010910/Siberia_3040010910_FABDIR1_.tif")
    #
    # Find.clip("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040011770/Siberia_3040011770_ModifiedDir_.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/buffer_clip_mask/SI/new_9230/Siberia_3020009320_3040011770.shp",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040011770/Siberia_3040011770_FABDIR1_.tif")
    #
    # result = ["/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_fdir1.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040017570/Siberia_3040017570_FABDIR1_.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010920/Siberia_3040010920_FABDIR1_.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040010910/Siberia_3040010910_FABDIR1_.tif",
    #           "/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia/Siberia_3040011770/Siberia_3040011770_FABDIR1_.tif"]
    # vrt_path = os.path.join("/datanode05/zhangbin/TBP_Stream/DATA/Siberia", "temp_merge.vrt")
    # vrt = gdal.BuildVRT(vrt_path, result)
    #
    # # 3. 将 VRT 转换为实际的 GeoTIFF 输出
    # gdal.Translate("/datanode05/zhangbin/TBP_Stream/DATA/Siberia/Siberia_fdir2.tif", vrt)
    # # 4. 关闭 VRT 对象
    # vrt = None
    # # （可选）删除中间 VRT 文件
    # os.remove(vrt_path)
    #

















    # """""""" 修复对向流 """"""""""

    Merge.merge_burnDEM("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica")
    # Merge.merge_burnDEM("/datanode05/zhangbin/TBP_Stream/DATA/Siberia")
    # Merge.merge_burnDEM("/datanode05/zhangbin/TBP_Stream/DATA/Asia")

    # 处理Africa对向流
    postprocess.divide_conuterflow_region(
        "/datanode05/zhangbin/TBP_Stream/DATA/Africa/counterflow/counterflowregions.tif",
        "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_BURN_DEM.tif",
        "/datanode05/zhangbin/TBP_Stream/DATA/Africa/Africa_fdir2.tif",
        "/datanode05/zhangbin/TBP_Stream/DATA/Africa/counterflow/counter_mask/")
    # postprocess.sbatch_PFD8_exorheic("/datanode05/zhangbin/TBP_Stream/DATA/Asia/counterflow/counter_mask")
    # postprocess.merge_counterflow_outdir("/datanode05/zhangbin/TBP_Stream/DATA/Asia/counterflow/counter_mask",
    #                                      "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_fdir2.tif")
    # postprocess.counterflow_process_check("/datanode05/zhangbin/TBP_Stream/DATA/Asia/Final_FDIR.tif")

    # 处理Asia对向流
    # postprocess.divide_conuterflow_region(
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Asia/counterflow/counterflowregions.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_BURN_DEM.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_fdir2.tif",
    #     "/datanode05/zhangbin/TBP_Stream/DATA/Asia/counterflow/counter_mask")
    # postprocess.sbatch_PFD8_exorheic("/datanode05/zhangbin/TBP_Stream/DATA/Asia/counterflow/counter_mask")
    # postprocess.merge_counterflow_outdir("/datanode05/zhangbin/TBP_Stream/DATA/Asia/counterflow/counter_mask",
    #                                      "/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_fdir2.tif")
    # postprocess.counterflow_process_check("/datanode05/zhangbin/TBP_Stream/DATA/Asia/Final_FDIR.tif")


    # 处理NA对向流

    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_BURN_DEM.tif","/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_BURN_DEM.png")
    postprocess.counterflow_process("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_fdir2.tif",
                                    "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/conuterflow/trail/counterflowregions.tif")
    postprocess.divide_conuterflow_region(
        "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/counterflow/counterflowregions.tif",
        "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_BURN_DEM.tif",
        "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_fdir2.tif",
        "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/counterflow/counter_mask")
    postprocess.sbatch_PFD8_exorheic("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/counterflow/counter_mask/")
    postprocess.merge_counterflow_outdir("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/counterflow/counter_mask/",
                                         "/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/NorthAmerica_fdir2.tif")
    # postprocess.counterflow_process_check("/datanode05/zhangbin/TBP_Stream/DATA/NorthAmerica/Final_FDIR.tif")


    split_cal_acc.cal_acc_work_flow_SA()


    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/FINAL_ACC.tif","/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/FINAL_ACC.png")
    # read_and_save_tif("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/Final_FDIR.tif","/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/Final_FDIR.png")
    # application_endorheic.get_SA()

    # a = "/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/region/"
    # for file in os.listdir(a):
    #     fileName = file.split('.')
    #     if not os.path.exists(os.path.join("/datanode05/zhangbin/TBP_Stream/DATA/SouthAmerica/cal_ACC/region_tree/",fileName[0]+"tree_.csv")):
    #         print(fileName)