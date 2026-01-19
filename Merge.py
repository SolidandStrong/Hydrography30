# -*- coding: utf-8 -*-
"""
@Time ： 2024/12/23 17:00
@Auth ：
@File ：Merge.py
@IDE ：PyCharm
"""
import os
import shutil

import Raster
from db import *
from process_sink import *
import whitebox
import numpy as np
wbt = whitebox.WhiteboxTools()

def merge_DEM(FABVenu,extent,outfile):
    """
    根据四至来筛选分块DEM进行合并
    :param FABVenu:
    :param extent:  [left_lon,top_lat,right_lon,bottom_lat] int/float
    :param outfile:
    :return:
    """
    files = []

    left_lon = int(extent[0])
    # if left_lon < extent[0] :
    #     left_lon -= 1
    top_lat = int(extent[1])
    right_lon = int(extent[2])
    bottom_lat = int(extent[3])

    for i in range(left_lon,right_lon+1):
        for j in range(bottom_lat,top_lat+1):
            if i >= 0 :
                flag2 = 'E'
            else:
                flag2 = 'W'
            if j >= 0 :
                flag1 = 'N'
            else:
                flag1 = 'S'

            fileName = flag1 + "{:02d}".format(abs(j)) + flag2 + "{:03d}".format(abs(i)) + '_FABDEM_V1-0.tif'

            filePath = os.path.join(FABVenu,fileName)
            if os.path.exists(filePath):
                files.append(filePath)
                print("****",filePath)
            else:
                print(fileName)
    # 构造 BuildVRTOptions（正确方式）
    vrt_options = gdal.BuildVRTOptions(
        separate=False,
        resolution='highest',
        VRTNodata=-9999  # 设定 NoData 值，可根据你的影像情况修改
    )

    # 构建 VRT

    vrt_path = os.path.join(os.path.dirname(outfile), 'mergedem.vrt')
    vrt_ds = gdal.BuildVRT(vrt_path, files, options=vrt_options)
    vrt_ds = None  # 保存并关闭

    # 输出 GeoTIFF + 压缩
    gdal.Translate(
        outfile,
        vrt_path,
        format="GTiff",
        creationOptions=[
            "COMPRESS=LZW",
            "TILED=YES",
            "BIGTIFF=IF_SAFER"
        ]
    )
    os.remove(vrt_path)
    # gdal.Warp(outfile, files, options=[
    #     "-co", "COMPRESS=LZW",  # 选择LZW压缩
    #     "-co", "BIGTIFF=YES"])  # 允许大文件])
    print("合并完成，输出文件：", outfile)

def merge_endorheic(Venu):
    """
    合并流域内内流区的流向:
    根据数据库中的起始坐标合并
    :param Venu:
    :return:
    """
    baseName = os.path.basename(Venu)

    DIRFILE = os.path.join(Venu,baseName + "_outDIR_.tif")
    fDIR = Raster.get_raster(DIRFILE)
    proj,geo,d_nodata = Raster.get_proj_geo_nodata(DIRFILE)

    fDIR[fDIR == d_nodata] = 255
    saveRaster = os.path.join(Venu,baseName + "_ModifiedDir.tif")

    db_name = os.path.join(Venu,baseName + "_burndb.db")
    print('Search the endorheic from the database')
    infos = query_data_byType(db_name, str(2))

    endorheicVenu = os.path.join(Venu,"endorheic") # 进入endorheic目录下
    for info in infos:
        sinkId = info[0]
        minX = info[-5]
        maxX = info[-4]
        minY = info[-3]
        maxY = info[-2]

        # 掩膜sink
        # maskDem = dem[minX:maxX + 1, minY:maxY + 1].copy()
        # maskSink = sink[minX:maxX + 1, minY:maxY + 1].copy()

        sinkVenu = os.path.join(endorheicVenu,str(sinkId))   # 进入待合并的sink目录
        sinkoutdir = os.path.join(sinkVenu,str(sinkId)+"_outdir_.tif")

        outdir = Raster.get_raster(sinkoutdir)
        _,_,sinknodata = Raster.get_proj_geo_nodata(sinkoutdir)

        maskrow, maskcol = outdir.shape
        # 写入流向
        for ti in range(maskrow):
            for tj in range(maskcol):
                if outdir[ti, tj] in [sinknodata]:
                    continue
                fDIR[minX + ti, minY + tj] = outdir[ti, tj]

    Raster.save_raster(saveRaster,fDIR,proj,geo,gdal.GDT_Byte,255)

def merge_continent1(venu):
    """
    将大洲内的大流域合并
    :param venu:
    :return:
    """
    baseName = os.path.basename(venu)
    dirNames = os.listdir(venu)

    inputfiles = []
    for dirName in dirNames:
        if len(dirName.split('.')) > 1:
            continue
        dirPath = os.path.join(venu,dirName)
        outDir = os.path.join(dirPath,dirName+'_FABDIR1_.tif')
        inputfiles.append(outDir)
        print(outDir)

    outfile = os.path.join(venu,baseName + '_fdir1.tif')
    gdal.Warp(outfile,inputfiles,options=[
        "-co", "COMPRESS=LZW",  # 选择LZW压缩
        "-co", "BIGTIFF=YES"])    # 允许大文件])
    print("合并完成，输出文件：", outfile)

def clip_modifiedDir(input_vector,input_raster,venu):
    from osgeo import gdal

    baseName = os.path.basename(venu)
    output_raster = os.path.join(venu,baseName+'_FABDIR1_.tif')
    if os.path.exists(output_raster):
        return
    # 打开输入栅格
    src_ds = gdal.Open(input_raster)

    # 执行裁剪操作
    gdal.Warp(
        output_raster,
        src_ds,
        format="GTiff",
        cutlineDSName=input_vector,
        cropToCutline=False,  # 确保边界经过的栅格不被裁剪掉
        dstNodata=0,  # 指定 NoData 值
        cutlineBlend=0,  # 设置为 0 确保边界保持精确
        creationOptions=["COMPRESS=LZW"]  # 可选的压缩选项
    )

    # 关闭数据集
    src_ds = None
    print("裁剪完成！")

def merge_sinks(venu):
    """
    遍历每个大洲的文件夹，讲含有endorheic文件夹的sink文件进行合并，得到大洲的sink.tif
    :param venu:  ./NorthAmerica
    :return:
    """
    baseName = os.path.basename(venu)
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
        temp_dir = os.path.join(popPath, 'endorheic')
        if os.path.exists(temp_dir):
            # 存储

            tempfils = os.listdir(temp_dir)
            for i in tempfils:
                # print(os.path.join(temp_dir, i, i + '_mask_.tif'))
                if os.path.join(temp_dir, i, i + '_mask_.tif') in result:
                    continue

                result.append(os.path.join(temp_dir, i, i + '_mask_.tif'))

        for temp_dirName in os.listdir(popPath):
            if len(temp_dirName.split('.')) > 1:
                continue
            if len(temp_dirName.split('_')) == 0:
                continue
            popPaths.append(os.path.join(popPath,temp_dirName))
    dataset = []
    for temp_path in result:
        temp_source = gdal.Open(temp_path)
        if temp_source == None:
            continue
        print(temp_source)
        dataset.append(temp_source)
    outfile = os.path.join(venu, baseName + '_depressions.tif')
    gdal.Warp(outfile, dataset, options=[
        "-co", "COMPRESS=LZW",  # 选择LZW压缩
        "-co", "BIGTIFF=YES"])  # 允许大文件])
    print("合并完成，输出文件：", outfile)

    # print(result)
    # with open("/datanode05/zhangbin/TBP_Stream/DATA/CHECK.csv",'w',newline='') as f:
    #     writer = csv.writer(f)
    #     writer.writerows(result)
    #     f.close()
    # os.chmod("/datanode05/zhangbin/TBP_Stream/DATA/CHECK.csv",0o777)

def merge_fill(Venu):
    """
    遍历每个大洲的文件夹，讲含有endorheic文件夹的sink文件进行合并，得到大洲的sink.tif
    :param venu:  ./NorthAmerica
    :return:
    """
    baseName = os.path.basename(Venu)
    result = [os.path.join(Venu,baseName+"_fdir1.tif")]
    venu = os.path.join(Venu,'buchong.1')
    dirNames = os.listdir(venu)
    popPaths = [os.path.join(venu, dirName) for dirName in dirNames]
    for popPath in popPaths:
        filenames = os.listdir(popPath)
        for filename in filenames:
            result.append(os.path.join(popPath,filename,filename+'_outdir_.tif'))

    outfile = os.path.join(Venu,baseName+"_fdir2.tif")
    # gdal.Warp(outfile, result, options=[
    #     "-co", "COMPRESS=LZW",  # 选择LZW压缩
    #     "-co", "BIGTIFF=YES"])  # 允许大文件])
    # print("合并完成，输出文件：", outfile)

    # 2. 构建虚拟栅格（VRT）
    vrt_path = os.path.join(Venu,"temp_merge.vrt")
    vrt = gdal.BuildVRT(vrt_path, result)

    # 3. 将 VRT 转换为实际的 GeoTIFF 输出
    gdal.Translate(outfile, vrt)
    # 4. 关闭 VRT 对象
    vrt = None

    # （可选）删除中间 VRT 文件
    os.remove(vrt_path)

def merge_burnDEM(venu):
    """
    合并该目录下所有的burn_DEM

    :param venu:
    :return:
    """
    baseName = os.path.basename(venu)
    result = []
    dirNames = os.listdir(venu)
    popPaths = [os.path.join(venu, dirName) for dirName in dirNames]
    while popPaths:
        popPath = popPaths.pop()

        dirName = os.path.basename(popPath)
        if len(dirName.split('.')) > 1:
            continue
        if len(dirName.split('_')) == 0:
            continue
        if not os.path.exists(os.path.join(popPath, dirName + '_ClipdDir.tif')):
            continue
        if os.path.exists(os.path.join(popPath, dirName + '_burnDEM.tif')):
            result.append(os.path.join(popPath, dirName + '_burnDEM.tif'))

        # for temp_dirName in os.listdir(popPath):
        #     if len(temp_dirName.split('.')) > 1:
        #         continue
        #     if len(temp_dirName.split('_')) == 0:
        #         continue
        #     popPaths.append(os.path.join(popPath, temp_dirName))
    # print(result)
    # 构造 BuildVRTOptions（正确方式）
    vrt_options = gdal.BuildVRTOptions(
        separate=False,
        resolution='highest',
        VRTNodata=0  # 设定 NoData 值，可根据你的影像情况修改
    )

    # 构建 VRT
    vrt_path = os.path.join(venu, "BurnDEM1.vrt")
    vrt_ds = gdal.BuildVRT(vrt_path, result, options=vrt_options)
    vrt_ds = None  # 保存并关闭

    # 输出 GeoTIFF + 压缩
    gdal.Translate(
        os.path.join(venu, "BURN_DEM.tif"),
        vrt_path,
        format="GTiff",
        creationOptions=[
            "COMPRESS=LZW",
            "TILED=YES",
            "BIGTIFF=IF_SAFER"
        ]
    )

    ref_raster = os.path.join(venu, "initial_dir.tif")  # 参考栅格（裁剪范围与行列号都与它一致）
    src_raster = os.path.join(venu, "BURN_DEM.tif")
    out_raster = os.path.join(venu, baseName+"_BURN_DEM.tif")

    fdir = Raster.get_raster(ref_raster)
    proj, geo, f_nodata = Raster.get_proj_geo_nodata(ref_raster)
    print(geo)

    dem = Raster.get_raster(src_raster)
    proj_1, geo_1, d_nodata = Raster.get_proj_geo_nodata(src_raster)
    print(geo_1)

    row, col = fdir.shape
    print(row, col)
    start_col = round((geo[0] - geo_1[0]) / geo_1[1])
    start_row = round((geo[3] - geo_1[3]) / geo_1[5])
    print(start_row, start_col)

    result = dem[start_row:start_row + row, start_col:start_col + col]
    Raster.save_raster(out_raster, result, proj, geo, gdal.GDT_Float32, d_nodata)

    print("✅ 栅格掩膜完成，结果保存至：", out_raster)


def merge_update_fdr(new_acc_dir,out_final_acc):
    """
    合并最终更新的acc。
    :param new_acc_dir:
    :param out_final_acc:
    :return:
    """
    files = [os.path.join(new_acc_dir,fileName) for fileName in os.listdir(new_acc_dir)]
    proj,geo,nodata = Raster.get_proj_geo_nodata(files[0])
    # 构造 BuildVRTOptions（正确方式）
    vrt_options = gdal.BuildVRTOptions(
        separate=False,
        resolution='highest',
        VRTNodata=nodata  # 设定 NoData 值，可根据你的影像情况修改
    )
    baseDir = os.path.dirname(out_final_acc)
    # 构建 VRT
    vrt_path = os.path.join(baseDir, "temp_merged3.vrt")
    vrt_ds = gdal.BuildVRT(vrt_path, files, options=vrt_options)
    vrt_ds = None  # 保存并关闭


    # 输出 GeoTIFF + 压缩
    gdal.Translate(
        out_final_acc,
        vrt_path,
        format="GTiff",
        creationOptions=[
            "COMPRESS=LZW",
            "TILED=YES",
            "BIGTIFF=IF_SAFER"
        ]
    )
    os.chmod(vrt_path,0o777)
    os.remove(vrt_path)

def mask_fdr_acc(ref_raster,src_raster,out_raster):

    # ref_raster = os.path.join(venu, "initial_dir.tif")  # 参考栅格（裁剪范围与行列号都与它一致）
    # src_raster = os.path.join(venu, "BURN_DEM.tif")
    # out_raster = os.path.join(venu, baseName + "_BURN_DEM.tif")

    fdir = Raster.get_raster(ref_raster)
    proj, geo, f_nodata = Raster.get_proj_geo_nodata(ref_raster)
    print(geo)
    print(fdir.shape)
    dem = Raster.get_raster(src_raster)
    proj_1, geo_1, d_nodata = Raster.get_proj_geo_nodata(src_raster)
    print(geo_1)

    print(dem.shape)
    # row, col = fdir.shape
    # print(row, col)
    # start_col = round((geo[0] - geo_1[0]) / geo_1[1])
    # start_row = round((geo[3] - geo_1[3]) / geo_1[5])
    # print(start_row, start_col)
    #
    # result = dem[start_row:start_row + row, start_col:start_col + col]
    # Raster.save_raster(out_raster, result, proj, geo, gdal.GDT_Float32, d_nodata)
    #
    # print("✅ 栅格掩膜完成，结果保存至：", out_raster)
if __name__ == "__main__":

    # merge_endorheic("/datanode05/zhangbin/TBP_Stream/DATA/Asia/Asia_4020050470/")
    # 参考栅格路径（保持行列号一致）
    # ref_raster = os.path.join(venu, baseName + "_fdir2.tif")  # 参考栅格（裁剪范围与行列号都与它一致）
    # # 原始需要裁剪的栅格
    # src_raster = os.path.join(venu, "BURN_DEM.tif")  # 被裁剪的原始影像
    # # 输出文件路径
    # out_raster = os.path.join(venu, baseName + "_BURN_DEM.tif")
    ref_raster = r'F:\青藏高原水体数据集\description\制图\内流区\对流\10\dem_10_.tif'  # 参考栅格（裁剪范围与行列号都与它一致）
    src_raster = r'F:\青藏高原水体数据集\description\制图\内流区\各大洲内流终点\mask_fdir.tif'
    out_raster = r'F:\青藏高原水体数据集\description\制图\内流区\各大洲内流终点\MSK.tif'

    fdir = Raster.get_raster(ref_raster)
    proj,geo,f_nodata = Raster.get_proj_geo_nodata(ref_raster)
    print(geo)

    dem = Raster.get_raster(src_raster)
    proj_1,geo_1,d_nodata = Raster.get_proj_geo_nodata(src_raster)
    print(geo_1)

    row,col = fdir.shape
    print(row,col)
    start_col = round((geo[0]-geo_1[0])/geo_1[1])
    start_row = round((geo[3]-geo_1[3])/geo_1[5])
    print(start_row,start_col)

    result = dem[start_row:start_row+row,start_col:start_col+col]
    Raster.save_raster(out_raster,result,proj,geo,gdal.GDT_Float32,d_nodata)