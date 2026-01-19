# -*- coding: utf-8 -*-
"""
@Time ： 2025/6/22 15:28
@Auth ：
@File ：process_FAB_SRTM.py
@IDE ：PyCharm

该文件下函数用于处理FABDEM中缺失的部分，使用SRTM补全并平滑
"""
import os
import shutil

import Merge
import Raster
from osgeo import gdal
import numpy as np
import matplotlib.pyplot as plt

import valid
# plt.rcParams["font.family"] = "Times New Roman"  # 设置全局字体
from pylab import mpl
import special_process

mpl.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams["font.family"] = "Times New Roman"  # 设置全局字体


def modify_SRTM_notda(SRTM_dir):
    """
    将SRTM的nodata更新为-9999，与FABDEM保持统一
    :param SRTM_dir:
    :return:
    """
    for file in os.listdir(SRTM_dir):

        if file.split('.')[-1] != 'tif':
            continue
        file_path = os.path.join(SRTM_dir,file)

        DEM = Raster.get_raster(file_path)
        proj,geo,d_nodata = Raster.get_proj_geo_nodata(file_path)

        DEM[DEM == d_nodata] = -9999

        Raster.save_raster(file_path,DEM,proj,geo,gdal.GDT_Float32,-9999)
        print("{:s} has been successfully processed.".format(file))

def merge_SRTM_FABDEM(SRTM_dir,FAB_dir,outFABDEM,outSRTM):
    """
    合并SRTM部分
    :param SRTM_dir:
    :param FAB_dir:
    :param outFABDEM:
    :param outSRTM:
    :return:
    """

    SRTMs = []
    FABs = []
    for file in os.listdir(SRTM_dir):
        if file.split('.')[-1] != 'tif':
            continue

        # infos = file.split('_')
        SRTMs.append(os.path.join(SRTM_dir,file))

    # print(result)
    # 构造 BuildVRTOptions（正确方式）
    vrt_options = gdal.BuildVRTOptions(
        separate=False,
        resolution='highest',
        VRTNodata=-9999  # 设定 NoData 值，可根据你的影像情况修改
    )

    # 构建 VRT
    vrt_path = os.path.join("SRTM.vrt")
    vrt_ds = gdal.BuildVRT(vrt_path, SRTMs, options=vrt_options)
    vrt_ds = None  # 保存并关闭

    # 输出 GeoTIFF + 压缩
    gdal.Translate(
        os.path.join(outSRTM),
        vrt_path,
        format="GTiff",
        creationOptions=[
            "COMPRESS=LZW",
            "TILED=YES",
            "BIGTIFF=IF_SAFER"
        ]
    )
    os.remove(vrt_path)

    # for file in os.listdir(FAB_dir):
    #     if file.split('.')[-1] != 'tif':
    #         continue
    #
    #     # infos = file.split('_')
    #     FABs.append(os.path.join(FAB_dir,file))
    #
    # # print(result)
    # # 构造 BuildVRTOptions（正确方式）
    # vrt_options = gdal.BuildVRTOptions(
    #     separate=False,
    #     resolution='highest',
    #     VRTNodata=-9999  # 设定 NoData 值，可根据你的影像情况修改
    # )
    #
    # # 构建 VRT
    # vrt_path = os.path.join("FABs.vrt")
    # vrt_ds = gdal.BuildVRT(vrt_path, FABs, options=vrt_options)
    # vrt_ds = None  # 保存并关闭
    #
    # # 输出 GeoTIFF + 压缩
    # gdal.Translate(
    #     os.path.join(outFABDEM),
    #     vrt_path,
    #     format="GTiff",
    #     creationOptions=[
    #         "COMPRESS=LZW",
    #         "TILED=YES",
    #         "BIGTIFF=IF_SAFER"
    #     ]
    # )
    # os.remove(vrt_path)

def clip_FABDEM(SRTM_file,FAB_file,out_path):
    from osgeo import gdal

    def align_raster_to_reference(src_path, ref_path, out_path):
        # 读取参考栅格信息
        ref_ds = gdal.Open(ref_path)
        gt = ref_ds.GetGeoTransform()
        proj = ref_ds.GetProjection()
        xres = gt[1]
        yres = -gt[5]
        xmin = gt[0]
        ymax = gt[3]
        xmax = xmin + ref_ds.RasterXSize * xres
        ymin = ymax - ref_ds.RasterYSize * yres

        # 使用 Warp 对齐
        gdal.Warp(
            out_path,
            src_path,
            format='GTiff',
            outputBounds=[xmin, ymin, xmax, ymax],
            xRes=xres,
            yRes=yres,
            dstSRS=proj,
            resampleAlg='bilinear',  # 可选：'near', 'bilinear', 'cubic'
            targetAlignedPixels=True,
            dstNodata=-9999
        )

        print(f"✅ 已将 {src_path} 对齐至 {ref_path}，输出为 {out_path}")

    align_raster_to_reference(FAB_file, SRTM_file, out_path)
    A=Raster.get_raster(out_path)
    print(A.shape)
    B = Raster.get_raster(SRTM_file)
    print(B.shape)

def cancha_model_smooth_connect(fabdem_path,srtm_path,out_path):
    import numpy as np
    from osgeo import gdal
    from scipy.ndimage import gaussian_filter

    # ===== 工具函数 =====
    def read_raster(path):
        ds = gdal.Open(path)
        arr = ds.ReadAsArray().astype(np.float32)
        nodata = ds.GetRasterBand(1).GetNoDataValue()
        gt = ds.GetGeoTransform()
        proj = ds.GetProjection()
        return arr, nodata, gt, proj

    def write_raster(path, arr, gt, proj, nodata):
        driver = gdal.GetDriverByName("GTiff")
        ds = driver.Create(path, arr.shape[1], arr.shape[0], 1, gdal.GDT_Float32)
        ds.SetGeoTransform(gt)
        ds.SetProjection(proj)
        band = ds.GetRasterBand(1)
        band.WriteArray(arr)
        band.SetNoDataValue(nodata)
        ds.FlushCache()

    # ===== 读取数据 =====


    fab, fab_nodata, gt, proj = read_raster(fabdem_path)
    srtm, srtm_nodata, _, _ = read_raster(srtm_path)
    fab = fab[:-1,:-1].copy()

    # ===== 有效区域 =====
    fab_valid = (fab != fab_nodata) & (~np.isnan(fab))
    srtm_valid = (srtm != srtm_nodata) & (~np.isnan(srtm))
    overlap = fab_valid & srtm_valid

    # ===== 计算残差 =====
    residual = np.full_like(fab, np.nan)
    residual[overlap] = fab[overlap] - srtm[overlap]

    # ===== 平滑残差（高斯滤波） =====
    smoothed_residual = gaussian_filter(np.nan_to_num(residual), sigma=5)

    # 掩掉原始没有值nodata的位置（不可信）
    smoothed_residual[~srtm_valid] = 0

    # ===== 用 smoothed residual 修正 SRTM =====
    srtm_adjusted = srtm + smoothed_residual

    # ===== 构建融合 DEM =====
    fused = np.full_like(fab, fab_nodata)
    fused[fab_valid] = fab[fab_valid]  # FABDEM 优先
    fused[~fab_valid & srtm_valid] = srtm_adjusted[~fab_valid & srtm_valid]  # 用调整后的 SRTM 填补

    # ===== 输出结果 =====
    write_raster(out_path, fused, gt, proj, fab_nodata)
    print(f"✅ 残差建模融合完成：{out_path}")


# ---------------------------- 残差评估代码 -----------------------------------
def PDF_CDF(modified_dem,fabdem):

    modified_dem_arr = Raster.get_raster(modified_dem)
    proj,geo,m_nodata = Raster.get_proj_geo_nodata(modified_dem)
    fabdem_arr = Raster.get_raster(fabdem)
    proj,geo,f_nodata = Raster.get_proj_geo_nodata(fabdem)
    fabdem_arr = fabdem_arr[:-1,:-1]
    # mask = (modified_dem_arr != m_nodata) & (fabdem_arr != f_nodata)
    #
    #
    modified_dem_arr = modified_dem_arr[modified_dem_arr!=m_nodata]
    # fabdem_arr = fabdem_arr[mask]
    fabdem_arr = fabdem_arr[fabdem_arr!=f_nodata]
    # diff = dem1 - dem2   # <-- 换成你实际的高程差数组

    diff = fabdem_arr#modified_dem_arr
    diff = diff.flatten()  # 展平成一维
    diff = diff[np.isfinite(diff)]  # 去除 NaN 和 inf

    diff1 = modified_dem_arr
    diff1 = diff1.flatten()  # 展平成一维
    diff1 = diff1[np.isfinite(diff1)]  # 去除 NaN 和 inf

    # print(np.max(diff),np.min(diff),np.mean(diff))

    fig, ax = plt.subplots()
    ax.hist(diff, bins=60, density=True,alpha=0.5, color="steelblue",label = 'Original FABDEM')
    ax.hist(diff1, bins=60, density=True, alpha=0.5, color="darkorange",label = 'Modified FABDEM')
    # ax.set_title("Probability Density Function (PDF) of Elevation Difference")
    # ax.set_xlabel("Elevation difference (m)")
    # ax.set_ylabel("Density")
    plt.legend()
    plt.savefig("/datanode05/zhangbin/FAB_hydrography/initial_data/SRTM/two_PDF.svg", transparent=True)
    # plt.savefig("/home/zhangbin/TBP_Stream/CODE/PDF.svg")
    plt.close()
    # plt.show()

    fig, ax = plt.subplots()
    sorted_diff = np.sort(diff)
    cdf = np.arange(1, len(diff) + 1) / len(diff)
    ax.plot(sorted_diff, cdf,color="steelblue",label = 'Original FABDEM')

    sorted_diff1 = np.sort(diff1)
    cdf1 = np.arange(1, len(diff1) + 1) / len(diff1)
    ax.plot(sorted_diff1, cdf1, color="darkorange",label = 'Modified FABDEM')

    # ax.set_title("Cumulative Distribution Function (CDF) of Elevation Difference")
    # ax.set_xlabel("Elevation difference (m)")
    # ax.set_ylabel("CDF")
    plt.legend()
    plt.savefig("/datanode05/zhangbin/FAB_hydrography/initial_data/SRTM/two_CSF.svg", transparent=True)
    # plt.savefig("/home/zhangbin/TBP_Stream/CODE/CDF.svg")
    # plt.show()
    plt.close()


def merge_AW3D(mask_dir,out_fdr_file):
    """
    将修复后的对向流栅格上游流向合并至原始流向进行修正。
    :param mask_dir:
    :param outfdir_file:
    :return:
    """


    paras = ["/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_fdr2.tif","/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/修正/AW3D30/AR1_fdr_mask.tif"]

    # for dirPath in os.listdir(mask_dir):
    #     if dirPath.split('.')[-1] != 'tif':
    #         continue
    #
    #     dirfile = os.path.join(mask_dir,dirPath)
    #     if os.path.exists(dirfile):
    #         paras.append(dirfile)
    #     # print(dirfile)

    venu = os.path.dirname(out_fdr_file)
    # 构造 BuildVRTOptions（正确方式）
    vrt_options = gdal.BuildVRTOptions(
        separate=False,
        resolution='highest',
        VRTNodata= 255 #-9999  # 设定 NoData 值，可根据你的影像情况修改
    )

    # 构建 VRT
    vrt_path = os.path.join(venu, "temp_merged6.vrt")
    # if os.path.exists(vrt_path):
    #     os.remove(vrt_path)
    vrt_ds = gdal.BuildVRT(vrt_path, paras, options=vrt_options)
    vrt_ds = None  # 保存并关闭

    # 输出 GeoTIFF + 压缩
    # os.path.join(venu, "Final_FDIR.tif")
    gdal.Translate(
        out_fdr_file,
        vrt_path,
        format="GTiff",
        creationOptions=[
            "COMPRESS=LZW",
            "TILED=YES",
            "BIGTIFF=IF_SAFER"
        ]
    )
    os.remove(vrt_path)

if __name__ == "__main__":

    # print('1')


    # ------------------------------------ 处理里海区域的DEM -------------------------------------------------
    # modify_SRTM_notda("/data/zhangbin/TBP_Stream_Data/SRTM")
    # merge_SRTM_FABDEM("/data/zhangbin/TBP_Stream_Data/SRTM","/data/zhangbin/TBP_Stream_Data/FABDEM","/data/zhangbin/TBP_Stream_Data/FAB.tif","/data/zhangbin/TBP_Stream_Data/SRTM.tif")
    #
    # clip_FABDEM("/data/zhangbin/TBP_Stream_Data/SRTM.tif","/data/zhangbin/TBP_Stream_Data/FAB.tif","/data/zhangbin/TBP_Stream_Data/clip_FAB.tif")
    #
    # cancha_model_smooth_connect("/data/zhangbin/TBP_Stream_Data/clip_FAB.tif","/data/zhangbin/TBP_Stream_Data/SRTM.tif","/data/zhangbin/TBP_Stream_Data/modified_FABDEM.tif")





    # ------------------------------------------------ 处理五大湖处的DEM -----------------------------------------------------------
    # modify_SRTM_notda("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatLakeSRTM/")
    # merge_SRTM_FABDEM("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatLakeSRTM/", "/data/zhangbin/TBP_Stream_Data/FABDEM",
    #                   "/data/zhangbin/TBP_Stream_Data/FAB.tif", "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_srtm.tif")
    #
    # A = Raster.get_raster("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_srtm.tif")
    # proj,geo,nodata = Raster.get_proj_geo_nodata("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_srtm.tif")
    #
    #
    # row1,col1 = valid.lonlat_to_pixel(geo,-87.1,48.1)
    # row2, col2 = valid.lonlat_to_pixel(geo,-85.9,46.9 )
    #
    # A[row1:row2,col1:col2] = 176
    # Raster.save_raster("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_srtm1.tif",A,proj,geo,gdal.GDT_Float32,nodata)
    #
    # # #
    # clip_FABDEM("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_srtm1.tif",
    #             "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/final_fdr/Basin/five_lake/big/dem_big.tif","/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_FAB.tif")
    # #
    # cancha_model_smooth_connect("/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_FAB.tif",
    #                             "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_srtm1.tif",
    #                             "/datanode05/xiefy/zhangb/FAB_Hydro/run_data/greatlake_modified_FABDEM.tif")

    # ---------------------------------------- 格陵兰岛的AW3D30-DEM -------------------------------------------------
    # outvenu = r'G:\datanode05\FABDEM_stream\initial_data\AW3D_Greenland\Arctic\Arctic'
    # for file in os.listdir('G:\datanode05\FABDEM_stream\initial_data\AW3D_Greenland\Arctic\新建文件夹'):
    #     venu = os.path.join(r'G:\datanode05\FABDEM_stream\initial_data\AW3D_Greenland\Arctic\新建文件夹',file)
    #     for dem in os.listdir(venu):
    #         if dem.split('.')[-1] != 'tif':
    #             continue
    #         if dem.split('.')[0].split('_')[-1] != 'DSM':
    #             continue
    #         if os.path.exists(os.path.join(outvenu,dem)):
    #             continue
    #         shutil.copyfile(os.path.join(venu,dem),os.path.join(outvenu,dem))
    #         print(dem)
    #     print(file)

    # A = Raster.get_raster(r'F:\青藏高原水体数据集\New_DATA\全局环\AR\AR1_fdr.tif')
    # proj,geo,nodata = Raster.get_proj_geo_nodata(r'F:\青藏高原水体数据集\New_DATA\全局环\AR\AR1_fdr.tif')
    #
    # A[A!=nodata] = 6
    # Raster.save_raster(r'F:\青藏高原水体数据集\New_DATA\全局环\AR\AR1_fdr_mask.tif',A,proj,geo,gdal.GDT_Byte,nodata)
    merge_AW3D("","/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/修正/AW3D30/AR1_fdr_6.tif")
    A = Raster.get_raster("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/修正/AW3D30/AR1_fdr_6.tif")
    proj,geo,nodata = Raster.get_proj_geo_nodata("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/修正/AW3D30/AR1_fdr_6.tif")

    A[A == 6] = 255
    Raster.save_raster("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/Arctic_fdr3.tif",A,proj,geo,gdal.GDT_Byte,nodata)


    # merge_AW3D("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/修正/AW3D30/Arctic/",
    #            "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/修正/AW3D30/AW3D30_Arctic.tif")
    # special_process.read_and_save_tif("/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/修正/AW3D30/AW3D30_Arctic.tif",
    #                                   "/datanode05/tangjj/zhangb/FAB_Hydro/run_data/Fdr/Arctic/修正/AW3D30/AW3D30_Arctic.png")
    # merge_AW3D("/datanode05/zhangbin/FAB_hydrography/initial_data/AW3D/Greenland","/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/GreenLand_AW3D.tif")
    # special_process.read_and_save_tif("/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/GreenLand_AW3D.tif","/datanode05/zhangbin/FAB_hydrography/initial_data/Continent_dem/GreenLand_AW3D.png")


    # --------------------------------------------- 绘制残差图 -----------------------------------------------------------
    # PDF_CDF("/datanode05/zhangbin/FAB_hydrography/initial_data/SRTM/modified_FABDEM.tif","/datanode05/zhangbin/FAB_hydrography/initial_data/SRTM/clip_FAB.tif")
    # PDF_CDF("/data/zhangbin/TBP_Stream_Data/modified_FABDEM.tif","/data/zhangbin/TBP_Stream_Data/clip_FAB.tif")