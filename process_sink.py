# -*- coding: utf-8 -*-
"""
@Time ： 2024/9/21 21:14
@Auth ：
@File ：process_sink.py
@IDE ：PyCharm
"""
import csv
import os.path
import numpy as np
# from osgeo import gdal
import Raster
import genral_functions
import process_sink
import sink
from genral_functions import *
from db import *
from multiprocessing import Pool
from functools import partial

# from C_funs import *
dmove = [(0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1), (-1, 0), (-1, 1)]
dmove_dic = {1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1), 16: (0, -1), 32: (-1, -1), 64: (-1, 0), 128: (-1, 1)}
cellSize = 30  # Size of cell


# gdal.DontUseExceptions()

def getInfo(db_name, sinkType):
    return query_data_byType(db_name, str(sinkType))


def soD8(mask, masknoData, startcell, newDir):
    row, col = mask.shape
    weight = np.zeros((row, col))
    weight[startcell[0], startcell[1]] = 1

    fDiss = np.zeros((row, col))
    vis = np.zeros((row, col))

    while True:

        if np.all(vis == 1) or np.all(weight == 0):
            break
        # print(vis)
        cells = np.where(weight == weight.max())
        cell = (cells[0][0], cells[1][0])

        # if vis[cell[0],cell[1]] == 1:
        #     break
        # print(cell)
        # print(weight)
        vis[cell[0], cell[1]] = 1
        weight[cell[0], cell[1]] = 0
        for k in range(8):
            centerCell = (cell[0] + dmove[k][0], cell[1] + dmove[k][1])

            if not check_boundary(centerCell[0], centerCell[1], row, col):
                continue
            if mask[centerCell[0], centerCell[1]] == masknoData:
                continue
            if vis[centerCell[0], centerCell[1]] == 1:
                continue

            # print(centerCell)
            # 以此栅格为中心栅格，寻找最大构建坡度的方向
            maxSlope = -1
            maxSlopeDir = -1
            for kk in range(8):
                nextCell = (centerCell[0] + dmove[kk][0], centerCell[1] + dmove[kk][1])
                if not check_boundary(nextCell[0], nextCell[1], row, col):
                    continue
                if mask[nextCell[0], nextCell[1]] == masknoData:
                    continue
                if vis[nextCell[0], nextCell[1]] == 0:
                    continue
                if kk in [0, 2, 4, 6]:
                    factor = 1
                else:
                    factor = math.sqrt(2)
                eDis = math.sqrt(abs(centerCell[0] - startcell[0]) * abs(centerCell[0] - startcell[0]) + abs(
                    centerCell[1] - startcell[1]) * abs(centerCell[1] - startcell[1]))
                fDis = fDiss[nextCell[0], nextCell[1]] + factor
                # print(eDis,fDis)
                likely = eDis / fDis
                # print('out',cell, centerCell, nextCell, likely)
                if likely > maxSlope:
                    # print(cell,centerCell,nextCell,likely)
                    maxSlope = likely
                    maxSlopeDir = 2 ** kk
            # print(maxSlope)
            weight[centerCell[0], centerCell[1]] = maxSlope
            newDir[centerCell[0], centerCell[1]] = maxSlopeDir
            vis[centerCell[0], centerCell[1]] = 1

    # print(newDir)
    return newDir


def iterflat(mask, masknoData, startcell, newDir):
    row, col = mask.shape

    vis = np.zeros((row, col))
    vis[startcell[0], startcell[1]] = 1

    popCells = [startcell]
    while popCells:
        popCell = popCells.pop()
        for k in range(8):
            nextCell = (popCell[0] + dmove[k][0], popCell[1] + dmove[k][1])
            if not check_boundary(nextCell[0], nextCell[1], row, col):
                continue
            if mask[nextCell[0], nextCell[1]] == masknoData:
                continue
            if vis[nextCell[0], nextCell[1]] == 1:
                continue

            newDir[nextCell[0], nextCell[1]] = 2 ** ((4 + k) % 8)
            vis[nextCell[0], nextCell[1]] = 1
            popCells.insert(0, nextCell)

    return newDir


def iterflat_pool(mask, masknoData, startcell, newDir, proj, geo, outpath, sinkId):
    row, col = mask.shape

    vis = np.zeros((row, col))
    vis[startcell[0], startcell[1]] = 1

    popCells = [startcell]
    while popCells:
        popCell = popCells.pop()
        for k in range(8):
            nextCell = (popCell[0] + dmove[k][0], popCell[1] + dmove[k][1])
            if not check_boundary(nextCell[0], nextCell[1], row, col):
                continue
            if mask[nextCell[0], nextCell[1]] == masknoData:
                continue
            if vis[nextCell[0], nextCell[1]] == 1:
                continue

            newDir[nextCell[0], nextCell[1]] = 2 ** ((4 + k) % 8)
            vis[nextCell[0], nextCell[1]] = 1
            popCells.insert(0, nextCell)

    Raster.save_raster(os.path.join(outpath, str(sinkId) + '.tif'), newDir, proj, geo, gdal.GDT_Byte, 0)


def queuePriorityFlood(mask, maskNodata, dem, startcell, newDir):
    """
    使用python库构建的栈和队列，加快速度
    :param mask:
    :param maskNodata:
    :param dem:
    :param startcell:
    :param newdir:
    :return:
    """

    row, col = mask.shape
    Vis = np.zeros((row, col))

    quickQ = genral_functions.QueueUsingQueue()
    maxQ = genral_functions.MaxPriorityQueue()
    minQ = genral_functions.MinPriorityQueue()

    maxQ.push(startcell[0], startcell[1], dem[startcell[0], startcell[1]])
    Vis[startcell[0], startcell[1]] = 1

    while (not quickQ.is_empty()) or (not maxQ.is_empty()) or (not minQ.is_empty()):

        while not quickQ.is_empty():
            # 先把中心栅格拿出来，
            tempOutCell = quickQ.dequeue()
            # print(tempOutCell)
            outCell = [tempOutCell]
            while outCell:
                popCell = outCell.pop()

                # print('Flat', popCell)
                popCellH = popCell[2]
                for k in range(8):
                    nextCell = (popCell[0] + dmove[k][0], popCell[1] + dmove[k][1])
                    if not check_boundary(nextCell[0], nextCell[1], row, col):
                        continue
                    if mask[nextCell[0], nextCell[1]] == maskNodata:
                        continue
                    nextCellH = dem[nextCell[0], nextCell[1]]

                    if Vis[nextCell[0], nextCell[1]] == 1:
                        continue
                    # 是，给流向加进去(分上下坡和平地加)，Vis赋值1；不是跳过
                    if nextCellH == popCellH:
                        newDir[nextCell[0], nextCell[1]] = 2 ** ((k + 4) % 8)
                        Vis[nextCell[0], nextCell[1]] = 1
                        outCell.insert(0, (nextCell[0],nextCell[1],nextCellH))
                    else:
                        # 如果不是平地也要加入相应的队列
                        if nextCellH < popCellH:
                            # 下坡,D8
                            newDir[nextCell[0], nextCell[1]] = reverseD8(nextCell, dem, Vis,mask)
                            maxQ.push(nextCell[0], nextCell[1], nextCellH)
                            Vis[nextCell[0], nextCell[1]] = 1
                        else:
                            # 上坡
                            minQ.push(nextCell[0], nextCell[1], nextCellH)
                            newDir[nextCell[0], nextCell[1]] = d8(nextCell, dem, Vis,mask)
                            Vis[nextCell[0], nextCell[1]] = 1
        # 下坡流向
        flag = False
        while not maxQ.is_empty():
            tempOutCell = maxQ.pop()
            tempOutCellH = tempOutCell[2]

            for k in range(8):
                nextCell = (tempOutCell[0] + dmove[k][0], tempOutCell[1] + dmove[k][1])
                if not check_boundary(nextCell[0], nextCell[1], row, col):
                    continue
                if mask[nextCell[0], nextCell[1]] == maskNodata:
                    continue
                if Vis[nextCell[0], nextCell[1]] == 1:
                    continue
                nextCellH = dem[nextCell[0], nextCell[1]]
                if nextCellH == tempOutCellH:
                    newDir[nextCell[0], nextCell[1]] = 2 ** ((k + 4) % 8)
                    Vis[nextCell[0], nextCell[1]] = 1
                    quickQ.enqueue(nextCell[0],nextCell[1],nextCellH)
                    # flatFlag = True
                    # break
                else:
                    # 如果不是平地也要加入相应的队列
                    if nextCellH < tempOutCellH:
                        # 下坡,D8
                        newDir[nextCell[0], nextCell[1]] = reverseD8(nextCell, dem, Vis,mask)
                        maxQ.push(nextCell[0],nextCell[1],nextCellH)
                        Vis[nextCell[0], nextCell[1]] = 1
                    else:
                        # 上坡
                        minQ.push(nextCell[0],nextCell[1],nextCellH)
                        newDir[nextCell[0], nextCell[1]] = d8(nextCell, dem, Vis,mask)
                        Vis[nextCell[0], nextCell[1]] = 1
            flag = True
            break  # 每次只进行一个中心栅格的邻域流向
        if flag:
            continue
        # 上坡流向
        while not minQ.is_empty():
            # 先把中心栅格拿出来，
            tempOutCell = minQ.pop()

            tempOutCellH = tempOutCell[2]

            for k in range(8):
                nextCell = (tempOutCell[0] + dmove[k][0], tempOutCell[1] + dmove[k][1])
                if not check_boundary(nextCell[0], nextCell[1], row, col):
                    continue
                if mask[nextCell[0], nextCell[1]] == maskNodata:
                    continue
                if Vis[nextCell[0], nextCell[1]] == 1:
                    continue
                nextCellH = dem[nextCell[0], nextCell[1]]
                if nextCellH == tempOutCellH:
                    newDir[nextCell[0], nextCell[1]] = 2 ** ((k + 4) % 8)
                    Vis[nextCell[0], nextCell[1]] = 1
                    quickQ.enqueue(nextCell[0],nextCell[1],nextCellH)
                    # break
                else:
                    # 如果不是平地也要加入相应的队列
                    if nextCellH < tempOutCellH:
                        # 下坡,D8
                        newDir[nextCell[0], nextCell[1]] = reverseD8(nextCell, dem, Vis,mask)
                        maxQ.push(nextCell[0],nextCell[1],nextCellH)
                        Vis[nextCell[0], nextCell[1]] = 1
                    else:
                        # 上坡
                        minQ.push(nextCell[0],nextCell[1],nextCellH)
                        newDir[nextCell[0], nextCell[1]] = d8(nextCell, dem, Vis,mask)
                        Vis[nextCell[0], nextCell[1]] = 1
            break  # 每次只进行一个中心栅格的邻域流向
    return newDir
def temp_queuePriorityFlood(mask, maskNodata, dem, startcell, newDir):
    """
    使用python库构建的栈和队列，加快速度
    :param mask:dem
    :param maskNodata:
    :param dem:
    :param startcell:
    :param newdir:
    :return:
    """

    row, col = mask.shape
    Vis = np.zeros((row, col))

    quickQ = genral_functions.QueueUsingQueue()
    maxQ = genral_functions.MaxPriorityQueue()
    minQ = genral_functions.MinPriorityQueue()

    maxQ.push(startcell[0], startcell[1], dem[startcell[0], startcell[1]])
    Vis[startcell[0], startcell[1]] = 1

    while (not quickQ.is_empty()) or (not maxQ.is_empty()) or (not minQ.is_empty()):

        while not quickQ.is_empty():
            # 先把中心栅格拿出来，
            tempOutCell = quickQ.dequeue()
            # print(tempOutCell)
            outCell = [tempOutCell]
            while outCell:
                popCell = outCell.pop()

                # print('Flat', popCell)
                popCellH = popCell[2]
                for k in range(8):
                    nextCell = (popCell[0] + dmove[k][0], popCell[1] + dmove[k][1])
                    if not check_boundary(nextCell[0], nextCell[1], row, col):
                        continue
                    if mask[nextCell[0], nextCell[1]] == maskNodata:
                        continue
                    nextCellH = dem[nextCell[0], nextCell[1]]

                    if Vis[nextCell[0], nextCell[1]] == 1:
                        continue
                    # 是，给流向加进去(分上下坡和平地加)，Vis赋值1；不是跳过
                    if nextCellH == popCellH:
                        newDir[nextCell[0], nextCell[1]] = 2 ** ((k + 4) % 8)
                        Vis[nextCell[0], nextCell[1]] = 1
                        outCell.insert(0, (nextCell[0],nextCell[1],nextCellH))
                    else:
                        # 如果不是平地也要加入相应的队列
                        if nextCellH < popCellH:
                            # 下坡,D8
                            newDir[nextCell[0], nextCell[1]] = reverseD8(nextCell, dem, Vis,mask)
                            maxQ.push(nextCell[0], nextCell[1], nextCellH)
                            Vis[nextCell[0], nextCell[1]] = 1
                        else:
                            # 上坡
                            minQ.push(nextCell[0], nextCell[1], nextCellH)
                            newDir[nextCell[0], nextCell[1]] = d8(nextCell, dem, Vis,mask)
                            Vis[nextCell[0], nextCell[1]] = 1
        # 下坡流向
        flag = False
        while not maxQ.is_empty():
            tempOutCell = maxQ.pop()
            tempOutCellH = tempOutCell[2]

            for k in range(8):
                nextCell = (tempOutCell[0] + dmove[k][0], tempOutCell[1] + dmove[k][1])
                if not check_boundary(nextCell[0], nextCell[1], row, col):
                    continue
                if mask[nextCell[0], nextCell[1]] == maskNodata:
                    continue
                if Vis[nextCell[0], nextCell[1]] == 1:
                    continue
                nextCellH = dem[nextCell[0], nextCell[1]]
                if nextCellH == tempOutCellH:
                    newDir[nextCell[0], nextCell[1]] = 2 ** ((k + 4) % 8)
                    Vis[nextCell[0], nextCell[1]] = 1
                    quickQ.enqueue(nextCell[0],nextCell[1],nextCellH)
                    # flatFlag = True
                    # break
                else:
                    # 如果不是平地也要加入相应的队列
                    if nextCellH < tempOutCellH:
                        # 下坡,D8
                        newDir[nextCell[0], nextCell[1]] = reverseD8(nextCell, dem, Vis,mask)
                        maxQ.push(nextCell[0],nextCell[1],nextCellH)
                        Vis[nextCell[0], nextCell[1]] = 1
                    else:
                        # 上坡
                        minQ.push(nextCell[0],nextCell[1],nextCellH)
                        newDir[nextCell[0], nextCell[1]] = d8(nextCell, dem, Vis,mask)
                        Vis[nextCell[0], nextCell[1]] = 1
            flag = True
            break  # 每次只进行一个中心栅格的邻域流向
        if flag:
            continue
        # 上坡流向
        while not minQ.is_empty():
            # 先把中心栅格拿出来，
            tempOutCell = minQ.pop()

            tempOutCellH = tempOutCell[2]

            for k in range(8):
                nextCell = (tempOutCell[0] + dmove[k][0], tempOutCell[1] + dmove[k][1])
                if not check_boundary(nextCell[0], nextCell[1], row, col):
                    continue
                if mask[nextCell[0], nextCell[1]] == maskNodata:
                    continue
                if Vis[nextCell[0], nextCell[1]] == 1:
                    continue
                nextCellH = dem[nextCell[0], nextCell[1]]
                if nextCellH == tempOutCellH:
                    newDir[nextCell[0], nextCell[1]] = 2 ** ((k + 4) % 8)
                    Vis[nextCell[0], nextCell[1]] = 1
                    quickQ.enqueue(nextCell[0],nextCell[1],nextCellH)
                    # break
                else:
                    # 如果不是平地也要加入相应的队列
                    if nextCellH < tempOutCellH:
                        # 下坡,D8
                        newDir[nextCell[0], nextCell[1]] = reverseD8(nextCell, dem, Vis,mask)
                        maxQ.push(nextCell[0],nextCell[1],nextCellH)
                        Vis[nextCell[0], nextCell[1]] = 1
                    else:
                        # 上坡
                        minQ.push(nextCell[0],nextCell[1],nextCellH)
                        newDir[nextCell[0], nextCell[1]] = d8(nextCell, dem, Vis,mask)
                        Vis[nextCell[0], nextCell[1]] = 1
            break  # 每次只进行一个中心栅格的邻域流向
    return newDir

def priorityFlood(mask, maskNodata, dem, startcell, newDir):
    """
    从startcell出发，开始使用PF算法重新计算流向.输入的数据都是与mask匹配的
    :param mask:
    :param dem:
    :param startcell:
    :return:
    """
    # print(mask.shape,dem.shape,newDir.shape)
    row, col = mask.shape

    Vis = np.zeros((row, col))
    # newDir = np.zeros((row,col))
    quickQ = np.zeros((row, col))
    maxQ = np.zeros((row, col))
    minQ = np.zeros((row, col))

    maxQ[startcell[0], startcell[1]] = dem[startcell[0], startcell[1]]
    print(maxQ[startcell[0], startcell[1]])

    Vis[startcell[0], startcell[1]] = 1

    while (not np.all(quickQ == 0)) or (not np.all(maxQ == 0)) or (not np.all(minQ == 0)):

        while not np.all(quickQ == 0):

            # 先把中心栅格拿出来，
            tempCells = np.where(quickQ != 0)

            tempOutCell = (tempCells[0][0], tempCells[1][0])
            # print('Quick',tempOutCell)
            outCell = [tempOutCell]
            quickQ[tempOutCell[0], tempOutCell[1]] = 0  # 把中心栅格删掉

            # 迭代平地
            while outCell:

                popCell = outCell.pop()

                # print('Flat', popCell)
                popCellH = dem[popCell[0], popCell[1]]
                for k in range(8):
                    nextCell = (popCell[0] + dmove[k][0], popCell[1] + dmove[k][1])
                    if not check_boundary(nextCell[0], nextCell[1], row, col):
                        continue
                    if mask[nextCell[0], nextCell[1]] == maskNodata:
                        continue
                    nextCellH = dem[nextCell[0], nextCell[1]]

                    if Vis[nextCell[0], nextCell[1]] == 1:
                        continue
                    # 是，给流向加进去(分上下坡和平地加)，Vis赋值1；不是跳过
                    if nextCellH == popCellH:
                        newDir[nextCell[0], nextCell[1]] = 2 ** ((k + 4) % 8)

                        Vis[nextCell[0], nextCell[1]] = 1
                        outCell.insert(0, nextCell)
                    else:
                        # 如果不是平地也要加入相应的队列
                        if nextCellH < popCellH:
                            # 下坡,D8
                            newDir[nextCell[0], nextCell[1]] = reverseD8(nextCell, dem, Vis,mask)
                            maxQ[nextCell[0], nextCell[1]] = nextCellH
                            Vis[nextCell[0], nextCell[1]] = 1
                        else:
                            # 上坡
                            minQ[nextCell[0], nextCell[1]] = nextCellH
                            newDir[nextCell[0], nextCell[1]] = d8(nextCell, dem, Vis,mask)
                            Vis[nextCell[0], nextCell[1]] = 1

        # 下坡流向
        flag = False
        flatFlag = False  # 如果遇到平地要立即进入快速栅格
        while not np.all(maxQ == 0):

            # 先把中心栅格拿出来，
            tempCells = np.where(maxQ == maxQ.max())


            tempOutCell = (tempCells[0][0], tempCells[1][0])
            # print(tempOutCell)
            # print('Down', tempOutCell)
            # outCell = [tempOutCell]
            maxQ[tempOutCell[0], tempOutCell[1]] = 0  # 把中心栅格删掉
            tempOutCellH = dem[tempOutCell[0], tempOutCell[1]]
            # print('Max',tempOutCell)
            for k in range(8):
                nextCell = (tempOutCell[0] + dmove[k][0], tempOutCell[1] + dmove[k][1])
                if not check_boundary(nextCell[0], nextCell[1], row, col):
                    continue
                if mask[nextCell[0], nextCell[1]] == maskNodata:
                    continue
                if Vis[nextCell[0], nextCell[1]] == 1:
                    continue
                nextCellH = dem[nextCell[0], nextCell[1]]
                if nextCellH == tempOutCellH:
                    newDir[nextCell[0], nextCell[1]] = 2 ** ((k + 4) % 8)
                    Vis[nextCell[0], nextCell[1]] = 1
                    quickQ[nextCell[0], nextCell[1]] = nextCellH

                    # flatFlag = True
                    # break
                else:
                    # 如果不是平地也要加入相应的队列
                    if nextCellH < tempOutCellH:
                        # 下坡,D8
                        newDir[nextCell[0], nextCell[1]] = reverseD8(nextCell, dem, Vis,mask)
                        maxQ[nextCell[0], nextCell[1]] = nextCellH
                        Vis[nextCell[0], nextCell[1]] = 1
                    else:
                        # 上坡
                        minQ[nextCell[0], nextCell[1]] = nextCellH
                        newDir[nextCell[0], nextCell[1]] = d8(nextCell, dem, Vis,mask)
                        Vis[nextCell[0], nextCell[1]] = 1

            flag = True
            break  # 每次只进行一个中心栅格的邻域流向
        if flag:
            continue
        # 上坡流向
        while not np.all(minQ == 0):
            # 先把中心栅格拿出来，

            tempCells = np.where(minQ == minQ[minQ != 0].min())
            tempOutCell = (tempCells[0][0], tempCells[1][0])
            # print('Up', tempOutCell)
            # outCell = [tempOutCell]
            minQ[tempOutCell[0], tempOutCell[1]] = 0  # 把中心栅格删掉
            tempOutCellH = dem[tempOutCell[0], tempOutCell[1]]
            # print('Min',tempOutCell)
            for k in range(8):
                nextCell = (tempOutCell[0] + dmove[k][0], tempOutCell[1] + dmove[k][1])
                if not check_boundary(nextCell[0], nextCell[1], row, col):
                    continue
                if mask[nextCell[0], nextCell[1]] == maskNodata:
                    continue
                if Vis[nextCell[0], nextCell[1]] == 1:
                    continue
                nextCellH = dem[nextCell[0], nextCell[1]]
                if nextCellH == tempOutCellH:
                    newDir[nextCell[0], nextCell[1]] = 2 ** ((k + 4) % 8)
                    Vis[nextCell[0], nextCell[1]] = 1
                    quickQ[nextCell[0], nextCell[1]] = nextCellH
                    # break
                else:
                    # 如果不是平地也要加入相应的队列
                    if nextCellH < tempOutCellH:
                        # 下坡,D8
                        newDir[nextCell[0], nextCell[1]] = reverseD8(nextCell, dem, Vis,mask)
                        maxQ[nextCell[0], nextCell[1]] = nextCellH
                        Vis[nextCell[0], nextCell[1]] = 1

                    else:
                        # 上坡
                        minQ[nextCell[0], nextCell[1]] = nextCellH
                        newDir[nextCell[0], nextCell[1]] = d8(nextCell, dem, Vis,mask)
                        Vis[nextCell[0], nextCell[1]] = 1
            break  # 每次只进行一个中心栅格的邻域流向

    return newDir


def priorityFlood_special(mask, maskNodata, dem, startcell, newDir):
    """
    从startcell出发，开始使用PF算法重新计算流向.输入的数据都是与mask匹配的.

    有低于 0的值

    :param mask:
    :param dem:
    :param startcell:
    :return:
    """
    # print(mask.shape,dem.shape,newDir.shape)
    row, col = mask.shape

    Vis = np.zeros((row, col))
    # newDir = np.zeros((row,col))
    quickQ = np.zeros((row, col))
    quickQ[:,:] = -9999
    maxQ = np.zeros((row, col))
    maxQ[:,:] = -9999
    minQ = np.zeros((row, col))
    minQ[:,:] = -9999

    maxQ[startcell[0], startcell[1]] = dem[startcell[0], startcell[1]]

    Vis[startcell[0], startcell[1]] = 1

    while (not np.all(quickQ == -9999)) or (not np.all(maxQ == -9999)) or (not np.all(minQ == -9999)):

        while not np.all(quickQ == -9999):

            # 先把中心栅格拿出来，
            tempCells = np.where(quickQ != -9999)


            tempOutCell = (tempCells[0][0], tempCells[1][0])
            # print(tempOutCell)
            # print('Quick',tempOutCell)
            outCell = [tempOutCell]
            quickQ[tempOutCell[0], tempOutCell[1]] = -9999  # 把中心栅格删掉

            # 迭代平地
            while outCell:

                popCell = outCell.pop()

                # print('Flat', popCell)
                popCellH = dem[popCell[0], popCell[1]]
                for k in range(8):
                    nextCell = (popCell[0] + dmove[k][0], popCell[1] + dmove[k][1])
                    if not check_boundary(nextCell[0], nextCell[1], row, col):
                        continue
                    if mask[nextCell[0], nextCell[1]] == maskNodata:
                        continue
                    nextCellH = dem[nextCell[0], nextCell[1]]

                    if Vis[nextCell[0], nextCell[1]] == 1:
                        continue
                    # 是，给流向加进去(分上下坡和平地加)，Vis赋值1；不是跳过
                    if nextCellH == popCellH:
                        newDir[nextCell[0], nextCell[1]] = 2 ** ((k + 4) % 8)

                        Vis[nextCell[0], nextCell[1]] = 1
                        outCell.insert(0, nextCell)
                    else:
                        # 如果不是平地也要加入相应的队列
                        if nextCellH < popCellH:
                            # 下坡,D8
                            newDir[nextCell[0], nextCell[1]] = reverseD8(nextCell, dem, Vis,mask)
                            maxQ[nextCell[0], nextCell[1]] = nextCellH
                            Vis[nextCell[0], nextCell[1]] = 1
                        else:
                            # 上坡
                            minQ[nextCell[0], nextCell[1]] = nextCellH
                            newDir[nextCell[0], nextCell[1]] = d8(nextCell, dem, Vis,mask)
                            Vis[nextCell[0], nextCell[1]] = 1

        # 下坡流向
        flag = False
        flatFlag = False  # 如果遇到平地要立即进入快速栅格
        while not np.all(maxQ == -9999):

            # 先把中心栅格拿出来，
            tempCells = np.where(maxQ == maxQ.max())

            tempOutCell = (tempCells[0][0], tempCells[1][0])
            # print(tempOutCell)
            # print('Down', tempOutCell)
            # outCell = [tempOutCell]
            maxQ[tempOutCell[0], tempOutCell[1]] = -9999  # 把中心栅格删掉
            tempOutCellH = dem[tempOutCell[0], tempOutCell[1]]
            # print('Max',tempOutCell)
            for k in range(8):
                nextCell = (tempOutCell[0] + dmove[k][0], tempOutCell[1] + dmove[k][1])
                if not check_boundary(nextCell[0], nextCell[1], row, col):
                    continue
                if mask[nextCell[0], nextCell[1]] == maskNodata:
                    continue
                if Vis[nextCell[0], nextCell[1]] == 1:
                    continue
                nextCellH = dem[nextCell[0], nextCell[1]]
                if nextCellH == tempOutCellH:
                    newDir[nextCell[0], nextCell[1]] = 2 ** ((k + 4) % 8)
                    Vis[nextCell[0], nextCell[1]] = 1
                    quickQ[nextCell[0], nextCell[1]] = nextCellH

                    # flatFlag = True
                    # break
                else:
                    # 如果不是平地也要加入相应的队列
                    if nextCellH < tempOutCellH:
                        # 下坡,D8
                        newDir[nextCell[0], nextCell[1]] = reverseD8(nextCell, dem, Vis,mask)
                        maxQ[nextCell[0], nextCell[1]] = nextCellH
                        Vis[nextCell[0], nextCell[1]] = 1
                    else:
                        # 上坡
                        minQ[nextCell[0], nextCell[1]] = nextCellH
                        newDir[nextCell[0], nextCell[1]] = d8(nextCell, dem, Vis,mask)
                        Vis[nextCell[0], nextCell[1]] = 1

            flag = True
            break  # 每次只进行一个中心栅格的邻域流向
        if flag:
            continue
        # 上坡流向
        while not np.all(minQ == -9999):
            # 先把中心栅格拿出来，

            tempCells = np.where(minQ == minQ[minQ != -9999].min())
            tempOutCell = (tempCells[0][0], tempCells[1][0])
            # print('Up', tempOutCell)
            # outCell = [tempOutCell]
            minQ[tempOutCell[0], tempOutCell[1]] = -9999  # 把中心栅格删掉
            tempOutCellH = dem[tempOutCell[0], tempOutCell[1]]
            # print('Min',tempOutCell)
            for k in range(8):
                nextCell = (tempOutCell[0] + dmove[k][0], tempOutCell[1] + dmove[k][1])
                if not check_boundary(nextCell[0], nextCell[1], row, col):
                    continue
                if mask[nextCell[0], nextCell[1]] == maskNodata:
                    continue
                if Vis[nextCell[0], nextCell[1]] == 1:
                    continue
                nextCellH = dem[nextCell[0], nextCell[1]]
                if nextCellH == tempOutCellH:
                    newDir[nextCell[0], nextCell[1]] = 2 ** ((k + 4) % 8)
                    Vis[nextCell[0], nextCell[1]] = 1
                    quickQ[nextCell[0], nextCell[1]] = nextCellH
                    # break
                else:
                    # 如果不是平地也要加入相应的队列
                    if nextCellH < tempOutCellH:
                        # 下坡,D8
                        newDir[nextCell[0], nextCell[1]] = reverseD8(nextCell, dem, Vis,mask)
                        maxQ[nextCell[0], nextCell[1]] = nextCellH
                        Vis[nextCell[0], nextCell[1]] = 1

                    else:
                        # 上坡
                        minQ[nextCell[0], nextCell[1]] = nextCellH
                        newDir[nextCell[0], nextCell[1]] = d8(nextCell, dem, Vis,mask)
                        Vis[nextCell[0], nextCell[1]] = 1
            break  # 每次只进行一个中心栅格的邻域流向

    return newDir


def main1(outVenu):
    # 整体代码运行
    baseName = os.path.basename(outVenu)
    db_name = os.path.join(outVenu, baseName + "_burndb_.db")
    dem_file = os.path.join(outVenu, baseName + "_burnDEM_.tif")
    dir_file = os.path.join(outVenu, baseName + "_burndir_.tif")
    sink_file = os.path.join(outVenu, baseName + "_burnsink_.tif")
    occ_file = os.path.join(outVenu, baseName + "_Occ_.tif")
    outDir_file = os.path.join(outVenu, baseName + "_outDIR_.tif")
    dem = Raster.get_raster(dem_file)
    dir = Raster.get_raster(dir_file)
    sink = Raster.get_raster(sink_file)
    dir = np.array(dir, dtype=np.int64)
    dir[dir == 255] = 0
    Occ = Raster.get_raster(occ_file)
    proj, geo, sinkNodata = Raster.get_proj_geo_nodata(sink_file)
    proj, geo, dirNodata = Raster.get_proj_geo_nodata(dir_file)
    proj, geo, demNodata = Raster.get_proj_geo_nodata(dem_file)
    _, _, O_Nodata = Raster.get_proj_geo_nodata(occ_file)
    # print('Search the flat sink from the database')
    # infos = query_data_byType(db_name, str(1))
    # for info in infos:
    #     # print(info)
    #     sinkId = info[0]
    #     outModifiedDirPath = os.path.join(outVenu, str(sinkId) + '.tif')
    #     if os.path.exists(outModifiedDirPath):
    #         continue
    #
    #     minX = info[-5]
    #     maxX = info[-4]
    #     minY = info[-3]
    #     maxY = info[-2]
    #     outCell = (info[3], info[4])
    #     outCellDir = dir[outCell[0], outCell[1]]
    #
    #     # 掩膜sink
    #     # print(minX,maxX,minY,maxY)
    #     maskDem = dem[minX:maxX + 1, minY:maxY + 1].copy()
    #     maskSink = sink[minX:maxX + 1, minY:maxY + 1].copy()
    #     maskSink[maskSink == sinkId] = 1
    #     maskSink[maskSink != 1] = 0
    #     row, col = maskSink.shape
    #
    #     # 记录出口点数据
    #     maskDir = dir[minX:maxX + 1, minY:maxY + 1].copy()
    #     maskDir[maskSink == 1] = 0
    #     maskDir[outCell[0] - minX, outCell[1] - minY] = outCellDir
    #
    #     A = iterflat(maskSink, 0, (outCell[0] - minX, outCell[1] - minY), maskDir)  # Python
    #     # A = C_iterFlat(maskSink,0,outCell[0] - minX,outCell[1] - minY,row,col)   # C语言
    #     maskDir[outCell[0] - minX, outCell[1] - minY] = outCellDir
    #     # 单独存储修正你后的流向
    #     # temp_geo = (geo[0] + geo[1] * minY, geo[1], geo[2], geo[3] + geo[5] * minX, geo[4], geo[5])
    #     # outModifiedDirPath = os.path.join(outVenu, str(sinkId) + '.tif')
    #     # Raster.save_raster(outModifiedDirPath, A, proj, temp_geo, gdal.GDT_Byte, 0)
    #
    #     # 写入原来的流向
    #     # C_updataDir(dir,minX,minY,A)
    #     for ti in range(row):
    #         for tj in range(col):
    #             if A[ti, tj] == 255:
    #                 continue
    #             dir[minX + ti, minY + tj] = A[ti, tj]
    #
    # print('Search the fake sink from the database')
    # infos = query_data_byType(db_name, str(3))
    # # print(infos)
    # nn = 0
    # for info in infos:
    #     # print(info)
    #     sinkId = info[0]
    #
    #     minX = info[-5]
    #     maxX = info[-4]
    #     minY = info[-3]
    #     maxY = info[-2]
    #     outCell = (info[3], info[4])
    #     outCellDir = dir[outCell[0], outCell[1]]
    #     # print(info[0])
    #     # 掩膜sink
    #     # print(minX,maxX,minY,maxY)
    #     # print(sinkId)
    #
    #     maskSink = sink[minX:maxX + 1, minY:maxY + 1].copy()
    #     maskSink[maskSink == sinkId] = 1
    #     maskSink[maskSink != 1] = 0  # 写成0 是因为好用：否则后面查最小值的时候会错误
    #     row, col = maskSink.shape
    #     maskDem = dem[minX:maxX + 1, minY:maxY + 1].copy()
    #     maskDem[maskSink == 0] = -9999
    #     # 记录出口点数据
    #     maskDir = dir[minX:maxX + 1, minY:maxY + 1].copy()  # np.zeros((row,col))
    #     maskDir[maskSink == 1] = 0
    #     maskDir[outCell[0] - minX, outCell[1] - minY] = outCellDir
    #
    #     # 计算淹没内的流向
    #     # A = priorityFlood(maskSink, 0, maskDem, (outCell[0] - minX, outCell[1] - minY), maskDir)  # Python
    #     A = queuePriorityFlood(maskSink, 0, maskDem, (outCell[0] - minX, outCell[1] - minY), maskDir)  # Python
    #     # A = C_PriorityFlood(maskSink, 0, maskDem, outCell[0] - minX, outCell[1] - minY,row,col)  # C语言
    #
    #     # 单独存储修正你后的流向
    #     temp_geo = (geo[0] + geo[1] * minY, geo[1], geo[2], geo[3] + geo[5] * minX, geo[4], geo[5])
    #
    #     # Raster.save_raster(outModifiedDirPath,A,proj,temp_geo,gdal.GDT_Byte,0)
    #
    #     # 写入原来的流向
    #     # C_updataDir(dir, minX, minY, A)
    #     for ti in range(row):
    #         for tj in range(col):
    #             if A[ti, tj] == 255:
    #                 continue
    #             dir[minX + ti, minY + tj] = A[ti, tj]

    ############################## 内流区特殊处理，需准备mask，fdir，dem ###########################################
    print('Search the endorheic from the database')
    infos = query_data_byType(db_name, str(2))
    # print(infos)

    for info in infos:
        # print(info)
        sinkId = info[0]
        minX = info[-5]
        maxX = info[-4]
        minY = info[-3]
        maxY = info[-2]
        outCell = (info[5], info[6])
        # outCellDir = dir[outCell[0], outCell[1]]
        # 掩膜sink
        maskDem = dem[minX:maxX + 1, minY:maxY + 1].copy()
        maskSink = sink[minX:maxX + 1, minY:maxY + 1].copy()
        maskSink[maskSink == sinkId] = 1
        maskSink[maskSink != 1] = 0  # 写成0 是因为好用：否则后面查最小值的时候会错误
        maskDem[maskSink == 0] = demNodata
        maskOcc = Occ[minX:maxX + 1, minY:maxY + 1].copy()
        maskOcc[maskSink == 0] = O_Nodata
        row, col = maskSink.shape

        # 记录出口点数据
        maskDir = dir[minX:maxX + 1, minY:maxY + 1].copy()  # np.zeros((row,col))
        maskDir[maskSink == 0] = 255
        # maskDir[outCell[0] - minX, outCell[1] - minY] = outCellDir

        # 计算淹没内的流向
        # A = priorityFlood(maskSink, 0, maskDem, (outCell[0] - minX, outCell[1] - minY), maskDir)  # Python
        # A = C_PriorityFlood(maskSink, 0, maskDem, outCell[0] - minX, outCell[1] - minY, row, col)  # C语言
        # 单独存储修正你后的流向

        # 存储掩膜文件
        temp_geo = (geo[0] + geo[1] * minY, geo[1], geo[2], geo[3] + geo[5] * minX, geo[4], geo[5])
        endor = os.path.join(outVenu, 'endorheic')
        if not os.path.exists(endor):
            os.mkdir(endor)
            os.chmod(endor, 0o777)
        sinkIdDir = os.path.join(endor, str(sinkId))
        if not os.path.exists(sinkIdDir):
            os.mkdir(sinkIdDir)
            os.chmod(sinkIdDir, 0o777)

        outmaskSinkPath = os.path.join(sinkIdDir, str(sinkId) + '_mask_.tif')
        outmaskDemPath = os.path.join(sinkIdDir, str(sinkId) + '_dem_.tif')
        outmaskDirPath = os.path.join(sinkIdDir, str(sinkId) + '_dir_.tif')
        outmaskOccPath = os.path.join(sinkIdDir, str(sinkId) + '_occ_.tif')
        Raster.save_raster(outmaskOccPath, maskOcc, proj, temp_geo, gdal.GDT_Float32, O_Nodata)
        Raster.save_raster(outmaskSinkPath, maskSink, proj, temp_geo, gdal.GDT_Byte, 0)
        Raster.save_raster(outmaskDirPath, maskDir, proj, temp_geo, gdal.GDT_Byte, 255)
        Raster.save_raster(outmaskDemPath, maskDem, proj, temp_geo, gdal.GDT_Float32, demNodata)
        # 存储掩膜文件

        # outModifiedDirPath = os.path.join(outVenu, str(sinkId) + '.tif')
        # Raster.save_raster(outModifiedDirPath, A, proj, temp_geo, gdal.GDT_Byte, 0)

        # 写入原来的流向
        # C_updataDir(dir, minX, minY, A)
        # for ti in range(row):
        #     for tj in range(col):
        #         if A[ti,tj] == 255:
        #             continue
        #         dir[minX+ti,minY+tj]=A[ti,tj]

    # outSink = r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\20240923\Sink2.tif'
    Raster.save_raster(outDir_file, dir, proj, geo, gdal.GDT_Byte, 0)
    # Raster.save_raster(outSink, sink, proj, geo, gdal.GDT_Float32, sinkNodata)
def exorheic_process(outVenu):
    """
    原main1的前半部分
    :param outVenu: NA 下的每个流域文件夹
    :return:
    """

    # 整体代码运行
    baseName = os.path.basename(outVenu)
    db_name = os.path.join(outVenu, baseName + "_burndb.db")
    if not os.path.exists(db_name):
        print(outVenu,'{:s} not exists'.format(db_name))
        return
    if os.path.exists(os.path.join(outVenu, baseName + "_ModifiedDir.tif")):
        return
    dem_file = os.path.join(outVenu, baseName + "_burnDEM.tif")
    dir_file = os.path.join(outVenu, baseName + "_burndir.tif")
    sink_file = os.path.join(outVenu, baseName + "_burnsink.tif")
    # occ_file = os.path.join(outVenu, baseName + "_Occ_.tif")
    outDir_file = os.path.join(outVenu, baseName + "_outDIR_.tif")
    dem = Raster.get_raster(dem_file)
    dir = Raster.get_raster(dir_file)
    sink = Raster.get_raster(sink_file)
    dir = np.array(dir, dtype=np.int64)
    dir[dir == 255] = 0
    # Occ = Raster.get_raster(occ_file)
    proj, geo, sinkNodata = Raster.get_proj_geo_nodata(sink_file)
    proj, geo, dirNodata = Raster.get_proj_geo_nodata(dir_file)
    proj, geo, demNodata = Raster.get_proj_geo_nodata(dem_file)
    # _, _, O_Nodata = Raster.get_proj_geo_nodata(occ_file)
    print('Search the flat sink from the database')
    infos = query_data_byType(db_name, str(1))
    for info in infos:
        # print(info)
        sinkId = info[0]
        outModifiedDirPath = os.path.join(outVenu, str(sinkId) + '.tif')
        if os.path.exists(outModifiedDirPath):
            continue

        minX = info[-5]
        maxX = info[-4]
        minY = info[-3]
        maxY = info[-2]
        outCell = (info[3], info[4])
        outCellDir = dir[outCell[0], outCell[1]]

        # 掩膜sink
        # print(minX,maxX,minY,maxY)
        maskDem = dem[minX:maxX + 1, minY:maxY + 1].copy()
        maskSink = sink[minX:maxX + 1, minY:maxY + 1].copy()
        maskSink[maskSink == sinkId] = 1
        maskSink[maskSink != 1] = 0
        row, col = maskSink.shape

        # 记录出口点数据
        maskDir = dir[minX:maxX + 1, minY:maxY + 1].copy()
        maskDir[maskSink == 1] = 0
        maskDir[outCell[0] - minX, outCell[1] - minY] = outCellDir

        A = iterflat(maskSink, 0, (outCell[0] - minX, outCell[1] - minY), maskDir)  # Python
        # A = C_iterFlat(maskSink,0,outCell[0] - minX,outCell[1] - minY,row,col)   # C语言
        maskDir[outCell[0] - minX, outCell[1] - minY] = outCellDir
        # 单独存储修正你后的流向
        # temp_geo = (geo[0] + geo[1] * minY, geo[1], geo[2], geo[3] + geo[5] * minX, geo[4], geo[5])
        # outModifiedDirPath = os.path.join(outVenu, str(sinkId) + '.tif')
        # Raster.save_raster(outModifiedDirPath, A, proj, temp_geo, gdal.GDT_Byte, 0)

        # 写入原来的流向
        # C_updataDir(dir,minX,minY,A)
        for ti in range(row):
            for tj in range(col):
                if A[ti, tj] == 255:
                    continue
                dir[minX + ti, minY + tj] = A[ti, tj]

    print('Search the fake sink from the database')
    infos = query_data_byType(db_name, str(3))
    # print(infos)
    nn = 0
    for info in infos:
        # print(info)
        sinkId = info[0]

        minX = info[-5]
        maxX = info[-4]
        minY = info[-3]
        maxY = info[-2]
        outCell = (info[3], info[4])
        outCellDir = dir[outCell[0], outCell[1]]
        # print(info[0])
        # 掩膜sink
        # print(minX,maxX,minY,maxY)
        # print(sinkId)

        maskSink = sink[minX:maxX + 1, minY:maxY + 1].copy()
        maskSink[maskSink == sinkId] = 1
        maskSink[maskSink != 1] = 0  # 写成0 是因为好用：否则后面查最小值的时候会错误
        row, col = maskSink.shape
        maskDem = dem[minX:maxX + 1, minY:maxY + 1].copy()
        maskDem[maskSink == 0] = -9999
        # 记录出口点数据
        maskDir = dir[minX:maxX + 1, minY:maxY + 1].copy()  # np.zeros((row,col))
        maskDir[maskSink == 1] = 0
        maskDir[outCell[0] - minX, outCell[1] - minY] = outCellDir

        # 计算淹没内的流向
        # A = priorityFlood(maskSink, 0, maskDem, (outCell[0] - minX, outCell[1] - minY), maskDir)  # Python
        A = queuePriorityFlood(maskSink, 0, maskDem, (outCell[0] - minX, outCell[1] - minY), maskDir)  # Python
        # A = C_PriorityFlood(maskSink, 0, maskDem, outCell[0] - minX, outCell[1] - minY,row,col)  # C语言

        # 单独存储修正你后的流向
        temp_geo = (geo[0] + geo[1] * minY, geo[1], geo[2], geo[3] + geo[5] * minX, geo[4], geo[5])

        # Raster.save_raster(outModifiedDirPath,A,proj,temp_geo,gdal.GDT_Byte,0)

        # 写入原来的流向
        # C_updataDir(dir, minX, minY, A)
        for ti in range(row):
            for tj in range(col):
                if A[ti, tj] == 255:
                    continue
                dir[minX + ti, minY + tj] = A[ti, tj]


    # 另外的代码，把内流的也作为外流处理
    print('Search the endorheic sink from the database')
    infos = query_data_byType(db_name, str(2))
    # print(infos)
    nn = 0
    for info in infos:
        # print(info)
        sinkId = info[0]

        minX = info[-5]
        maxX = info[-4]
        minY = info[-3]
        maxY = info[-2]
        outCell = (info[3], info[4])
        outCellDir = dir[outCell[0], outCell[1]]
        # print(info[0])
        # 掩膜sink
        # print(minX,maxX,minY,maxY)
        # print(sinkId)

        maskSink = sink[minX:maxX + 1, minY:maxY + 1].copy()
        maskSink[maskSink == sinkId] = 1
        maskSink[maskSink != 1] = 0  # 写成0 是因为好用：否则后面查最小值的时候会错误
        row, col = maskSink.shape
        maskDem = dem[minX:maxX + 1, minY:maxY + 1].copy()
        maskDem[maskSink == 0] = -9999
        # 记录出口点数据
        maskDir = dir[minX:maxX + 1, minY:maxY + 1].copy()  # np.zeros((row,col))
        maskDir[maskSink == 1] = 0
        maskDir[outCell[0] - minX, outCell[1] - minY] = outCellDir

        # 计算淹没内的流向
        # A = priorityFlood(maskSink, 0, maskDem, (outCell[0] - minX, outCell[1] - minY), maskDir)  # Python
        A = queuePriorityFlood(maskSink, 0, maskDem, (outCell[0] - minX, outCell[1] - minY), maskDir)  # Python
        # A = C_PriorityFlood(maskSink, 0, maskDem, outCell[0] - minX, outCell[1] - minY,row,col)  # C语言

        # 单独存储修正你后的流向
        temp_geo = (geo[0] + geo[1] * minY, geo[1], geo[2], geo[3] + geo[5] * minX, geo[4], geo[5])

        # Raster.save_raster(outModifiedDirPath,A,proj,temp_geo,gdal.GDT_Byte,0)

        # 写入原来的流向
        # C_updataDir(dir, minX, minY, A)
        for ti in range(row):
            for tj in range(col):
                if A[ti, tj] == 255:
                    continue
                dir[minX + ti, minY + tj] = A[ti, tj]
    Raster.save_raster(outDir_file, dir, proj, geo, gdal.GDT_Byte, 0)
def endorheic_process(outVenu):
    """
    原main1的后半部分，生成endorheic文件夹
    :param outVenu:
    :return:
    """
    # 整体代码运行
    baseName = os.path.basename(outVenu)
    db_name = os.path.join(outVenu, baseName + "_burndb.db")
    if not os.path.exists(db_name):
        print(outVenu,'{:s} not exists'.format(db_name))
        return
    dem_file = os.path.join(outVenu, baseName + "_burnDEM.tif")
    dir_file = os.path.join(outVenu, baseName + "_burndir.tif")
    sink_file = os.path.join(outVenu, baseName + "_burnsink.tif")
    occ_file = os.path.join(outVenu, baseName + "_Occ.tif")
    outDir_file = os.path.join(outVenu, baseName + "_outDIR.tif")
    dem = Raster.get_raster(dem_file)
    dir = Raster.get_raster(dir_file)
    sink = Raster.get_raster(sink_file)
    dir = np.array(dir, dtype=np.int64)
    dir[dir == 255] = 0
    Occ = Raster.get_raster(occ_file)
    proj, geo, sinkNodata = Raster.get_proj_geo_nodata(sink_file)
    proj, geo, dirNodata = Raster.get_proj_geo_nodata(dir_file)
    proj, geo, demNodata = Raster.get_proj_geo_nodata(dem_file)
    _, _, O_Nodata = Raster.get_proj_geo_nodata(occ_file)

    ############################## 内流区特殊处理，需准备mask，fdir，dem ###########################################
    print('Search the endorheic from the database')
    infos = query_data_byType(db_name, str(2))
    # print(infos)
    con = [['record','newtype']]
    outfile = os.path.join(outVenu, baseName + "_update.csv")
    for info in infos:
        # print(info)
        sinkId = info[0]
        con.append([sinkId,3])
        minX = info[-5]
        maxX = info[-4]
        minY = info[-3]
        maxY = info[-2]
        outCell = (info[5], info[6])
        # outCellDir = dir[outCell[0], outCell[1]]
        # 掩膜sink
        maskDem = dem[minX:maxX + 1, minY:maxY + 1].copy()
        maskSink = sink[minX:maxX + 1, minY:maxY + 1].copy()
        maskSink[maskSink == sinkId] = 1
        maskSink[maskSink != 1] = 0  # 写成0 是因为好用：否则后面查最小值的时候会错误
        # maskDem[maskSink == 0] = demNodata  # 不建议打开，当洼地内部有山顶时，会导致dem缺失这部分，后续计算会出现错误
        maskOcc = Occ[minX:maxX + 1, minY:maxY + 1].copy()
        maskOcc[maskSink == 0] = O_Nodata
        row, col = maskSink.shape

        # 记录出口点数据
        maskDir = dir[minX:maxX + 1, minY:maxY + 1].copy()  # np.zeros((row,col))
        maskDir[maskSink == 0] = 255
        # maskDir[outCell[0] - minX, outCell[1] - minY] = outCellDir

        # 计算淹没内的流向
        # A = priorityFlood(maskSink, 0, maskDem, (outCell[0] - minX, outCell[1] - minY), maskDir)  # Python
        # A = C_PriorityFlood(maskSink, 0, maskDem, outCell[0] - minX, outCell[1] - minY, row, col)  # C语言
        # 单独存储修正你后的流向

        # 存储掩膜文件
        temp_geo = (geo[0] + geo[1] * minY, geo[1], geo[2], geo[3] + geo[5] * minX, geo[4], geo[5])
        endor = os.path.join(outVenu, 'endorheic')
        if not os.path.exists(endor):
            os.mkdir(endor)
            os.chmod(endor, 0o777)
        sinkIdDir = os.path.join(endor, str(sinkId))
        if not os.path.exists(sinkIdDir):
            os.mkdir(sinkIdDir)
            os.chmod(sinkIdDir, 0o777)

        outmaskSinkPath = os.path.join(sinkIdDir, str(sinkId) + '_mask_.tif')
        outmaskDemPath = os.path.join(sinkIdDir, str(sinkId) + '_dem_.tif')
        outmaskDirPath = os.path.join(sinkIdDir, str(sinkId) + '_dir_.tif')
        outmaskOccPath = os.path.join(sinkIdDir, str(sinkId) + '_occ_.tif')
        Raster.save_raster(outmaskOccPath, maskOcc, proj, temp_geo, gdal.GDT_Float32, O_Nodata)
        Raster.save_raster(outmaskSinkPath, maskSink, proj, temp_geo, gdal.GDT_Byte, 0)
        Raster.save_raster(outmaskDirPath, maskDir, proj, temp_geo, gdal.GDT_Byte, 255)
        Raster.save_raster(outmaskDemPath, maskDem, proj, temp_geo, gdal.GDT_Float32, demNodata)
        # # 存储掩膜文件


    # with open(outfile,'w',newline="") as f:
    #     writer = csv.writer(f)
    #     writer.writerows(con)
    #     f.close()
def copy(Venu,B):
    """
    将endorheic的mask.tif赋值到mask文件夹，进行目视检查
    :param Venu:
    :param B:
    :return:
    """

    import shutil
    import os
    if not os.path.exists(B):
        os.mkdir(B)
        os.chmod(B,0o777)
    venudirs = os.listdir(Venu)
    runVenus = []
    for venudir in venudirs:
        venu = os.path.join(Venu, venudir)  # 二级根目录：存放具体文件_run_datas

        baseName = os.path.basename(venu)
        runVenus.append(os.path.join(venu,baseName+"_mask_.tif"))


        shutil.copy(os.path.join(venu,baseName+"_mask_.tif"),B)



### 直接使用Priority-Flood+D8计算流向

def construct_mask_sink(dem,proj,geo,dem_nodata,out_mask_sink_file,nodata=0):
    """
    构建mask_sink
    :param dem:
    :param proj:
    :param geo:
    :param dem_nodata:
    :param nodata:
    :return:
    """
    row,col = dem.shape

    mask_sink = np.zeros((row,col),dtype = np.int8)
    mask_sink[dem != dem_nodata] = 1

    Raster.save_raster(out_mask_sink_file,mask_sink,proj,geo,gdal.GDT_Byte,nodata)
    print('mask_sink construct successfully')

def find_outlet(dem,d_nodata):
    row,col = dem.shape
    print(row,col)

    extent_cells = []
    f1 = False
    f2 = False
    seedCell = [-1,-1]
    for i in range(row):
        for j in range(col):
            if dem[i,j] == d_nodata:
                continue
            for k in range(8):
                nextCell = (i+dmove[k][0],j+dmove[k][1])
                if not check_boundary(nextCell[0],nextCell[1],row,col):
                    continue
                if dem[nextCell[0],nextCell[1]] == d_nodata:
                    seedCell = [i,j]
                    f1 = True
                    f2 = True
                    break
            if f1:
                break
        if f2:
            break
    # dfs
    print(seedCell)
    outCell = [-1,-1,math.inf]

    vis = np.zeros((row,col))
    popCells = [seedCell]

    while popCells:
        popCell = popCells.pop()
        if outCell[2] >= dem[popCell[0],popCell[1]]:
            outCell = [popCell[0],popCell[1],dem[popCell[0],popCell[1]]]

        if vis[popCell[0],popCell[1]] == 1:
            continue
        vis[popCell[0], popCell[1]] = 1
        for k in range(8):
            nextCell = (popCell[0]+dmove[k][0],popCell[1]+dmove[k][1])
            if not check_boundary(nextCell[0],nextCell[1],row,col):
                continue

            if dem[nextCell[0],nextCell[1]] == d_nodata:
                continue

            flag = False
            for kk in range(8):
                nnCell = (nextCell[0]+dmove[kk][0],nextCell[1]+dmove[kk][1])
                if not check_boundary(nnCell[0],nnCell[1],row,col):
                    flag = True
                    break
                if dem[nnCell[0],nnCell[1]] == d_nodata:
                    flag = True
                    break

            if flag:
                if vis[nextCell[0],nextCell[1]] == 0:
                    popCells.append(nextCell)


    return outCell
def priority_D8_cal_dir(venu):
    """
    直接使用Priority-Flood+D8计算流向，避免平行河网/
    要准备mask_sink，即全局的掩膜。mask_sink 为1，nodata=0
    mask_dem.
    边界上的最低点
    mask_dir、仅保留出水口的流向，其余赋值为0
    :param dem_file:
    :param burn_dir_file:
    :param mask_file:
    :return:
    """
    baseName = os.path.basename(venu)   # 不加\
    dem_file = os.path.join(venu,baseName+'_burnDEM.tif')
    burn_dir_file = os.path.join(venu,baseName+'_burndir.tif')
    mask_file = os.path.join(venu,baseName+'_mask_sink_.tif')
    out_modified_dir_file = os.path.join(venu,baseName+'_ModifiedDir.tif')

    dem = Raster.get_raster(dem_file)
    row,col = dem.shape
    proj,geo,d_nodata = Raster.get_proj_geo_nodata(dem_file)
    if not os.path.exists(burn_dir_file):
        sink.Cal_dir(dem_file,burn_dir_file)
    fdir = Raster.get_raster(burn_dir_file)

    if not os.path.exists(mask_file):
        construct_mask_sink(dem,proj,geo,d_nodata,mask_file)

    mask_sink = Raster.get_raster(mask_file)

    # 寻找出水口
    outlet_cell = find_outlet(dem,d_nodata)
    print('Find outlet successfully')
    print(outlet_cell)
    print(fdir[outlet_cell[0],outlet_cell[1]])

    mask_dir = np.zeros((row,col))
    mask_dir[outlet_cell[0],outlet_cell[1]] = fdir[outlet_cell[0],outlet_cell[1]]

    mask_dem = dem.copy()
    A = queuePriorityFlood(mask_sink,0,mask_dem,(outlet_cell[0],outlet_cell[1]),mask_dir)

    for ti in range(row):
        for tj in range(col):
            if A[ti, tj] == 255:
                continue
            fdir[ti, tj] = A[ti, tj]
    Raster.save_raster(out_modified_dir_file,fdir,proj,geo,gdal.GDT_Byte,0)
    os.chmod(out_modified_dir_file,0o777)


### 直接使用Priority-Flood+D8计算流向


def sbatch_main1(venu):
    """
    输入NorthAmerica路径
    :param venu:
    :return:
    """
    dirNames = os.listdir(venu)
    for dirName in dirNames:
        # if dirName in ['Siberia_3040017570']:
        #     continue
        # main1(dirName)
        print('************',dirName,'***********')
        if len(dirName.split('.')) > 1:
            continue
        outvenu = os.path.join(venu,dirName)
        try:
            if os.path.exists(os.path.join(outvenu,dirName+'_ModifiedDir.tif')):
                continue
            endorheic_process(outvenu)   # 生成endorheic文件夹
        except Exception as e:
            print(e)

        try:
            if os.path.exists(os.path.join(outvenu,dirName+'_ModifiedDir_.tif')):
                continue
            if os.path.exists(os.path.join(venu,dirName,'endorheic')):
                copy(os.path.join(venu,dirName,'endorheic'),os.path.join(venu,dirName,'mask'))  # 生成mask文件夹
        except Exception as e:
            print(e)

        # try:
        #     exorheic_process(outvenu)
        # except Exception as e:
        #     print(e)
def Allmain(dem_file, dir_file, sink_file, infosDic, deleteIds, outVenu):
    # 分幅处理
    # 整体代码运行

    dem = Raster.get_raster(dem_file)
    dir = Raster.get_raster(dir_file)
    sink = Raster.get_raster(sink_file)
    dir = np.array(dir, dtype=np.int64)
    dir[dir == 255] = 0
    proj, geo, sinkNodata = Raster.get_proj_geo_nodata(sink_file)
    proj, geo, dirNodata = Raster.get_proj_geo_nodata(dir_file)
    proj, geo, demNodata = Raster.get_proj_geo_nodata(dem_file)

    if 1 in infosDic:
        print('Search the flat sink from the database')
        infos = infosDic[1]
        for info in infos:
            # print(info)
            sinkId = info[0]
            outModifiedDirPath = os.path.join(outVenu, str(sinkId) + '.tif')
            if os.path.exists(outModifiedDirPath):
                continue
            minX = info[-5]
            maxX = info[-4]
            minY = info[-3]
            maxY = info[-2]
            outCell = (info[3], info[4])
            outCellDir = dir[outCell[0], outCell[1]]

            # 掩膜sink
            # print(minX,maxX,minY,maxY)
            maskDem = dem[minX:maxX + 1, minY:maxY + 1].copy()
            maskSink = sink[minX:maxX + 1, minY:maxY + 1].copy()
            maskSink[maskSink == sinkId] = 1
            maskSink[maskSink != 1] = 0
            row, col = maskSink.shape

            # 记录出口点数据
            maskDir = dir[minX:maxX + 1, minY:maxY + 1].copy()
            maskDir[maskSink == 1] = 0
            maskDir[outCell[0] - minX, outCell[1] - minY] = outCellDir

            A = iterflat(maskSink, 0, (outCell[0] - minX, outCell[1] - minY), maskDir)  # Python
            # A = C_iterFlat(maskSink,0,outCell[0] - minX,outCell[1] - minY,row,col)   # C语言
            maskDir[outCell[0] - minX, outCell[1] - minY] = outCellDir
            # 单独存储修正你后的流向
            # temp_geo = (geo[0] + geo[1] * minY, geo[1], geo[2], geo[3] + geo[5] * minX, geo[4], geo[5])
            # outModifiedDirPath = os.path.join(outVenu, str(sinkId) + '.tif')
            # Raster.save_raster(outModifiedDirPath, A, proj, temp_geo, gdal.GDT_Byte, 0)

            # 写入原来的流向
            # C_updataDir(dir,minX,minY,A)
            for ti in range(row):
                for tj in range(col):
                    if A[ti, tj] == 0:
                        continue
                    dir[minX + ti, minY + tj] = A[ti, tj]

    if 3 in infosDic:
        print('Search the fake sink from the database')
        infos = infosDic[3]
        # print(infos)
        nn = 0
        for info in infos:

            # print(info)
            sinkId = info[0]
            if nn == 1:
                print('Fake sinkid', sinkId)
                nn += 1
                continue
            if sinkId == 5883:
                nn += 1
            outModifiedDirPath = os.path.join(outVenu, str(sinkId) + '.tif')
            if os.path.exists(outModifiedDirPath):
                continue
            minX = info[-5]
            maxX = info[-4]
            minY = info[-3]
            maxY = info[-2]
            outCell = (info[3], info[4])
            outCellDir = dir[outCell[0], outCell[1]]

            # 掩膜sink
            # print(minX,maxX,minY,maxY)
            # print(sinkId)
            maskDem = dem[minX:maxX + 1, minY:maxY + 1].copy()
            maskSink = sink[minX:maxX + 1, minY:maxY + 1].copy()
            maskSink[maskSink == sinkId] = 1
            maskSink[maskSink != 1] = 0  # 写成0 是因为好用：否则后面查最小值的时候会错误
            row, col = maskSink.shape

            # 记录出口点数据
            maskDir = dir[minX:maxX + 1, minY:maxY + 1].copy()  # np.zeros((row,col))

            maskDir[maskSink == 1] = 0
            maskDir[outCell[0] - minX, outCell[1] - minY] = outCellDir

            # 计算淹没内的流向
            A = priorityFlood(maskSink, 0, maskDem, (outCell[0] - minX, outCell[1] - minY), maskDir)  # Python
            # A = C_PriorityFlood(maskSink, 0, maskDem, outCell[0] - minX, outCell[1] - minY,row,col)  # C语言

            # 单独存储修正你后的流向
            temp_geo = (geo[0] + geo[1] * minY, geo[1], geo[2], geo[3] + geo[5] * minX, geo[4], geo[5])

            # Raster.save_raster(outModifiedDirPath,A,proj,temp_geo,gdal.GDT_Byte,0)

            # 写入原来的流向
            # C_updataDir(dir, minX, minY, A)
            for ti in range(row):
                for tj in range(col):
                    if A[ti, tj] == 0:
                        continue
                    dir[minX + ti, minY + tj] = A[ti, tj]

    if 2 in infosDic:
        print('Search the endorheic from the database')
        infos = infosDic[2]
        # print(infos)
        for info in infos:
            # print(info)
            sinkId = info[0]
            minX = info[-5]
            maxX = info[-4]
            minY = info[-3]
            maxY = info[-2]
            outCell = (info[5], info[6])
            # outCellDir = dir[outCell[0], outCell[1]]

            # 掩膜sink
            maskDem = dem[minX:maxX + 1, minY:maxY + 1].copy()
            maskSink = sink[minX:maxX + 1, minY:maxY + 1].copy()
            maskSink[maskSink == sinkId] = 1
            maskSink[maskSink != 1] = 0  # 写成0 是因为好用：否则后面查最小值的时候会错误
            row, col = maskSink.shape

            # 记录出口点数据
            maskDir = dir[minX:maxX + 1, minY:maxY + 1].copy()  # np.zeros((row,col))

            maskDir[maskSink == 1] = 0
            # maskDir[outCell[0] - minX, outCell[1] - minY] = outCellDir

            # 计算淹没内的流向
            A = priorityFlood(maskSink, 0, maskDem, (outCell[0] - minX, outCell[1] - minY), maskDir)  # Python
            # A = C_PriorityFlood(maskSink, 0, maskDem, outCell[0] - minX, outCell[1] - minY, row, col)  # C语言
            # 单独存储修正你后的流向
            temp_geo = (geo[0] + geo[1] * minY, geo[1], geo[2], geo[3] + geo[5] * minX, geo[4], geo[5])

            # 写入原来的流向
            # C_updataDir(dir, minX, minY, A)
            for ti in range(row):
                for tj in range(col):
                    if A[ti, tj] == 0:
                        continue
                    dir[minX + ti, minY + tj] = A[ti, tj]

    # 处理边界污染的sink，使其流向赋为空值
    for deleteId in deleteIds:
        dir[dir == deleteId] = 0
    row, col = dir.shape
    dir[0, :] = 0
    dir[row - 1, :] = 0
    dir[0, :] = 0
    dir[col - 1, :] = 0
    demFile = os.path.basename(dem_file)
    outModifiedDirPath = os.path.join(outVenu, 'Modify_' + demFile)
    # outSink = r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\20240923\Sink2.tif'
    Raster.save_raster(outModifiedDirPath, dir, proj, geo, gdal.GDT_Byte, 0)
    # Raster.save_raster(outSink, sink, proj, geo, gdal.GDT_Float32, sinkNodata)




if __name__ == '__main__':
    # # PF示例数据
    # mask = np.array(
    #     [[0, 0, 0, 1, 1, 1], [0, 1, 1, 1, 1, 0], [0, 1, 1, 0, 0, 0], [1, 1, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0]])
    # newDir = np.array(
    #     [[4, 4, 1, 0, 0, 0], [1, 0, 0, 0, 0, 64], [1, 0, 0, 64, 32, 16], [0, 0, 16, 32, 64, 16], [4, 8, 8, 16, 64, 32]])
    # # mask = np.array([[0,0,1,0],[1,1,1,1],[0,1,1,0],[0,1,1,0]])
    #
    # dem = np.array(
    #     [[14, 62, 59, 2, 3, 4],
    #      [88, 5, 4, 4, 5, 31],
    #      [6, 5, 5, 6, 6, 74],
    #      [5, 5, 6, 7, 22,49],
    #      [4, 5, 6, 22, 19,27]])
    # # newDir = np.array([[0,0,64,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
    # A = priorityFlood(mask, 0, dem, (4,1), newDir)
    # print(A)

    # process2示例数据
    # mask = np.array(
    #     [[1, 1, 1, 1, 1], [1, 1, 1, 1, 0], [1, 1, 1, 1, 1], [1, 1, 1, 1, 0], [0, 1, 0, 0, 0]])
    # newDir = np.array(
    #     [[4, 4, 1, 0, 0], [1, 0, 0, 0, 64], [1, 0, 0, 64, 32], [0, 0, 16, 32, 64], [4, 8, 8, 16, 64]])
    # # mask = np.array([[0,0,1,0],[1,1,1,1],[0,1,1,0],[0,1,1,0]])
    #
    # dem = np.array(
    #     [[6,6,6,7,7],
    #      [9,5,8,6,9],
    #      [10,4,3,4,7],
    #      [10,7,10,7,9],
    #      [1,8,1,1,1]])
    # A = priorityFlood(mask, 0, dem, (2,2), newDir)
    # print(A)

    # main
    # dem_file = r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\DEM.tif'
    # dir_file = r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\Dir1.tif'
    # sink_file = r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\sink1.tif' # r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\sink1.tif'
    # A = r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\20240923\D2.tif'
    # outVenu = r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\outDir'
    # t1 = time.time()
    # # process3(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\sinkdb.db', dem_file, dir_file, sink_file,A)
    # # process1(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\sinkdb.db', dem_file, dir_file, sink_file,
    # #          A)
    # main(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\sinkdb.db', dem_file, dir_file, sink_file,
    #          A,outVenu)
    # t2 = time.time()
    # print('Time consuming {:2} h'.format((t2 - t1) / 60 / 60))

    # tbp
    # dem_file = "/datanode05/zhangbin/TBP_Stream/TBP/TBP_FABDEM.tif"
    # dir_file = "/datanode05/zhangbin/TBP_Stream/TBP/TBP_FABDEM_Dir.tif"
    # sink_file = "/datanode05/zhangbin/TBP_Stream/DATA/20240818/TBP_sink.tif"
    # A = "/datanode05/zhangbin/TBP_Stream/DATA/20240818/TBP_Dir_Modified.tif"
    # sinkdb = "/datanode05/zhangbin/TBP_Stream/DATA/20240818/sinkdb.db"
    # outVenu = "/datanode05/zhangbin/TBP_Stream/DATA/20240818/modifiedDir"
    # t1 = time.time()
    # # process3(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\sinkdb.db', dem_file, dir_file, sink_file,A)
    # # process1(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\rundata_20240604\sinkdb.db', dem_file, dir_file, sink_file,
    # #          A)
    # main(sinkdb, dem_file, dir_file, sink_file,
    #          A,outVenu)
    # t2 = time.time()
    # print('Time consuming {:2} h'.format((t2 - t1) / 60 / 60))

    # maskvenu = '/datanode05/zhangbin/TBP_Stream/DATA/sinkMask'
    # demvenu = '/datanode05/zhangbin/TBP_Stream/DATA/sinkDEM'
    # dirvenu = '/datanode05/zhangbin/TBP_Stream/DATA/sinkDir'
    # outDirPath = '/datanode05/zhangbin/TBP_Stream/DATA/20240818/modifiedDir'
    # db_name = '/datanode05/zhangbin/TBP_Stream/DATA/20240818/sinkdb.db'
    # main(maskvenu,demvenu,dirvenu,outDirPath,db_name)

    # maskfile = r'F:\青藏高原水体数据集\DATA\TBP_endorheic\Mask_7758_.tif'
    # demfile = r'F:\青藏高原水体数据集\DATA\TBP_endorheic\Dem7758_.tif'
    # dirfile = r'F:\青藏高原水体数据集\DATA\TBP_endorheic\Dir_7758_.tif'
    # mask = Raster.get_raster(maskfile)
    # proj,geo,masknodata = Raster.get_proj_geo_nodata(maskfile)
    # dem = Raster.get_raster(demfile)
    # dir = Raster.get_raster(dirfile)
    # A = priorityFlood(mask,masknodata,dem,(2637, 675),dir)
    # Raster.save_raster(r'F:\青藏高原水体数据集\DATA\TBP_endorheic\NewDir_7758_.tif',A,proj,geo,gdal.GDT_Byte,255)

    maskfile = r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\globalRegion\C_trail\sink_N27E085_FABDEM_V1-0.tif'
    demfile = r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\globalRegion\C_trail\N27E085_FABDEM_V1-0.tif'
    dirfile = r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\globalRegion\C_trail\dir_N27E085_FABDEM_V1-0.tif'
    # mask = Raster.get_raster(maskfile)
    # proj, geo, masknodata = Raster.get_proj_geo_nodata(maskfile)
    # dem = Raster.get_raster(demfile)
    # dir = Raster.get_raster(dirfile)
    # A = priorityFlood(mask, masknodata, dem, (2637, 675), dir)
    # Raster.save_raster(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\globalRegion\trail\NewDir.tif', A, proj, geo, gdal.GDT_Byte, 255)
    # main1(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\globalRegion\C_trail\sinkdb.db',
    #       demfile, dirfile, maskfile,
    #       r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\globalRegion\C_trail\NewDir.tif',
    #       r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\globalRegion\C_trail')

    # pool_main(r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\globalRegion\C_trail\sinkdb.db',
    #           demfile,dirfile,maskfile,
    #           r'F:\青藏高原水体数据集\DATA\青藏高原水体研究\globalRegion\C_trail')

    pass
