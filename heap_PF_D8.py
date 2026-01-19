# -*- coding: utf-8 -*-
"""
@Time ： 2025/11/17 22:27
@Auth ：
@File ：heap_PF_D8.py
@IDE ：PyCharm
"""
import heapq
from collections import deque
import numpy as np

from genral_functions import check_boundary, reverseD8, d8

dmove=[(0,1),(1,1),(1,0),(1,-1),(0,-1),(-1,-1),(-1,0),(-1,1)]
dmove_dic = {1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1), 16: (0, -1), 32: (-1, -1), 64: (-1, 0), 128: (-1, 1)}

def optimized_flow_repair(dem, mask,maskNodata, startcell):
    """

    :param dem:
    :param mask:
    :param maskNodata:
    :param startcell:  ！！！必须一定是tuple形式！！！！！
    :return:
    """

    row, col = dem.shape
    Vis = np.zeros((row, col), dtype=np.uint8)
    newDir = np.zeros((row, col), dtype=np.int16)

    # 三队列
    quickQ = deque()
    maxQ = []      # ( -height, r, c ) 最大值优先
    minQ = []      # (  height, r, c ) 最小值优先

    # 初始化
    Vis[startcell] = 1
    startH = dem[startcell]
    heapq.heappush(maxQ, (-startH, startcell[0], startcell[1]))

    # 主循环：无全矩阵扫描 🚀
    while quickQ or maxQ or minQ:

        # ------------------------------------------------
        # 1. 平地优先（和你原逻辑一致）
        # ------------------------------------------------
        while quickQ:
            r, c = quickQ.popleft()
            baseH = dem[r, c]

            for k in range(8):
                nr = r + dmove[k][0]
                nc = c + dmove[k][1]
                if not check_boundary(nr, nc, row, col):
                    continue
                if mask[nr, nc] == maskNodata:
                    continue
                if Vis[nr, nc]:
                    continue

                nh = dem[nr, nc]

                if nh == baseH:  # 平地
                    newDir[nr, nc] = 2 ** ((k + 4) % 8)
                    Vis[nr, nc] = 1
                    quickQ.append((nr, nc))
                else:
                    if nh < baseH:     # 下坡
                        newDir[nr, nc] = reverseD8((nr, nc), dem, Vis,mask)
                        heapq.heappush(maxQ, (-nh, nr, nc))
                    else:               # 上坡
                        newDir[nr, nc] = d8((nr, nc), dem, Vis,mask)
                        heapq.heappush(minQ, (nh, nr, nc))
                    Vis[nr, nc] = 1

        # ------------------------------------------------
        # 2. 下坡优先（取最大值）
        # ------------------------------------------------
        if maxQ:
            _, r, c = heapq.heappop(maxQ)
            baseH = dem[r, c]

            for k in range(8):
                nr = r + dmove[k][0]
                nc = c + dmove[k][1]
                if not check_boundary(nr, nc, row, col):
                    continue
                if mask[nr, nc] == maskNodata:
                    continue
                if Vis[nr, nc]:
                    continue

                nh = dem[nr, nc]

                if nh == baseH:
                    newDir[nr, nc] = 2 ** ((k + 4) % 8)
                    Vis[nr, nc] = 1
                    quickQ.append((nr, nc))
                else:
                    if nh < baseH:
                        newDir[nr, nc] = reverseD8((nr, nc), dem, Vis,mask)
                        heapq.heappush(maxQ, (-nh, nr, nc))
                    else:
                        newDir[nr, nc] = d8((nr, nc), dem, Vis,mask)
                        heapq.heappush(minQ, (nh, nr, nc))
                    Vis[nr, nc] = 1
            continue  # 与原逻辑一致（每次只处理一个下坡中心）

        # ------------------------------------------------
        # 3. 上坡（取最小值）
        # ------------------------------------------------
        if minQ:
            baseH, r, c = heapq.heappop(minQ)

            for k in range(8):
                nr = r + dmove[k][0]
                nc = c + dmove[k][1]
                if not check_boundary(nr, nc, row, col):
                    continue
                if mask[nr, nc] == maskNodata:
                    continue
                if Vis[nr, nc]:
                    continue

                nh = dem[nr, nc]

                if nh == baseH:
                    newDir[nr, nc] = 2 ** ((k + 4) % 8)
                    Vis[nr, nc] = 1
                    quickQ.append((nr, nc))
                else:
                    if nh < baseH:
                        newDir[nr, nc] = reverseD8((nr, nc), dem, Vis,mask)
                        heapq.heappush(maxQ, (-nh, nr, nc))
                    else:
                        newDir[nr, nc] = d8((nr, nc), dem, Vis,mask)
                        heapq.heappush(minQ, (nh, nr, nc))
                    Vis[nr, nc] = 1

    return newDir
