#!/usr/bin/env python
#-*- coding:utf-8 -*-
# datetime:2021/1/23 22:00
import os

import pandas as pd
import numpy as np
# 忽略警告
import warnings
warnings.filterwarnings("ignore")
import re

AminoAcids = ["A", "C", "D", "E", "F", "G","H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
AminoAcids2 = []
for i in range(len(AminoAcids)):
    for j in range(len(AminoAcids)):
        AminoAcids2.append(AminoAcids[i]+AminoAcids[j])


class A_A:
    def __init__(self):
        train_data_pos = pd.read_table('../Sample/train/pos', sep=',', names=['Amino', 'Label'])
        test_data_pos = pd.read_table('../Sample/test/pos', sep=',', names=['Amino', 'Label'])

        self.amino_pos = pd.concat([train_data_pos, test_data_pos], axis=0, ignore_index=True)

    # 获取A_A矩阵 n取0,1,2,3,4,5,6
    def matrix(self, n):
        M = [[0]*20 for _ in range(20)]
        amino = self.amino_pos

        for i in range(20):
            for j in range(i, 20):
                p1 = AminoAcids[i]+'.'*n+AminoAcids[j]
                p2 = AminoAcids[j]+'.'*n+AminoAcids[i]
                for k in range(len(amino)):
                    if len(re.findall(p1, amino['Amino'][k])):
                        M[i][j] += 1
                    if i != j and len(re.findall(p2, amino['Amino'][k])):
                        M[j][i] += 1
        return M

    def diff(self, n):
        M = self.matrix(n)
        pattern = []
        rowSum = []
        for i in range(20):
            s = 0
            for j in range(20):
                s += M[i][j]
            rowSum.append(s)
        colSum = []
        for j in range(20):
            s = 0
            for i in range(20):
                s += M[i][j]
            colSum.append(s)
        Sum = sum(rowSum)
        for i in range(20):
            for j in range(20):
                # if M[i][j] < 5:
                #     continue
                p_ij = M[i][j]/Sum
                p_i_ = rowSum[i]/Sum
                p__j = colSum[j]/Sum
                diff = (p_ij - (p_i_ * p__j)) / (np.sqrt((p_i_ * p__j * (1 - p_i_) * (1 - p__j)) / Sum))
                if diff >= 1.96:
                    w = self.weight(p_ij, p_i_, p__j)
                    if w > 0:
                        pattern.append([AminoAcids[i]+'.'*n+AminoAcids[j], w, M[i][j]])
        return pattern

    def weight(self, p_ij, p_i_, p__j):
        if p_ij == p_i_:
            return 0
        w = np.log(p_ij / (p_i_ * p__j)) - np.log((p_i_ - p_ij) / (p_i_ * (1 - p__j)))
        return round(w, 2)

    def run(self):
        pattern = []
        for n in range(7):
            pattern.extend(self.diff(n))
        return pattern


class A_AB:
    def __init__(self):
        train_data_pos = pd.read_table('../Sample/train/pos', sep=',', names=['Amino', 'Label'])
        test_data_pos = pd.read_table('../Sample/test/pos', sep=',', names=['Amino', 'Label'])

        self.amino_pos = pd.concat([train_data_pos, test_data_pos], axis=0, ignore_index=True)

    # 获取A_AB矩阵 n取0,1,2,3,4,5
    def matrix(self, n):  # 20*400
        M = [[0]*400 for _ in range(20)]
        amino = self.amino_pos

        for i in range(20):
            for j in range(400):
                p = AminoAcids[i]+'.'*n+AminoAcids2[j]
                for k in range(len(amino)):
                    if len(re.findall(p, amino['Amino'][k])):
                        M[i][j] += 1
        return M

    def diff(self, n):
        M = self.matrix(n)
        pattern = []
        rowSum = []
        for i in range(20):
            s = 0
            for j in range(400):
                s += M[i][j]
            rowSum.append(s)
        colSum = []
        for j in range(400):
            s = 0
            for i in range(20):
                s += M[i][j]
            colSum.append(s)
        Sum = sum(rowSum)
        for i in range(20):
            for j in range(400):
                # if M[i][j] < 5:
                #     continue
                p_ij = M[i][j]/Sum
                p_i_ = rowSum[i]/Sum
                p__j = colSum[j]/Sum
                diff = (p_ij - (p_i_ * p__j)) / (np.sqrt((p_i_ * p__j * (1 - p_i_) * (1 - p__j)) / Sum))
                if diff >= 1.96:
                    w = self.weight(p_ij, p_i_, p__j)
                    if w > 0:
                        pattern.append([AminoAcids[i]+'.'*n+AminoAcids2[j], w, M[i][j]])
        return pattern

    def weight(self, p_ij, p_i_, p__j):
        if p_ij == p_i_:
            return 0
        w = np.log(p_ij / (p_i_ * p__j)) - np.log((p_i_ - p_ij) / (p_i_ * (1 - p__j)))
        return round(w, 2)

    def run(self):
        pattern = []
        for n in range(6):
            pattern.extend(self.diff(n))
        return pattern


class AB_A:
    def __init__(self):
        train_data_pos = pd.read_table('../Sample/train/pos', sep=',', names=['Amino', 'Label'])
        test_data_pos = pd.read_table('../Sample/test/pos', sep=',', names=['Amino', 'Label'])

        self.amino_pos = pd.concat([train_data_pos, test_data_pos], axis=0, ignore_index=True)

    # 获取A_AB矩阵 n取0,1,2,3,4,5
    def matrix(self, n):  # 400*20
        M = [[0]*20 for _ in range(400)]
        amino = self.amino_pos

        for i in range(400):
            for j in range(20):
                p = AminoAcids2[i]+'.'*n+AminoAcids[j]
                for k in range(len(amino)):
                    if len(re.findall(p, amino['Amino'][k])):
                        M[i][j] += 1
        return M

    def diff(self, n):
        M = self.matrix(n)
        pattern = []
        rowSum = []
        for i in range(400):
            s = 0
            for j in range(20):
                s += M[i][j]
            rowSum.append(s)
        colSum = []
        for j in range(20):
            s = 0
            for i in range(400):
                s += M[i][j]
            colSum.append(s)
        Sum = sum(rowSum)
        for i in range(400):
            for j in range(20):
                # if M[i][j] < 5:
                #     continue
                p_ij = M[i][j]/Sum
                p_i_ = rowSum[i]/Sum
                p__j = colSum[j]/Sum
                diff = (p_ij - (p_i_ * p__j)) / (np.sqrt((p_i_ * p__j * (1 - p_i_) * (1 - p__j)) / Sum))
                if diff >= 1.96:
                    w = self.weight(p_ij, p_i_, p__j)
                    if w > 0:
                        pattern.append([AminoAcids2[i]+'.'*n+AminoAcids[j], w, M[i][j]])
        return pattern

    def weight(self, p_ij, p_i_, p__j):
        if p_ij == p_i_:
            return 0
        w = np.log(p_ij / (p_i_ * p__j)) - np.log((p_i_ - p_ij) / (p_i_ * (1 - p__j)))
        return round(w, 2)

    def run(self):
        pattern = []
        for n in range(6):
            pattern.extend(self.diff(n))
        return pattern


if __name__ == '__main__':
    print("Start extract feature")
    if not os.path.exists("../Evocleave V2.0"):
        os.makedirs("../Evocleave V2.0")
    pattern_A_A = A_A().run()
    dc = pd.DataFrame(data=pattern_A_A, columns=["AminoT2", "WeightC", "Count"])
    dc.to_csv('../Evocleave V2.0/A_A-dc.csv', index=False)

    pattern_A_AB = A_AB().run()
    dc = pd.DataFrame(data=pattern_A_AB, columns=["AminoT2", "WeightC", "Count"])
    dc.to_csv('../Evocleave V2.0/A_AB-dc.csv', index=False)

    pattern_AB_A = AB_A().run()
    dc = pd.DataFrame(data=pattern_AB_A, columns=["AminoT2", "WeightC", "Count"])
    dc.to_csv('../Evocleave V2.0/AB_A-dc.csv', index=False)

    print("Extract feature is Saved")