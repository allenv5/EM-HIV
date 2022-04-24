#!/usr/bin/env python
# -*- coding:utf-8 -*-
# datetime:2020/12/8 10:49

import os

import numpy as np
import pandas as pd
import re
import warnings

from sklearn import svm

warnings.filterwarnings("ignore")

# 二十种不同的氨基酸
AminoAcids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
ChemicalFeatures = [("C", "M"), ("A", "G", "P"), ("I", "L", "V"), ("D", "E"),
                    ("H", "K", "R"), ("F", "W", "Y"), ("N", "Q"), ("S", "T")]
PosLabel = '1'  # 可切割标记为1
NegLabel = '-1'  # 不可切割标记为-1


def __feature_oc(dc, amino_data):
    print("Use orthogonal coding scheme to construct feature vector")
    orthogonal_features = dict()
    for x, amino in enumerate(AminoAcids[:-1]):
        num = [-1] * 19
        num[x] = 1
        orthogonal_features[amino] = num[:]
    orthogonal_features["Y"] = [-1] * 19
    train_x_data = []
    train_y_data = []
    train = []
    for index in range(len(amino_data)):
        feature = []

        for x in amino_data["Amino"][index]:
            feature.extend(orthogonal_features[x])

        label = amino_data['Label'][index]
        train_x_data.append(feature)
        train_y_data.append(label)
    train.append(train_x_data)
    train.append(train_y_data)
    return train


def __feature_cc(dc, amino_data):
    print("Use chemical coding schemes to construct feature vectors")
    chemical_features = dict()
    for x, amino in enumerate(ChemicalFeatures[:-1]):
        num = [-1] * 7
        num[x] = 1
        chemical_features[amino] = num[:]
    chemical_features[ChemicalFeatures[-1]] = [-1] * 7

    train_x_data = []
    train_y_data = []
    train = []
    for index in range(len(amino_data)):
        feature = []

        for x in amino_data["Amino"][index]:
            for group in ChemicalFeatures:
                if x in group:
                    feature.extend(chemical_features[group])

        label = amino_data['Label'][index]
        train_x_data.append(feature)
        train_y_data.append(label)
    train.append(train_x_data)
    train.append(train_y_data)
    return train


def __feature_dc(dc, amino_data):
    print("Use Evocleave V2.0 schemes to construct feature vectors")
    train_x_data = []
    train_y_data = []
    train = []
    for index in range(len(amino_data)):
        feature = []

        for j in range(len(dc)):
            if re.search(dc['AminoT2'][j], amino_data['Amino'][index]):
                feature.append(dc['WeightC'][j])
            else:
                feature.append(0)
        label = amino_data['Label'][index]
        train_x_data.append(feature)
        train_y_data.append(label)
    train.append(train_x_data)
    train.append(train_y_data)
    return train


def __feature_oc_cc(dc, amino_data):
    print("Use orthogonal coding schemes and chemical coding schemes to construct feature vectors")
    orthogonal_features = dict()
    for x, amino in enumerate(AminoAcids[:-1]):
        num = [-1] * 19
        num[x] = 1
        orthogonal_features[amino] = num[:]
    orthogonal_features["Y"] = [-1] * 19

    chemical_features = dict()
    for x, amino in enumerate(ChemicalFeatures[:-1]):
        num = [-1] * 7
        num[x] = 1
        chemical_features[amino] = num[:]
    chemical_features[ChemicalFeatures[-1]] = [-1] * 7

    train_x_data = []
    train_y_data = []
    train = []
    for index in range(len(amino_data)):
        feature = []

        for x in amino_data["Amino"][index]:
            feature.extend(orthogonal_features[x])

        for x in amino_data["Amino"][index]:
            for group in ChemicalFeatures:
                if x in group:
                    feature.extend(chemical_features[group])

        label = amino_data['Label'][index]
        train_x_data.append(feature)
        train_y_data.append(label)
    train.append(train_x_data)
    train.append(train_y_data)
    return train


def __feature_cc_dc(dc, amino_data):
    print("Use chemical coding schemes and Evocleave coding schemes to construct feature vectors")
    chemical_features = dict()
    for x, amino in enumerate(ChemicalFeatures[:-1]):
        num = [-1] * 7
        num[x] = 1
        chemical_features[amino] = num[:]
    chemical_features[ChemicalFeatures[-1]] = [-1] * 7

    train_x_data = []
    train_y_data = []
    train = []
    for index in range(len(amino_data)):
        feature = []

        for x in amino_data["Amino"][index]:
            for group in ChemicalFeatures:
                if x in group:
                    feature.extend(chemical_features[group])

        for j in range(len(dc)):
            if re.search(dc['AminoT2'][j], amino_data['Amino'][index]):
                feature.append(dc['WeightC'][j])
            else:
                feature.append(0)

        label = amino_data['Label'][index]
        train_x_data.append(feature)
        train_y_data.append(label)

    train.append(train_x_data)
    train.append(train_y_data)
    return train


def __feature_oc_dc(dc, amino_data):
    print("Use orthogonal coding schemes and Evocleave coding schemes to construct feature vectors")
    orthogonal_features = dict()
    for x, amino in enumerate(AminoAcids[:-1]):
        num = [-1] * 19
        num[x] = 1
        orthogonal_features[amino] = num[:]
    orthogonal_features["Y"] = [-1] * 19

    train_x_data = []
    train_y_data = []
    train = []
    for index in range(len(amino_data)):
        feature = []

        for x in amino_data["Amino"][index]:
            feature.extend(orthogonal_features[x])

        for j in range(len(dc)):
            if re.search(dc['AminoT2'][j], amino_data['Amino'][index]):
                feature.append(dc['WeightC'][j])
            else:
                feature.append(0)

        label = amino_data['Label'][index]
        train_x_data.append(feature)
        train_y_data.append(label)

    train.append(train_x_data)
    train.append(train_y_data)
    return train


def __feature_mul(dc, amino_data):
    print("Use orthogonal coding schemes, chemical coding schemes "
          "and Evocleave V2.0 coding schemes to construct feature vectors")
    orthogonal_features = dict()
    for x, amino in enumerate(AminoAcids[:-1]):
        num = [-1] * 19
        num[x] = 1
        orthogonal_features[amino] = num[:]
    orthogonal_features["Y"] = [-1] * 19

    chemical_features = dict()
    for x, amino in enumerate(ChemicalFeatures[:-1]):
        num = [-1] * 7
        num[x] = 1
        chemical_features[amino] = num[:]
    chemical_features[ChemicalFeatures[-1]] = [-1] * 7

    train_x_data = []
    train_y_data = []
    train = []
    for index in range(len(amino_data)):
        feature = []

        for x in amino_data["Amino"][index]:
            feature.extend(orthogonal_features[x])

        for x in amino_data["Amino"][index]:
            for group in ChemicalFeatures:
                if x in group:
                    feature.extend(chemical_features[group])

        feature_evo = []
        for j in range(len(dc)):
            if re.search(dc['AminoT2'][j], amino_data['Amino'][index]):
                feature_evo.append(dc['WeightC'][j])
            else:
                feature_evo.append(0)
        feature.extend(feature_evo)

        label = amino_data['Label'][index]
        train_x_data.append(feature)
        train_y_data.append(label)
    train.append(train_x_data)
    train.append(train_y_data)
    return train


def __extract_dc():
    ls = []
    for tp in ("A_A", "A_AB", "AB_A"):
        try:
            ls.append(pd.read_csv('../Evocleave V2.0/{}-dc.csv'.format(tp)))
        except Exception as e:
            print(e)
            exit()
    dc = pd.concat(ls, axis=0, ignore_index=True)
    dc = dc[dc["Count"] >= 5]
    dc.index = range(len(dc))
    return dc


def __train_test_model(feature_method, t, c1, beta, dc, input_file, train_pos_seq, train_neg_seq, test_seq, cv=-1):
    print("AsBaggingSVM is running")
    try:
        if cv == -1:
            if not os.path.exists("{}/results".format(input_file)):
                os.makedirs("{}/results".format(input_file))
            f = open('{}/results/result.txt'.format(input_file), 'w+', encoding='utf-8')
        else:
            if not os.path.exists("{}/{}/results".format(input_file, str(cv))):
                os.makedirs("{}/{}/results".format(input_file, str(cv)))
            f = open('{}/{}/results/result.txt'.format(input_file, str(cv)), 'w+', encoding='utf-8')
    except Exception as e:
        print(e)
        exit()

    test_data = feature_method(dc, test_seq)
    x_test, y_test = test_data
    y_test = list(y_test)

    clf = svm.SVC(kernel='linear', probability=True, random_state=42, max_iter=2000,
                  class_weight={-1: c1 / beta, 1: c1})

    scores = []

    for _ in range(t):
        # sampling with replacement
        train_neg_seq_s = train_neg_seq.sample(len(train_pos_seq), replace=True)
        train_neg_seq_s.index = range(len(train_neg_seq_s))

        train_seq = pd.concat([train_pos_seq, train_neg_seq_s], axis=0, ignore_index=True)

        train_data = feature_method(dc, train_seq)

        x_train, y_train = train_data

        clf.fit(x_train, y_train)
        test_score = []
        for score in clf.predict_proba(x_test):
            test_score.append(round(score[1], 5))
        scores.append(test_score)

    s = list(np.mean(np.array(scores), axis=0))

    f.write("Amino label 1\n")
    for inx in range(len(s)):
        f.write("{} {} {}\n".format(test_seq["Amino"][inx], y_test[inx], s[inx]))
    f.close()
    print("The result of AsBaggingSVM method is saved")


def run(feature, t, c1, beta, input_file, cv):
    code_feature = {
        "0": __feature_oc,
        "1": __feature_cc,
        "2": __feature_dc,
        "3": __feature_oc_cc,
        "4": __feature_cc_dc,
        "5": __feature_oc_dc,
        "6": __feature_mul
    }
    t, c1, beta, cv = eval(t), eval(c1), eval(beta), eval(cv)
    if c1 <= 0 or beta <= 0:
        print(" c1>=0 and beta>0")
        exit()
    if feature not in code_feature.keys():
        print("the code of feature selection is error")
        print("code 0: feature_oc")
        print("code 1: feature_cc")
        print("code 2: feature_dc")
        print("code 3: feature_oc_cc")
        print("code 4: feature_cc_dc")
        print("code 5: feature_oc_dc")
        print("code 6: feature_mul")
        return
    else:
        feature_method = code_feature[feature]
        dc = __extract_dc()
        if cv == -1:
            try:
                frame_train_p = pd.read_table('{}/train/pos'.format(input_file), sep=',',
                                              names=['Amino', 'Label'])
                frame_train_n = pd.read_table('{}/train/neg'.format(input_file), sep=',',
                                              names=['Amino', 'Label'])
                frame_train = pd.concat([frame_train_p, frame_train_n], axis=0)
                frame_train.index = range(len(frame_train))

                frame_test_p = pd.read_table('{}/test/pos'.format(input_file), sep=',',
                                             names=['Amino', 'Label'])
                frame_test_n = pd.read_table('{}/test/neg'.format(input_file), sep=',',
                                             names=['Amino', 'Label'])
                frame_test = pd.concat([frame_test_p, frame_test_n], axis=0)
                frame_test.index = range(len(frame_test))

                __train_test_model(feature_method, t, c1, beta, dc, input_file,
                                   frame_train_p, frame_train_n, frame_test)
            except Exception as e:
                print(e)
                exit()
        else:
            for i in range(1, cv + 1):
                try:
                    frame_train_p = pd.read_table('{}/{}/train/pos'.format(input_file, str(i)), sep=',',
                                                  names=['Amino', 'Label'])
                    frame_train_n = pd.read_table('{}/{}/train/neg'.format(input_file, str(i)), sep=',',
                                                  names=['Amino', 'Label'])
                    frame_train = pd.concat([frame_train_p, frame_train_n], axis=0)
                    frame_train.index = range(len(frame_train))

                    frame_test_p = pd.read_table('{}/{}/test/pos'.format(input_file, str(i)), sep=',',
                                                 names=['Amino', 'Label'])
                    frame_test_n = pd.read_table('{}/{}/test/neg'.format(input_file, str(i)), sep=',',
                                                 names=['Amino', 'Label'])
                    frame_test = pd.concat([frame_test_p, frame_test_n], axis=0)
                    frame_test.index = range(len(frame_test))

                    __train_test_model(feature_method, t, c1, beta, dc, input_file,
                                       frame_train_p, frame_train_n, frame_test, i)
                except Exception as e:
                    print(e)
                    exit()