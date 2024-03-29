#!/usr/bin/env python
# -*- coding:utf-8 -*-
# datetime:2022/4/22 15:51


import argparse
import AsBaggingSVM as em

parser = argparse.ArgumentParser(description="The AsBaggingSVM method")
parser.add_argument("--feature", "-f", help="feature_code", default="6")
parser.add_argument("--bootstraps", "-t", help="number of bootstraps", default="30")
parser.add_argument("--c1", "-c1", help="", default="8")
parser.add_argument("--beta", "-beta", help="", default="2")
parser.add_argument("--input", "-i", help="input folder", default="../Sample")
parser.add_argument("--cross_validation", "-cv", help="K value of cross validation", default="-1")
args = parser.parse_args()

if __name__ == '__main__':
    try:
        em.run(args.feature, args.bootstraps, args.c1, args.beta, args.input, args.cross_validation)
    except Exception as e:
        print(e)
