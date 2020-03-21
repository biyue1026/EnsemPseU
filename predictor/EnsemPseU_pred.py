#!/usr/bin/env python
# _*_coding:utf-8_*_


import itertools
import pickle
import argparse
import os,sys,re
import numpy as np
from collections import Counter
import joblib

def Binary(sequences):
    AA = 'ACGU'
    binary_feature = []
    for seq in sequences:
        binary = []
        for aa in seq:
            for aa1 in AA:
                tag = 1 if aa == aa1 else 0
                binary.append(tag)
        binary_feature.append(binary)
    return binary_feature



def NCP(sequences):
    chemical_property = {
        'A': [1, 1, 1],
        'C': [0, 1, 0],
        'G': [1, 0, 0],
        'U': [0, 0, 1], }
    ncp_feature = []
    for seq in sequences:
        ncp = []
        for aaindex, aa in enumerate(seq):
            ncp = ncp + chemical_property.get(aa, [0, 0, 0])
        ncp_feature.append(ncp)
    return ncp_feature


def ND(sequences):
    nd_feature = []
    for seq in sequences:
        nd = []
        for aaindex, aa in enumerate(seq):
            nd.append(seq[0: aaindex + 1].count(seq[aaindex]) / (aaindex + 1))
        nd_feature.append(nd)
    return nd_feature



def ENAC(sequences):
    AA = 'ACGU'
    enac_feature = []
    window = 5
    for seq in sequences:
        l = len(seq)
        enac= []
        for i in range(0, l):
            if i < l and i + window <= l:
                count = Counter(seq[i:i + window])
                for key in count:
                    count[key] = count[key] / len(seq[i:i + window])
                for aa in AA:
                    enac.append(count[aa])
        enac_feature.append(enac)
    return enac_feature



def Kmer123(sequences):
    AA = 'ACGU'
    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i
    kmer123_feature = []

    for seq in sequences:
        kmer1 = [0] * 4
        for j in range(len(seq)):
            kmer1[AADict[seq[j]]] = kmer1[AADict[seq[j]]] + 1
        if sum(kmer1) != 0:
            kmer1 = [i / sum(kmer1) for i in kmer1]

        kmer2 = [0] * 16
        for j in range(len(seq) - 2 + 1):
            kmer2[AADict[seq[j]] * 4 + AADict[seq[j + 1]]] = kmer2[AADict[seq[j]] * 4 + AADict[seq[j + 1]]] + 1
        if sum(kmer2) != 0:
            kmer2 = [i / sum(kmer2) for i in kmer2]

        kmer3 = [0] * 64
        for j in range(len(seq) - 3 + 1):
            kmer3[AADict[seq[j]] * 16 + AADict[seq[j + 1]] * 4 + AADict[seq[j + 2]]] = kmer3[AADict[seq[j]] * 16 + AADict[seq[j + 1]] * 4 + AADict[seq[j + 2]]] + 1
        if sum(kmer3) != 0:
            kmer3 = [i / sum(kmer3) for i in kmer3]

        kmer = kmer1+kmer2 + kmer3
        kmer123_feature.append(kmer)
    return kmer123_feature



def Kmer4(sequences):
    AA = 'ACGU'
    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i
    kmer4 = []
    for seq in sequences:
        tmpCode = [0] * 256
        for j in range(len(seq) - 4 + 1):
            tmpCode[AADict[seq[j]] * 64 + AADict[seq[j + 1]] * 16 + AADict[seq[j + 2]] * 4 + AADict[seq[j + 3]]] = tmpCode[AADict[seq[j]] * 64 + AADict[seq[j + 1]] * 16 + AADict[seq[j + 2]] * 4 + AADict[seq[j + 3]]] + 1
        if sum(tmpCode) != 0:
            tmpCode = [i / sum(tmpCode) for i in tmpCode]
        kmer4.append(tmpCode)
    return kmer4



def Kmer5(sequences):
    AA = 'ACGU'
    kmer5 = []
    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i
    for seq in sequences:
        tmpCode = [0] * 1024
        for j in range(len(seq) - 5 + 1):
            tmpCode[AADict[seq[j]] * 256 + AADict[seq[j + 1]] * 64 + AADict[seq[j + 2]] * 16 + AADict[seq[j + 3]] * 4 +
                    AADict[seq[j + 4]]] = tmpCode[
                                              AADict[seq[j]] * 256 + AADict[seq[j + 1]] * 64 + AADict[seq[j + 2]] * 16 +
                                              AADict[seq[j + 3]] * 4 + AADict[seq[j + 4]]] + 1
        if sum(tmpCode) != 0:
            tmpCode = [i / sum(tmpCode) for i in tmpCode]
        kmer5.append(tmpCode)

    return kmer5



def read_fasta(inputfile):
    if os.path.exists(inputfile) == False:
        print('Error: file " %s " does not exist.' % inputfile)
        sys.exit(1)
    with open(inputfile) as f:
        record = f.readlines()
    if re.search('>',record[0]) == None:
        print('Error: the input file " %s " must be fasta format!' % inputfile)
        sys.exit(1)

    data = {}
    for line in record:
        if line.startswith('>'):
            name = line.replace('>','').split('\n')[0]
            data[name] = ''
        else:
            data[name] += line.replace('\n','')
    return data




def extract_features(data):
    sequences = data
    basic1 = np.array(Binary(sequences))
    basic2 = np.array(Kmer123(sequences))
    basic3 = np.array(Kmer4(sequences))
    basic4 = np.array(Kmer5(sequences))
    basic5 = np.array(ENAC(sequences))
    basic6 = np.array(NCP(sequences))
    basic7 = np.array(ND(sequences))
    feature_vector = np.concatenate((basic1,basic2,basic3,basic4,basic5,basic6,basic7), axis=1)
    return feature_vector


def select_feature(feature_vector,type):
    if os.path.exists('selector/selector_'+type) == False:
        s_feature_vector = feature_vector
    else:
        selector = []
        with open('selector/selector_'+type) as fff:
            for line in fff:
                selector.append(int(line.strip()))
        s_feature_vector = []
        for i in feature_vector:
            i = i[selector]
            s_feature_vector.append(i)
        s_feature_vector = np.array(s_feature_vector)
    return s_feature_vector





def predict_site(data,outputfile,type):
    bs = joblib.load('model/'+type+'.pkl')
    vector = extract_features(data.values())
    s_vector = select_feature(vector,type)
    predictions = bs.predict(s_vector)
    probability = ['%.5f' % float(i) for i in predictions]
    name = list(data.keys())
    seq = list(data.values())
    with open(outputfile,'w') as f:
        for i in range(len(data)):
            if float(probability[i]) > 0.5:
                f.write(probability[i] + '\t')
                f.write(name[i] + '\t')
                f.write(seq[i] + '\n')
            else:
                f.write(probability[i] + '\t')
                f.write(name[i] + '\t')
                f.write(seq[i] + '\n')
    return None





def main():
    parser = argparse.ArgumentParser(description='EnsemPseU: Identifying pseudouridine sites based on a voted ensemble method')
    parser.add_argument('--species',dest='type',type=str,required=True,help='Select the species that you want to predict.'
                    'You can choose "H", "M" and "S", which represent H. sapiens, M. musculus and S. cerevisiae, respectively ')
    parser.add_argument('--input',dest='inputfile',type=str,required=True,help='Query RNA sequences to be predicted in fasta format.'
                        'Please note that the length of input must be 21-bp when the species is H. sapiens or M. musculus, and 31-bp when the species is S. cerevisiae')
    parser.add_argument('--output',dest='outputfile',type=str,required=False,help='Save the prediction results in csv format.')
    args = parser.parse_args()

    inputfile = args.inputfile
    outputfile = args.outputfile
    spe = args.type
    data = read_fasta(inputfile)

    if outputfile != None:
        predict_site(data,outputfile,spe)
        print('Output are saved in ' + outputfile + '.')
    else:
        default_output = 'output'
        predict_site(data, default_output,spe)
        print('Output are saved in the current directory'+ '.')



if __name__ == "__main__":
    main()


