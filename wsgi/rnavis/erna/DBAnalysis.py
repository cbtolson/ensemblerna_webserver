##########################################################################################################
#Functions for the analysis of dot-bracket structures
    #encodeStructsNested(db_file)
    #getNestFreq(bin_mat_2d)
##########################################################################################################
import re
from itertools import groupby, count
from operator import itemgetter
import numpy as np
from collections import Counter

##########################################################################################################
#Checks parentheses
#Input: string
#Output: indices of unpaired bracekts
##########################################################################################################
def _parChecker(symbolString):
    s = []
    balanced = True
    index = 0
    ind = []
    while index < len(symbolString):
        symbol = symbolString[index]
        if symbol == "(":
            s.append(index)
        else:
            if not s:
                balanced = False
                ind = ind+[index]
            else:
                s = s[:-1]
        index = index + 1
    if balanced and not s:
        return True
    else:
        return ind+s


##########################################################################################################
#Finds inner loops
#Input: right index, left index, string
#Output: left index, right index, length
##########################################################################################################
def _findLoop(left, right, s):
    if right+1 > len(s)-1 or left-1 < 0:
        loop = [left, right, s[left:right].count("(")*2]
        return loop
    elif s[right+1] == '(' or s[left-1] == ')':
        loop = [left, right, s[left:right].count("(")*2]
        return loop
    elif s[right+1] == '.' or s[left-1] == '.':
        loop = [left, right, s[left:right].count("(")*2]
        return  loop
    else:
        loop = _findLoop(left-1, right+1, s)
        return loop


##########################################################################################################
#Function converting dot-bracket notation files into 2D matrices with neded secondary structure information
#Input: dot-bracket file header
#Output: 2D matrix, list of sequences
##########################################################################################################
def encodeStructsNested(db_file, rg):

    #Read structure
    file = open(db_file+".db","r")
    lines = file.readlines()
    seq = [x.strip() for x in lines]

    #initialize variables
    n = len(seq)
    bin_mat_2d = [[] for i in range(n)]
    number = 0
    ig = itemgetter(0, -1)

    #loop through each structure in ensemble
    for i in seq:
        #initialize variables
        ig = itemgetter(0, -1)
        s = i[rg[0]:((rg[len(rg)-1])+1)]
        s = s.replace('.', '')
        pattern = []
        
        if s == '':
            arr2d = [0]
        else:
            #remove unmatched parentheses
            p = _parChecker(s)
            if p != True:
                if(isinstance(p,int)):
                    p = [p]
                if len(p)>0:
                    s = list(s)
                    g = [s.pop(x) for x in p[::-1]]
                    s = ''.join(s)

            #find nested loops
            while s.find('(') > -1:
                left = [x.start() for x in re.finditer('\((\.)*\)', s)][0]
                right = [x.end()-1 for x in re.finditer('\((\.)*\)', s)][0]
                loop = _findLoop(left, right, s)
                pattern += [loop]
                s = ''.join(['.'  if x in range((loop[0]),(loop[1]+1)) else s[x] for x in range(len(s))])
            pattern = sorted(pattern, key=lambda x: x[1])
        
            #remove small loops
            pattern = [pattern[x] for x in range(len(pattern)) if pattern[x][2] > 4][::-1]
            
            #encode into vector
            vec = []
            while len(pattern)>0:
                left = pattern[0][0]
                right = pattern[0][1]
                count = [x for x in range(len(pattern)) if pattern[x][0] >= left and pattern[x][1] <= right]
                for x in count[::-1]:
                    pattern.pop(x)
                vec += [len(count)-1]
            vec = vec[::-1]
    
            #format vector
            if(all(x<2 for x in vec)):
                arr2d = [len(vec)]
            else:
                vec = [x if x!=0 else x+1 for x in vec]
                arr2d = [0]+vec
                        
        bin_mat_2d[number] = arr2d
        number = number+1
            
    #reshape array
    maxlen = max([len(i) for i in bin_mat_2d])
    for i in range(len(bin_mat_2d)):
        diff = maxlen-len(bin_mat_2d[i])
        if diff>0:
            for j in range(diff):
                bin_mat_2d[i].append(0)

    #close file
    file.close()

    return{'norm2d_':bin_mat_2d, 'db_':seq}


##########################################################################################################
#Function to get frequency of nested structures
#Input: 2D matrix
#Output: lists for vector representation, frequency, sequency, representative for structure clusters
##########################################################################################################
def getNestFreq(bin_mat_2d):
    
    #initialize variables
    seqs = [[] for i in range(len(bin_mat_2d))]

    #convert to string
    maxlen = max([len(i) for i in bin_mat_2d])
    for i in range(len(bin_mat_2d)):
        seqs[i] = ','.join(map(str, bin_mat_2d[i]))

    #get frequency
    d = Counter(seqs)
    d = sorted(d.items())
    n = [x[0] for x in d]
    m = [x[1] for x in d]

    #get arrays
    arr = [[] for i in range(len(n))]
    s = [[] for i in range(len(n))]
    n = list(n)
    for i in range(len(n)):
        arr[i] = bin_mat_2d[seqs.index(n[i])]

    return{'nest':n, 'freq':m, 'seqs':seqs, 'arr':arr}













