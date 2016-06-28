##########################################################################################################
#Functions generating dot-bracket structures
    #getMapSeqs
    #getMapStruct
    #getStruct
    #getStructSHAPE
##########################################################################################################
import os
import numpy
import string
import re
import subprocess
import csv
import scipy.cluster as spcl
import scipy.spatial as spsp
from threading import Thread
from threading import BoundedSemaphore

##########################################################################################################
#Function to determine structures with SHAPE constraints
#Input: directiory, fasta file, shape file, outfile
#Output: db file
##########################################################################################################
def _con2db(dir, con_file, outfile):
    
    #initialize variables
    db = None

    #open input files
    with open(dir+con_file, 'r') as myfile:
        cfile = myfile.read()
    myfile.close()

    #open output files
    with open(dir+outfile, 'w') as dbfile:
        writer = csv.writer(dbfile)

        #loop through structures
        line = re.split('\r|\n|\r\n', cfile)
        line = [x for x in filter(None, line)]
        for i in range(0,len(line)):
            if 'Seq' in line[i]:
                #print previous structure
                if(db is not None):
                    db = ''.join(db)
                    writer.writerow([db])

                #start new structure
                arr = line[i].split()
                numOfChar = int(arr[0])
                db = '.' * numOfChar
                db = list(db)
            else:
                #create dot bracket
                arr = line[i].split()
                left = int(arr[0])
                right = int(arr[4])
                if right != 0:
                    db[left-1] = ")"
                    db[right-1] = "("

        db = ''.join(db)
        writer.writerow([db])

    return True


##########################################################################################################
#Function to get base pairing probability
#Input: directory, fasta file header, number of sequences in map
#Output: 1D array
##########################################################################################################
def _getBasePairing(dir, j, vec, se, i, rsp, threadLimiter, SCRIPT_FOLDER):
    #limit thread
    threadLimiter.acquire()
    
    #set datapath
    os.environ["DATAPATH"] = SCRIPT_FOLDER+'RNAstructure/data_tables/'
    
    try:
        #creates partition function
        cmd = SCRIPT_FOLDER+'RNAstructure/exe/partition ' + dir + 'temp_'+str(j+1)+'.seq ' + dir + 'temp_'+str(j+1)+'.pfs'+rsp
        subprocess.check_output(cmd, shell=True)
        
        #check parition function (for first time)
        if os.path.exists(dir + '/temp_'+str(j+1)+'.pfs') == False:
            vec[i] = None
            se[i] = None
            return True
        
        #gets base pairing probabilites
        cmd = SCRIPT_FOLDER+'RNAstructure/exe/ProbabilityPlot ' + dir + 'temp_'+str(j+1)+'.pfs ' + dir + 'temp_'+str(j+1)+'.txt -t'
        subprocess.check_output(cmd, shell=True)

        #check base pairing probabilities (for first time)
        with open(dir+'temp_'+str(j+1)+'.txt', 'r') as myfile:
            prob_file = myfile.read()
        myfile.close()
        prob = re.split('\r|\n|\r\n', prob_file)
        prob = filter(None, prob)
        prob = [x  for x in filter(None, prob)]
        if len(prob) < 3:
            vec[i] = None
            se[i] = None
            return True

        #calculate base pairing probabilities
        table = numpy.loadtxt(dir + 'temp_'+str(j+1)+'.txt', skiprows=2)
        values = table[:,2]
        groups = table[:,0]
        order = numpy.argsort(groups)
        groups = groups[order]
        values = values[order]
        values.cumsum(out=values)
        index = numpy.ones(len(groups), 'bool')
        index[:-1] = groups[1:] != groups[:-1]
        values = values[index]
        groups = groups[index]
        values[1:] = values[1:] - values[:-1]
        diff = numpy.setdiff1d(list(range(1,(len(vec)+1))), groups.tolist())
        nt = numpy.append(groups, diff)
        prob = numpy.append(values, numpy.zeros(len(diff)))
        bp2 = numpy.append([nt], [prob], axis=0)
        bp2 = bp2[numpy.lexsort(numpy.fliplr(bp2).T)]
        bp2 = bp2[1,:]
        vec[i] = [x if x == 0 else 10**(-x) for x in bp2]

        #get shannon entropy
        se[i] = -sum([x if x == 0 else x*numpy.log(x) for x in bp2])/numpy.log(len(vec))
    finally:
        #limit thread
        threadLimiter.release()


##########################################################################################################
#Function to determine sequences for the map of conformational space
#Input: directory, fasta file header, number of sequences in map
#Output: 1D array
##########################################################################################################
def getMapSeqs(dir, fasta_file, n, rg, rsp, SCRIPT_FOLDER, thmax=1):
    
    #Read structure from fasta file
    lines = fasta_file.split("\n")
    seq = lines.pop(0)
    
    #initialize parameters
    threadLimiter = BoundedSemaphore(thmax)
    vec = [[] for i in range(len(rg)+1)]
    se = [[] for i in range(len(rg)+1)]
    mutseqs = ['' for i in range(len(rg)+1)]
    threads = [[] for i in range(len(rg)+1)]
    indices = [0]+[i+1 for i in rg]
    
    #loop through each nucleotide in the sequence
    for i in range(len(rg)+1):
        #initialize parameters
        mut = seq.upper()
        
        mutseqs[i] = mut
        if(i < len(rg)):
            #initialize parameters
            j = rg[i]
                        
            #create vector of nucleotides not in WT
            if mut[j] == "A":
                nucs = "U"
            elif mut[j] == "U":
                nucs = "A"
            elif mut[j] == "G":
                nucs = "C"
            else:
                nucs = "G"

            #create single point mutant
            mutseqs[i] = mut[:j] + re.sub('.', nucs, mut[j]) + mut[j+1:]
        else:
            j = -1
        
        #create seq file
        cmd = 'echo \;seq\ > ' + dir + 'temp_'+str(j+1)+'.seq'
        subprocess.check_output(cmd, shell=True)
        cmd = 'echo \Seq\ >> ' + dir + 'temp_'+str(j+1)+'.seq'
        subprocess.check_output(cmd, shell=True)
        cmd = 'echo ' + mutseqs[i] + '1 >> ' + dir + 'temp_'+str(j+1)+'.seq'
        subprocess.check_output(cmd, shell=True)
        
        #get base pairing probability and shannon entropy
        threads[i] = Thread(target=_getBasePairing, args=(dir, j, vec, se, i, rsp, threadLimiter, SCRIPT_FOLDER))
        threads[i].start()
    
    for i in range(len(threads)):
        threads[i].join()

    #filter sequences that don't fold
    vec = vec[:]
    se = se[:]
    ind = [i for i in range(len(se))]
    ind = [x for x in ind if(vec[x] is not None and se[x] is not None)]
    if(len(ind)<n):
        cmd = 'rm -f ' + dir + '/*'
        subprocess.check_output(cmd, shell=True)
        raise ValueError("Map cannot be created. Not enough map sequences are properly folding.")
    sep = [se[x] for x in ind if(vec[x] is not None and se[x] is not None)]

    #filter low Shannon entopy
    p = numpy.percentile(sep, 25)
    ip = [i for i in ind if se[i]>p]
    if len(ip) > n:
        ind = ip
    if(len(ind)<n):
        cmd = 'rm -f ' + dir + '*'
        subprocess.check_output(cmd, shell=True)
        raise ValueError("Map cannot be created. Not enough map sequences are properly folding.")
    vec = numpy.array([vec[i] for i in ind])
    mutseqs = [mutseqs[i] for i in ind]
    indices = [indices[i] for i in ind]

    #hierarchical cluster
    count = 0
    d = spsp.distance.pdist(vec)
    hc = spcl.hierarchy.linkage(d)
    tree = spcl.hierarchy.fcluster(hc, n, 'maxclust')
    num = len(set(tree))
    mapseqs = ['' for i in range(num)]
    mapinds = ['' for i in range(num)]
    for i in range(1,num+1):
        clust = [j for j in range(len(tree)) if tree[j]==i]
        ind = clust[0]
        mapseqs[count] = mutseqs[ind]
        mapinds[count] = indices[ind]
        count = count+1
    
    #removes temporary files
    cmd = 'rm ' + dir + 'temp*.txt'
    subprocess.check_output(cmd, shell=True)
    cmd = 'rm ' + dir + 'temp*.seq'
    subprocess.check_output(cmd, shell=True)

    return{'mapseqs':mapseqs, 'mapinds':mapinds}

##########################################################################################################
#Function to determine structures for the map of conformational space
#Input: 1D array, outfile
#Output: db file
##########################################################################################################
def getMapStruct(dir, mapinds, outfile, SCRIPT_FOLDER):
    
    #set datapath
    os.environ["DATAPATH"] = SCRIPT_FOLDER+'RNAstructure/data_tables/'
    
    #removes exsiting file
    if(os.path.exists(dir + outfile)):
        cmd = 'rm ' + dir + outfile
        subprocess.check_output(cmd, shell=True)

    #loop through each of map sequences
    for i in mapinds:
        
        #runs stochastic without SHAPE contraints
        cmd = SCRIPT_FOLDER+'RNAstructure/exe/stochastic ' + dir + 'temp_'+str(i)+'.pfs ' + dir + 'temp.con -e 1000'
        subprocess.check_output(cmd, shell=True)
        
        #converts to dot bracket notation
        _con2db(dir, 'temp.con', 'temp.txt')
        cmd = 'cat ' + dir + 'temp.txt >> ' + dir + outfile
        subprocess.check_output(cmd, shell=True)

##########################################################################################################
#Function to determine structures
#Input: directory, fasta file header, cluster number
#Output: db file
##########################################################################################################
def getStruct(dir, outfile, SCRIPT_FOLDER):
    
    #set datapath
    os.environ["DATAPATH"] = SCRIPT_FOLDER+'RNAstructure/data_tables/'

    #runs stochastic without SHAPE contraints
    cmd = SCRIPT_FOLDER+'RNAstructure/exe/stochastic ' + dir + 'temp.pfs ' + dir + 'temp.con -e 1000'
    subprocess.check_output(cmd, shell=True)
    
    #converts to dot bracket notation
    _con2db(dir, 'temp.con', outfile)
    
    #removes temporary files
    cmd = 'rm ' + dir + 'temp*'
    subprocess.check_output(cmd, shell=True)

##########################################################################################################
#Function to determine structures with SHAPE constraints
#Input: directiory, fasta file, shape file, outfile
#Output: db file
##########################################################################################################
def getStructSHAPE(dir, outfile, SCRIPT_FOLDER):
    
    #set datapath
    os.environ["DATAPATH"] = SCRIPT_FOLDER+'RNAstructure/data_tables/'
    
    #runs stochastic with SHAPE contraints
    cmd = SCRIPT_FOLDER+'RNAstructure/exe/stochastic ' + dir + 'temp.pfs ' + dir + 'temp.con -e 1000'
    subprocess.check_output(cmd, shell=True)

    #converts to dot bracket notation
    _con2db(dir, 'temp.con', outfile)

    #removes temporary files
    cmd = 'rm ' + dir + 'temp*'
    subprocess.check_output(cmd, shell=True)



