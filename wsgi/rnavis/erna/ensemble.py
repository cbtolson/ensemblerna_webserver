##########################################################################################################
#Visualization of RNA ensembles
    #getMap
    #getMapDB
    #getRef
    #getRefDB
    #main
##########################################################################################################
import sys
import random
from rnavis.erna.ErrorCheck import *
from rnavis.erna.DBStructs import *
from rnavis.erna.DBAnalysis import *
from rnavis.erna.PlotVis import *

##########################################################################################################
#Function to create map of conformational space
#Input: fasta sequence, file header, output directory, map size
#Output: map positions, map 2D matrix, map sequences, map dot bracket structures
##########################################################################################################
def _getMap(fasta, header, dir, size, plotm, rg, rsp, thmax, SCRIPT_FOLDER):
    
    #initialize variables
    header = re.split('/', header)
    header = header[len(header)-1]
    dir = dir
    
    #get map structure
    map = getMapSeqs(dir, fasta, size, rg, rsp, SCRIPT_FOLDER, thmax)
    mapseqs = map['mapseqs']
    mapinds = map['mapinds']
    getMapStruct(dir, mapinds, header+'_map.db', SCRIPT_FOLDER)
    
    #encode map structure
    bin_mat_2d = encodeStructsNested(dir+header+'_map', rg)
    map2d = bin_mat_2d['norm2d_']
    mapdb = bin_mat_2d['db_']
    
    #get map frequency
    nestfreq = getNestFreq(map2d)
    mapfreq = nestfreq['freq']
    mapnest = nestfreq['nest']
    mapseqs = nestfreq['seqs']
    maparr = nestfreq['arr']
    
    #plot map
    mappos = plotMap(maparr, mapfreq, mapnest, mapseqs, mapdb, map2d, dir+header+'_map', plotm)
    
    #return map
    return{'mappos':mappos, 'mapnest':mapnest, 'mapseqs':mapseqs, 'mapdb':mapdb, 'map2d':map2d, 'maparr':maparr}

##########################################################################################################
#Function to create map of conformational space from dot-bracket
#Input: dot-bracket file, file header, output directory, map size
#Output: map positions, map 2D matrix, map sequences, map dot bracket structures
##########################################################################################################
def _getMapDB(mapdb, header, dir, size, plotm, rg):
    
    #initialize variables
    header = re.split('/', header)
    header = header[len(header)-1]
    dir = dir
    
    #get map structure
    cmd = 'cat ' + mapdb + ' > ' + dir+header+'_map.db'
    subprocess.check_output(cmd, shell=True)
    
    #encode map structure
    bin_mat_2d = encodeStructsNested(dir+header+'_map', rg)
    map2d = bin_mat_2d['norm2d_']
    mapdb = bin_mat_2d['db_']
    
    #get map frequency
    nestfreq = getNestFreq(map2d)
    mapfreq = nestfreq['freq']
    mapnest = nestfreq['nest']
    mapseqs = nestfreq['seqs']
    maparr = nestfreq['arr']
    
    #plot map
    mappos = plotMap(maparr, mapfreq, mapnest, mapseqs, mapdb, map2d, dir+header+'_map', plotm)
    
    #return map
    return{'mappos':mappos, 'mapnest':mapnest, 'mapseqs':mapseqs, 'mapdb':mapdb, 'map2d':map2d, 'maparr':maparr}


##########################################################################################################
#Function to create reference visualization
#Input: map object, fasta sequence, file header,output directory, shape data
#Output: csv, db, pdf, png
##########################################################################################################
def _getRef(map, fasta, header, dir, rg, plotint, SCRIPT_FOLDER, shape=None):
    
    #initialize variables
    header = re.split('/', header)
    header = header[len(header)-1]
    dir = dir
    
    #get structure
    if shape is None:
        getStruct(dir, header+'.db', SCRIPT_FOLDER)
    else:
        getStructSHAPE(dir, header+'.db', SCRIPT_FOLDER)

    #read structure
    bin_mat_2d = encodeStructsNested(dir+header, rg)
    ref2d = bin_mat_2d['norm2d_']
    refdb = bin_mat_2d['db_']
    
    #get reference frequency
    nestfreq = getNestFreq(ref2d)
    reffreq = nestfreq['freq']
    refnest = nestfreq['nest']
    refseqs = nestfreq['seqs']
    refarr = nestfreq['arr']

    #get map
    mappos = map['mappos']
    mapnest = map['mapnest']
    mapseqs = map['mapseqs']
    mapdb = map['mapdb']
    map2d = map['map2d']
    maparr = map['maparr']
    
    #plot reference
    ref = plotRef(reffreq, refnest, refarr, mappos, mapnest, mapseqs, mapdb, maparr, map2d, dir+header, refdb, refseqs, rg)
    #plot interactive
    if plotint == 'T':
        structs = ref['structs']
        diversity = ref['diversity']
        freq = ref['freq']
        plotInteractive(dir, header, mappos, freq, structs, diversity, maparr, rg)

    return True

##########################################################################################################
#Function to create reference visualization with dot bracket
#Input: map object, dot bracket, file header,output directory
#Output: csv, db, pdf, png
##########################################################################################################
def _getRefDB(map, db, header, dir, rg, plotint, md=None):
    
    #initialize variables
    header = re.split('/', header)
    header = header[len(header)-1]
    dir = dir
    
    #remove map temporary files
    if md is None:
        cmd = 'rm ' + dir + '/temp*'
        subprocess.check_output(cmd, shell=True)
    
    #get reference structure
    cmd = 'cat ' + db + ' > ' + dir+header+'.db'
    subprocess.check_output(cmd, shell=True)

    #read structure
    bin_mat_2d = encodeStructsNested(dir+header, rg)
    ref2d = bin_mat_2d['norm2d_']
    refdb = bin_mat_2d['db_']
    
    #get reference frequency
    nestfreq = getNestFreq(ref2d)
    reffreq = nestfreq['freq']
    refnest = nestfreq['nest']
    refseqs = nestfreq['seqs']
    refarr = nestfreq['arr']
    
    #get map
    mappos = map['mappos']
    mapnest = map['mapnest']
    mapseqs = map['mapseqs']
    mapdb = map['mapdb']
    map2d = map['map2d']
    maparr = map['maparr']
    
    #plot reference
    ref = plotRef(reffreq, refnest, refarr, mappos, mapnest, mapseqs, mapdb, maparr, map2d, dir+header, refdb, refseqs, rg)

    #plot interactive
    if plotint == 'T':
        structs = ref['structs']
        diversity = ref['diversity']
        freq = ref['freq']
        plotInteractive(dir, header, mappos, freq, structs, diversity, maparr, rg)
    
    return True


#runs main part of script
def ensembleRNA(fasta_file, outdir, shape_file, db_file, map_file, map_dbfile, size, nucrg, maxd, SCRIPT_FOLDER):

    #set seed
    random.seed(113)

    #preset variables
    sint = None
    slope = None
    tempcalc = None
    plotm = 'T'
    plotint = 'T'
    thmax = 2

    #initialize variables
    header = fasta_file.split('/')[-1]
    header = header.split('.')[0]

    #edit outdir
    if(outdir[-1] == '/'):
        outdir = outdir[:-1]

    #read fasta file
    fasta = checkFasta(fasta_file)

    #check inputs for errors
    outdir = checkDir(outdir)
    if nucrg == None:
        rg = range(0,len(fasta),1)
    else:
        rg = checkRange(nucrg, len(fasta))
    checkSize(size, len(fasta), rg)
    checkThmax(thmax)

    #check RNAstructure parameters
    rsp = ''
    ssp = ''
    if maxd is not None:
        rsp = rsp+' -md '+str(maxd)
    if tempcalc is not None:
        rsp = rsp+' -t '+str(tempcalc)
    if sint is not None:
        ssp = ssp+' -si '+str(sint)
    if sint is not None:
        ssp = ssp+' -sm '+str(slope)

    #check reference
    if db_file is None:
        #choose shape function (shape or none)
        if shape_file is None:
            shape = None
            checkRNAStruct(fasta, outdir, rsp, SCRIPT_FOLDER)
        else:
            shape = 1
            checkSHAPE(shape_file, outdir)
            checkRNAStructSHAPE(fasta, outdir, rsp, ssp, SCRIPT_FOLDER)
    else:
        #check dot bracket
        checkDB(db_file, outdir, len(fasta))

    #choose map function (db or fasta)
    if map_dbfile is not None:
        #check dot bracket
        checkDB(map_dbfile, outdir, len(fasta))
        
        #get map
        map = _getMapDB(map_dbfile, header, outdir, size, plotm, rg)
    else:
        #read map file
        if map_file is not None:
            mfasta = checkFasta(map_file, len(fasta))
        else:
            mfasta = fasta

        #get map
        map = _getMap(mfasta, header, outdir, size, plotm, rg, rsp, thmax, SCRIPT_FOLDER)

    #choose reference function (db or fasta)
    if db_file is not None:
        #visualize reference
        if map_dbfile is None:
            ref = _getRefDB(map, db_file, header, outdir, rg, plotint)
        else:
            ref = _getRefDB(map, db_file, header, outdir, rg, plotint, md=1)
    else:
        #visualize reference
        ref = _getRef(map, fasta, header, outdir, rg, plotint, SCRIPT_FOLDER, shape)

    return(True)


