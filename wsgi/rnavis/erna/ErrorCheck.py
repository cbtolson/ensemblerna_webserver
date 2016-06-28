##########################################################################################################
#Functions to check input errors
    #checkDir
    #checkFasta
    #checkRNAStructCMD
    #checkRNAStruct
    #checkRNAStructSHAPE
    #checkSHAPE
    #checkSize
    #checkThmax
    #checkDB
##########################################################################################################
import os
import re
import subprocess
    
##########################################################################################################
#Function to check output directory
#Input: output directory
#Output: True/False
##########################################################################################################
def checkDir(outdir):
    #check if directory exists
    if os.path.isdir(outdir) == False:
        raise ValueError("ERROR: " + outdir + " does not exist.")
    
    if outdir[-1:] != '/':
        outdir = outdir+'/'

    return outdir

##########################################################################################################
#Function to check fasta file
#Input: fasta file
#Output: sequence
##########################################################################################################
def checkFasta(filename, length=None):

    #check if file exists
    if os.path.exists(filename) == False:
        raise ValueError("ERROR: " + filename + " does not exist.")

    #read in file
    with open(filename, 'r') as myfile:
        fasta_file = myfile.read()
    myfile.close()
    ff = fasta_file.split("\n")
    
    #Error for wrong file structure
    if(ff[0].startswith(">")!=True):
        raise ValueError("ERROR: Check fasta file formatting.")
    
    #Error for wrong file content
    ffl = ff.pop(0)
    ffl = ''.join(ff)
    ffl = ''.join(ffl.split())
    if(all(c in "ACGTUXacgtux" for c in ffl)!=True):
        raise ValueError("ERROR: Check fasta file formatting.")
    
    #Error for length
    if(len(ffl)>2500 or len(ffl)<1):
        raise ValueError("ERROR: Check sequence length.")

    #Error for mismatching lengths
    if length is not None:
        if(len(ffl) != length):
            raise ValueError("ERROR: Reference and map sequences must have same length.")

    return ffl


##########################################################################################################
#Function to check RNAstructure
#Input: fasta file, shape file, output directory
#Output: True/False
##########################################################################################################
def checkRNAStruct(fasta, dir, rsp, SCRIPT_FOLDER):
    
    #set datapath
    os.environ["DATAPATH"] = SCRIPT_FOLDER+'RNAstructure/data_tables/'
    
    #create seq file
    cmd = 'echo \;seq\ > ' + dir + '/temp.seq'
    subprocess.check_output(cmd, shell=True)
    cmd = 'echo \Seq\ >> ' + dir + '/temp.seq'
    subprocess.check_output(cmd, shell=True)
    cmd = 'echo ' + fasta + '1 >> ' + dir + '/temp.seq'
    subprocess.check_output(cmd, shell=True)

    #check partition function
    cmd = SCRIPT_FOLDER+'RNAstructure/exe/partition ' + dir + '/temp.seq ' + dir + '/temp.pfs ' + rsp
    out = subprocess.check_output(cmd, shell=True)
    if os.path.exists(dir + '/temp.pfs') == False:
        cmd = 'rm -f ' + dir + '/*'
        subprocess.check_output(cmd, shell=True)
        raise ValueError("ERROR: Partition function cannot be generated. Please check fasta file inputs")
    
    cmd = SCRIPT_FOLDER+'RNAstructure/exe/ProbabilityPlot ' + dir + '/temp.pfs ' + dir + '/temp.txt -t'
    subprocess.check_output(cmd, shell=True)
    with open(dir+'/temp.txt', 'r') as myfile:
        prob_file = myfile.read()
    myfile.close()
    prob = re.split('\r|\n|\r\n', prob_file)
    prob = [x  for x in filter(None, prob)]
    if len(prob) < 3:
        cmd = 'rm -f ' + dir + '/*'
        subprocess.check_output(cmd, shell=True)
        raise ValueError("ERROR: Partition function cannot be generated. Please check fasta file inputs")

    return True


##########################################################################################################
#Function to check RNAstructure
#Input: fasta file, shape file, output directory
#Output: True/False
##########################################################################################################
def checkRNAStructSHAPE(fasta_file, dir, rsp, ssp, SCRIPT_FOLDER):
    
    #set datapath
    os.environ["DATAPATH"] = SCRIPT_FOLDER+'RNAstructure/data_tables/'

    #Read structure from fasta file
    lines = fasta_file.split("\n")
    seq = lines.pop(0)

    #create seq file
    cmd = 'echo \;seq\ > ' + dir + '/temp.seq'
    subprocess.check_output(cmd, shell=True)
    cmd = 'echo \Seq\ >> ' + dir + '/temp.seq'
    subprocess.check_output(cmd, shell=True)
    cmd = 'echo ' + seq + '1 >> ' + dir + '/temp.seq'
    subprocess.check_output(cmd, shell=True)
    
    #check partition function
    cmd = SCRIPT_FOLDER+'RNAstructure/exe/partition ' + dir + 'temp.seq ' + dir + 'temp.pfs -sh ' + dir + 'shapefile.shape' + rsp + ssp
    out = subprocess.check_output(cmd, shell=True)
    if os.path.exists(dir + 'temp.pfs') == False:
        cmd = 'rm -f ' + dir + '/*'
        subprocess.check_output(cmd, shell=True)
        raise ValueError("ERROR: Partition function cannot be generated. Please check fasta file inputs")

    cmd = SCRIPT_FOLDER+'RNAstructure/exe/ProbabilityPlot ' + dir + 'temp.pfs ' + dir + 'temp.txt -t'
    subprocess.check_output(cmd, shell=True)
    with open(dir+'temp.txt', 'r') as myfile:
        prob_file = myfile.read()
    myfile.close()
    prob = re.split('\r|\n|\r\n', prob_file)
    prob = [x  for x in filter(None, prob)]
    if len(prob) < 3:
        cmd = 'rm -f ' + dir + '/*'
        subprocess.check_output(cmd, shell=True)
        raise ValueError("ERROR: Partition function cannot be generated. Please check fasta file inputs")
    
    #removes temporary files
    cmd = 'rm ' + dir + 'shapefile.shape'
    subprocess.check_output(cmd, shell=True)

    return True

##########################################################################################################
#Function to check shape file
#Input: shape file
#Output: True/False
##########################################################################################################
def checkSHAPE(filename, dir):
    
    #check if file exists
    if os.path.exists(filename) == False:
        raise ValueError("ERROR: " + filename + " does not exist.")
    
    #Error for empty file
    if os.stat(filename).st_size == 0:
        raise ValueError("ERROR: Check shape file.")

    #read in file
    with open(filename, 'r') as myfile:
        shape_file = myfile.read()
    myfile.close()
    sf = re.split('\r|\n|\r\n', shape_file)

    #Error for wrong content
    sfl = ''.join(sf)
    sfl = ''.join(sfl.split())
    if(all(c in "1234567890.-" for c in sfl)!=True):
        raise ValueError("ERROR: Check shape file formatting.")

    #Error for wrong file structure
    sf = re.split('\r|\n|\r\n', shape_file)
    sf = [x  for x in filter(None, sf)]
    if(all(len(x.split())==2 for x in sf)==True):
       if(all(x.split()[0].isdigit() for x in sf)!=True):
            raise ValueError("Check shape file formatting.")
    else:
       raise ValueError("ERROR: Check shape file formatting.")
       
    #create shape file
    cmd = 'cat ' + filename + ' > ' + dir+'/shapefile.shape'
    subprocess.check_output(cmd, shell=True)

    #return true if checks are passed
    return True
       
##########################################################################################################
#Function to check db file
#Input: db file
#Output: True/False
##########################################################################################################
def checkDB(filename, dir, length=None):
       
    #check if file exists
    if os.path.exists(filename) == False:
        raise ValueError("ERROR: " + filename + " does not exist.")
    
    #Error for empty file
    if os.stat(filename).st_size == 0:
        raise ValueError("ERROR: Check dot-bracket file.")
       
    #read in file
    with open(filename, 'r') as myfile:
       db_file = myfile.read()
    myfile.close()
       
    #Error for wrong file structure
    dbf = db_file.strip()
    dbf = re.split('\r|\n|\r\n', dbf)
    dbf = [x  for x in filter(None, dbf)]
    df = ''.join(dbf)
    if(all(c in "()." for c in df)!=True):
        cmd = 'rm -f ' + dir + '/*'
        subprocess.check_output(cmd, shell=True)
        raise ValueError("ERROR: Check dot-bracket formatting.")

    #Error for mismatching lengths
    if length is not None:
        if(all(len(x)==length for x in dbf)!=True):
            cmd = 'rm -f ' + dir + '/*'
            subprocess.check_output(cmd, shell=True)
            raise ValueError("ERROR: Reference sequence and dot-bracket structures must be the same length")

    #return true if checks are passed
    return True

##########################################################################################################
#Function to check size
#Input: size, sequence length
#Output: True/False
##########################################################################################################
def checkSize(size, n, rg):
    if size < 1:
        raise ValueError("ERROR: Map size must be larger than 0.")
    if size > n:
        raise ValueError("ERROR: Map size cannot be larger than sequence length")
    if size > len(rg):
        raise ValueError("ERROR: Map size cannot be larger than range")

    return True

##########################################################################################################
#Function to thread maximum
#Input: thread maximum
#Output: True/False
##########################################################################################################
def checkThmax(thmax):
    if thmax < 1:
        raise ValueError("ERROR: Number of threads must be larger than 0.")
    
    return True

##########################################################################################################
#Function to check range
#Input: output directory
#Output: True/False
##########################################################################################################
def checkRange(nucrg, length):

    #check range length
    if nucrg[0] < 1:
        raise ValueError("ERROR: Start of range cannot be smaller than 1.")
    if nucrg[1] > length:
        raise ValueError("ERROR: End of range cannot exceed the length of the RNA.")

    #check range
    if nucrg[1]<nucrg[0]:
        raise ValueError("ERROR: End of range cannot begin before start of range.")

    #calculate range
    rg = range(nucrg[0]-1, nucrg[1], 1)

    return rg

