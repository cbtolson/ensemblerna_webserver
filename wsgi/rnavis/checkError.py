##########################################################################################################
#Functions to check input errors
    #getTaskID
    #checkFasta
    #checkSHAPE
    #checkRange
    #checkSize
    #checkMaxD
    #checkDB
##########################################################################################################
import os
import re
import uuid
from flask import Flask
from pymongo import MongoClient
from gridfs import *

##########################################################################################################
#Function for getting task id
##########################################################################################################
def getTaskID(app):
    
    #initialize variables
    task_id = uuid.uuid4()
    tmpdir = os.path.join(app.config['TEMP_FOLDER'], str(task_id))
    
    #check if task id exists
    client = MongoClient(os.environ['OPENSHIFT_MONGODB_DB_URL'])
    db = client['ensemblerna']
    fs = GridFS(db)
    while((fs.exists(filename="interactive.html", task_id=task_id) == True) or (os.path.isdir(tmpdir) == True)):
        task_id = uuid.uuid4()
        tmpdir = os.path.join(app.config['TEMP_FOLDER'], str(task_id))
    client.close()

    return str(task_id)


##########################################################################################################
#Function to check fasta file
#Input: fasta file
#Output: sequence
##########################################################################################################
def checkFasta(filename, length=None):

    #check if file exists
    if os.path.exists(filename) == False:
        error = "ERROR: " + filename + " does not exist."
        return{'tf':error}

    #read in file
    with open(filename, 'r') as myfile:
        fasta_file = myfile.read()
    myfile.close()
    ff = fasta_file.split("\n")
    
    #Error for wrong file structure
    if(ff[0].startswith(">")!=True):
        error = "ERROR: Check fasta file formatting."
        return{'tf':error}
    
    #Error for wrong file content
    ffl = ff.pop(0)
    ffl = ''.join(ff)
    ffl = ''.join(ffl.split())
    if(all(c in "ACGTUXacgtux" for c in ffl)!=True):
        error = "ERROR: Check fasta file formatting."
        return{'tf':error}
    
    #Error for length
    if(len(ffl)>2500 or len(ffl)<1):
        error = "ERROR: Check sequence length."
        return{'tf':error}

    #Error for mismatching lengths
    if length is not None:
        if(len(ffl) != length):
            error = "ERROR: Reference and map sequences must have same length."
            return{'tf':error}

    #Check maximum length
    if len(ffl) > 400:
        error = "At this time the EnsembleRNA Webserver is limited to sequences shorter than 400 nts. For longer RNAs, please download the EnsembleRNA Python Package."
        return{'tf':error}

    return{'tf':True, 'length':len(ffl)}


##########################################################################################################
#Function to check shape file
#Input: shape file
#Output: True/False
##########################################################################################################
def checkSHAPE(filename):
    
    #check if file exists
    if os.path.exists(filename) == False:
        error = "ERROR: " + filename + " does not exist."
        return(error)
    
    #Error for empty file
    if os.stat(filename).st_size == 0:
        error = "ERROR: Check shape file."
        return(error)

    #read in file
    with open(filename, 'r') as myfile:
        shape_file = myfile.read()
    myfile.close()
    sf = re.split('\r|\n|\r\n', shape_file)

    #Error for wrong content
    sfl = ''.join(sf)
    sfl = ''.join(sfl.split())
    if(all(c in "1234567890.-" for c in sfl)!=True):
        error = "ERROR: Check shape file formatting."
        return(error)

    #Error for wrong file structure
    sf = re.split('\r|\n|\r\n', shape_file)
    sf = [x  for x in filter(None, sf)]
    if(all(len(x.split())==2 for x in sf)==True):
       if(all(x.split()[0].isdigit() for x in sf)!=True):
            error = "Check shape file formatting."
            return(error)
    else:
        error = "ERROR: Check shape file formatting."
        return(error)

    #return true if checks are passed
    return True
       
##########################################################################################################
#Function to check db file
#Input: db file
#Output: True/False
##########################################################################################################
def checkDB(filename, length=None):
       
    #check if file exists
    if os.path.exists(filename) == False:
        error = "ERROR: " + filename + " does not exist."
        return(error)
    
    #Error for empty file
    if os.stat(filename).st_size == 0:
        error = "ERROR: Check dot-bracket file."
        return(error)
       
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
        error = "ERROR: Check dot-bracket formatting."
        return(error)

    #Error for mismatching lengths
    if length is not None:
        if(all(len(x)==length for x in dbf)!=True):
            error = "ERROR: Reference sequence and dot-bracket structures must be the same length"
            return(error)

    #return true if checks are passed
    return True

##########################################################################################################
#Function to check size
#Input: size, sequence length
#Output: True/False
##########################################################################################################
def checkSize(size, n, rg):
    if size < 1:
        error = "ERROR: Map size must be larger than 0."
        return(error)
    if size > n:
        error = "ERROR: Map size cannot be larger than sequence length"
        return(error)
    if size > len(rg):
        error = "ERROR: Map size cannot be larger than range"
        return(error)
    if size > 50:
        error = "At this time the EnsembleRNA Webserver is limited to a map size of 50 sequences. For larger map sizes, please download the EnsembleRNA Python Package."
        return(error)

    return size

##########################################################################################################
#Function to thread maximum
#Input: thread maximum
#Output: True/False
##########################################################################################################
def checkMaxD(maxd):
    
    #check maximum distance
    if maxd < 1:
        error = "ERROR: Maximum distance between pairs must be larger than 0."
        return(error)
    
    return maxd

##########################################################################################################
#Function to check range
#Input: output directory
#Output: True/False
##########################################################################################################
def checkRange(nucrg, length):

    #check range length
    if nucrg[0] < 1:
        error = "ERROR: Start of range cannot be smaller than 1."
        return(error)
    if nucrg[1] > length:
        error = "ERROR: End of range cannot exceed the length of the RNA."
        return(error)

    #check range
    if nucrg[1]<nucrg[0]:
        error = "ERROR: End of range cannot begin before start of range."
        return(error)

    #calculate range
    rg = range(nucrg[0]-1, nucrg[1], 1)

    return rg

