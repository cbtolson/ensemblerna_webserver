##########################################################################################################
#CELERY TASK
##########################################################################################################
import os
import subprocess
import smtplib
from pymongo import MongoClient
from gridfs import *
import gzip
from email.mime.text import MIMEText
from celery import Celery
from celery.result import AsyncResult
from rnavis.erna.ensemble import ensembleRNA

#create celery app
capp = Celery('tasks', backend='amqp', broker=os.getenv('CLOUDAMQP_URI'), include=['rnavis.tasks'])

#celery task
@capp.task(trail=True)
def runEnsembleRNA(task_id, ref_fasta, ref_shape, ref_db, map_fasta, map_db, rgstart, rgend, size, maxd, email, TEMP_FOLDER, SCRIPT_FOLDER):
    
    #initialize variables
    fasta_file = os.path.join(TEMP_FOLDER, ref_fasta)
    outdir = os.path.join(TEMP_FOLDER,task_id)
    
    #add reference flags
    if ref_db is not None:
        db_file = os.path.join(TEMP_FOLDER, ref_db)
    else:
        db_file = None
    if ref_shape is not None:
        shape_file = os.path.join(TEMP_FOLDER, ref_shape)
    else:
        shape_file = None
    
    #add map flags
    if map_db is not None:
        map_dbfile = os.path.join(TEMP_FOLDER, map_db)
    else:
        map_dbfile = None
    if map_fasta is not None:
        map_file = os.path.join(TEMP_FOLDER, map_fasta)
    else:
        map_file = None
    
    #add advanced flags
    nucrg = [int(rgstart), int(rgend)]
    size = int(size)
    if maxd is not None:
        maxd = int(maxd)
    else:
        maxd = None
    
    #run ensemblerna
    try:
        output = ensembleRNA(fasta_file, outdir, shape_file, db_file, map_file, map_dbfile, size, nucrg, maxd, SCRIPT_FOLDER)
    except BaseException as exception:
        output = str(exception)
    
    #check output
    if isinstance(output, str):
        
        if email is not None and email != '' and email != 'None':
            #email message
            TEXT = "Sorry, an error occured and we were unable to generate your visualization. Please check your input files or contact us at laederachlab@gmail.com (Please do not reply to this email)\n\n"+output
        
            #send email
            gmail_user = "ensemblerna@gmail.com"
            gmail_pwd = "ribosnitch"
            FROM = "ensemblerna@gmail.com"
            TO = [email]
            SUBJECT = "EnsembleRNA Output is Ready"
            message = """\From: %s\nTo: %s\nSubject: %s\n\n%s """ % (FROM, ", ".join(TO), SUBJECT, TEXT)
            try:
                server = smtplib.SMTP("smtp.gmail.com", 587)
                server.ehlo()
                server.starttls()
                server.login(gmail_user, gmail_pwd)
                server.sendmail(FROM, TO, message)
                server.close()
            except:
                return "ERROR: Could not send email. Visualization is accessible at\n\nhttp://127.0.0.1:5000/output/"

        return(output)
    else:
        #clean temporary folder
        if map_file is not None or map_dbfile is not None:
            cmd = 'rm '+TEMP_FOLDER+str(task_id)+'/map_* '
            subprocess.check_output(cmd, shell=True)
        if shape_file is not None or db_file is not None:
            cmd = 'rm '+TEMP_FOLDER+str(task_id)+'/ref_* '
            subprocess.check_output(cmd, shell=True)
        cmd = 'rm '+TEMP_FOLDER+str(task_id)+'/*fa '
        subprocess.check_output(cmd, shell=True)

        #zip folder
        cmd = 'tar -zcf ' + os.path.join(TEMP_FOLDER,(str(task_id)+'.tar.gz')) + ' ' + TEMP_FOLDER + str(task_id)
        subprocess.check_output(cmd, shell=True)

        #check MongoDB
        client = MongoClient(os.environ['OPENSHIFT_MONGODB_DB_URL'])
        db = client['ensemblerna']
        fs = GridFS(db)
        f = gzip.open(os.path.join(TEMP_FOLDER,(str(task_id)+".tar.gz")), 'rb')
        input_gz = fs.put(f, filename=(str(task_id)+".tar.gz"), task_id=task_id)
        f = open(os.path.join(TEMP_FOLDER,task_id+"/interactive.html"), 'rb')
        input_html = fs.put(f, filename="interactive.html", task_id=task_id)
        client.close()

        #delete temporary and upload files
        cmd = 'rm -r ' + os.path.join(TEMP_FOLDER, task_id+'*')
        subprocess.check_output(cmd, shell=True)

        #email output
        if email is not None and email != '' and email != 'None':
                
            #email message
            TEXT = "\n\nYour visualization is ready. To view your output follow the link below. Thank you for using EnsembleRNA. (Please do not reply to this email)\n\n"+"http://ensemblerna-cbtolson.rhcloud.com/output/"+str(task_id)+"\n\n\n\n"
    
            #send email
            gmail_user = "ensemblerna@gmail.com"
            gmail_pwd = "ribosnitch"
            FROM = "ensemblerna@gmail.com"
            TO = [email]
            SUBJECT = "EnsembleRNA Output is Ready"
            message = """\From: %s\nTo: %s\nSubject: %s\n\n%s """ % (FROM, ", ".join(TO), SUBJECT, TEXT)
            try:
                server = smtplib.SMTP("smtp.gmail.com", 587)
                server.ehlo()
                server.starttls()
                server.login(gmail_user, gmail_pwd)
                server.sendmail(FROM, TO, message)
                server.close()
            except:
                return "Could not send email. Visualization is accessible at\n\nhttp://http://ensemblerna-cbtolson.rhcloud.com/output/"+str(task_id)

        return(True)
