##########################################################################################################
#Functions for viewing html templates
##########################################################################################################
import subprocess
import os
from flask import Flask, render_template, request, redirect, url_for, send_from_directory, make_response, g
from flask_util_js import FlaskUtilJs
from validate_email import validate_email
from werkzeug import secure_filename
from pymongo import MongoClient
from gridfs import *
import gzip
from rnavis.forms import RnaForm
from rnavis.checkError import *
from rnavis.tasks import runEnsembleRNA

#create Flask app
app = Flask(__name__)
fujs = FlaskUtilJs(app)

#initialize variables
app.config['SCRIPT_FOLDER'] = os.path.join(os.environ['OPENSHIFT_REPO_DIR'],'wsgi/static/')
app.config['STATIC_FOLDER'] = '/static/'
app.config['UPLOAD_FOLDER'] = os.environ['OPENSHIFT_TMP_DIR']
app.config['TEMP_FOLDER'] = os.environ['OPENSHIFT_TMP_DIR']

#set up callback requests
def after_this_request(func):
    if not hasattr(g, 'call_after_request'):
        g.call_after_request = []
    g.call_after_request.append(func)
    return func

#after request callback
@app.after_request
def per_request_callbacks(response):
    for func in getattr(g, 'call_after_request', ()):
        response = func(response)
    return response

#get interactive
@app.route('/outhtml')
def outhtml(task_id=None):
    
    #restrict access
    if request.referrer is None:
        return render_template("403.html")
    if('/output' not in request.referrer):
        return render_template("403.html")

    #initialize variable
    task_id = request.args.get('task_id')
    
    #retrieve from MongoDB
    client = MongoClient(os.environ['OPENSHIFT_MONGODB_DB_URL'])
    db = client.ensemblerna
    fs = GridFS(db)
    fileobj = fs.get_last_version(filename="interactive.html", task_id=task_id)
    response = app.response_class(fileobj, mimetype="text/html", direct_passthrough=True)
    client.close()

    #get interactive
    return response

#get zip
@app.route('/ensemble_output')
def ensemble_output(task_id=None):
    
    #restrict access
    if request.referrer is None:
        return render_template("403.html")
    if('/output' not in request.referrer):
        return render_template("403.html")
    
    #initialize variable
    task_id = request.args.get('task_id')
    
    #retrieve from MongoDB
    client = MongoClient(os.environ['OPENSHIFT_MONGODB_DB_URL'])
    db = client['ensemblerna']
    fs = GridFS(db)
    fileobj = fs.get_last_version(filename=(str(task_id)+".tar.gz"), task_id=task_id)
    client.close()

    #write file
    with gzip.open(app.config['TEMP_FOLDER']+(str(task_id)+'.tar.gz'), 'wb') as f:
        f.write(fileobj.read())

    @after_this_request
    def remove_file(response):
        os.remove(app.config['TEMP_FOLDER']+(str(task_id)+'.tar.gz'))
        return(response)

    #get zip
    return send_from_directory(app.config['TEMP_FOLDER'], (str(task_id)+'.tar.gz'))


#display output
@app.route('/output/<string:task_id>')
def output(task_id=None):
    return render_template("output.html", task_id=task_id)


#display errorpage
@app.route('/errorpage')
def errorpage(error=None, task_id=None):
    
    #restrict access
    if request.referrer is None:
        return render_template("403.html")
    if('/inputs' not in request.referrer and '/loading' not in request.referrer):
        return render_template("403.html")

    #initialize variables
    task_id = request.args.get('task_id')
    error = request.args.get('error')

    #delete temporary and upload files
    if(os.path.exists(os.path.join(app.config['TEMP_FOLDER'], task_id)) == True):
        cmd = 'rm -r ' + os.path.join(app.config['TEMP_FOLDER'], task_id+'*')
        subprocess.check_output(cmd, shell=True)

    #get error message
    return render_template("errorpage.html", error=error)


#run ensemble rna while loading page displayed
@app.route('/running', methods=('GET', 'POST'))
def running(tf=None, task_id=None):
    
    #restrict access
    if request.referrer is None:
        return render_template("403.html")
    if('/loading' not in request.referrer):
        return render_template("403.html")

    #GET request
    if request.method == 'GET':
        #initialize variables
        tf = request.args.get('tf')
        task_id = request.args.get('task_id')

        #go to output page
        if tf == 'True':
            return redirect(url_for("output", task_id=task_id))
        
        #go to error page
        if 'does not exist' in tf:
            return redirect(url_for("errorpage", error="Visualization not found. Do not navigate away from loading page unless you have supplied an email address.", task_id=task_id))
        return redirect(url_for("errorpage", error=tf, task_id=task_id))

    #POST request
    return render_template("403.html")


#display loading
@app.route('/loading/<string:task_id>/')
def loading(task_id='startload', ref_fasta=None, ref_shape=None, ref_db=None, map_fasta=None, map_db=None, rgstart=None, rgend=None, size=None, maxd=None, upload_id=None, email=None):
    
    #initialize variables
    ref_fasta = request.args.get('ref_fasta')
    ref_shape = request.args.get('ref_shape')
    ref_db = request.args.get('ref_db')
    map_fasta = request.args.get('map_fasta')
    map_db = request.args.get('map_db')
    rgstart = request.args.get('rgstart')
    rgend = request.args.get('rgend')
    size = request.args.get('size')
    maxd = request.args.get('maxd')
    upload_id = request.args.get('upload_id')
    email = request.args.get('email')
    
    #restrict access
    if request.referrer is None:
        return render_template("403.html")
    if('/inputs' not in request.referrer and '/loading' not in request.referrer):
        return render_template("403.html")

    #get task id
    task_id = upload_id

    #run ensemblerna
    out = runEnsembleRNA.delay(task_id, ref_fasta, ref_shape, ref_db, map_fasta, map_db, rgstart, rgend, size, maxd, email, app.config['TEMP_FOLDER'], app.config['SCRIPT_FOLDER'])
    
    #go to loading page
    return render_template("loading.html", task_id=task_id, celery_id=out.id)

#display inputs
@app.route('/inputs', methods=('GET', 'POST'))
def inputs():
    form = RnaForm()
    error = None
    if request.method == 'POST':
        
        #get upload id
        upload_id = getTaskID(app)
        
        #check required file
        ref_fasta = request.files['ref_fasta']
        if ref_fasta.filename =='':
            return render_template("inputs.html", form=form, error="ERROR: Reference FASTA required.")
        
        #make upload folder
        cmd = 'mkdir ' + os.path.join(app.config['UPLOAD_FOLDER'], upload_id)
        subprocess.check_output(cmd, shell=True)
        
        #check reference fasta
        fn = os.path.join(upload_id, secure_filename(ref_fasta.filename))
        ref_fasta.save(os.path.join(app.config['UPLOAD_FOLDER'], fn))
        ref_fasta = checkFasta(os.path.join(app.config['UPLOAD_FOLDER'], fn))
        if ref_fasta['tf'] != True:
            return render_template("inputs.html", form=form, error=ref_fasta['tf'])
        length = ref_fasta['length']
        ref_fasta = fn
        
        #check reference shape
        ref_shape = request.files['ref_shape']
        if ref_shape.filename =='':
            ref_shape = None
        else:
            fn = os.path.join(upload_id, 'ref_'+secure_filename(ref_shape.filename))
            ref_shape.save(os.path.join(app.config['UPLOAD_FOLDER'], fn))
            tf = checkSHAPE(os.path.join(app.config['UPLOAD_FOLDER'], fn))
            if tf != True:
                return render_template("inputs.html", form=form, error=tf)
            ref_shape = fn
        
        #check reference db
        ref_db = request.files['ref_db']
        if ref_db.filename =='':
            ref_db = None
        else:
            fn = os.path.join(upload_id, 'ref_'+secure_filename(ref_db.filename))
            ref_db.save(os.path.join(app.config['UPLOAD_FOLDER'], fn))
            tf = checkDB(os.path.join(app.config['UPLOAD_FOLDER'], fn), length)
            if tf != True:
                return render_template("inputs.html", form=form, error=tf)
            ref_db = fn
        
        #check map fasta
        map_fasta = request.files['map_fasta']
        if map_fasta.filename =='':
            map_fasta = None
        else:
            fn = os.path.join(upload_id, 'map_'+secure_filename(map_fasta.filename))
            map_fasta.save(os.path.join(app.config['UPLOAD_FOLDER'], fn))
            map_fasta = checkFasta(os.path.join(app.config['UPLOAD_FOLDER'], fn), length)
            if map_fasta['tf'] != True:
                return render_template("inputs.html", form=form, error=map_fasta['tf'])
            map_fasta = fn
        
        #check map db
        map_db = request.files['map_db']
        if map_db.filename =='':
            map_db = None
        else:
            fn = os.path.join(upload_id, 'map_'+secure_filename(map_db.filename))
            map_db.save(os.path.join(app.config['UPLOAD_FOLDER'], fn))
            tf = checkDB(os.path.join(app.config['UPLOAD_FOLDER'], fn), length)
            if tf != True:
                return render_template("inputs.html", form=form, error=tf)
            map_db = fn
        
        #check range
        start = request.form.get('ref_rgstart')
        end = request.form.get('ref_rgend')
        if start == '':
            start = 1
        elif(start.isdigit()==False):
            return render_template("inputs.html", form=form, error="ERROR: Range start must be a positive integer.")
        if end == '':
            end = length
        elif(end.isdigit()==False):
            return render_template("inputs.html", form=form, error="ERROR: Range end must be a positive integer.")
        rg = checkRange([int(start), int(end)], length)
        if isinstance(rg, str):
            return render_template("inputs.html", form=form, error=rg)
        
        #check size
        size = request.form.get('ref_size')
        if size == '':
            size = 10
        elif(size.isdigit()==False):
            return render_template("inputs.html", form=form, error="ERROR: Size must be a positive integer.")
        size = checkSize(int(size), length, rg)
        if isinstance(size, str):
            return render_template("inputs.html", form=form, error=size)
        
        #check maxd
        maxd = request.form.get('ref_maxd')
        if maxd == '':
            maxd = None
        elif(maxd.isdigit()==False):
            return render_template("inputs.html", form=form, error="ERROR: Maximum pairing distance must be a positive integer.")
        else:
            maxd = checkMaxD(int(maxd))
            if isinstance(maxd, str):
                return render_template("inputs.html", form=form, error=maxd)
    
        #check email
        email = request.form.get('email')
        if email == '':
            email = None
        elif(validate_email(email)!=True):
            return render_template("inputs.html", form=form, error="ERROR: Email address must be valid.")
        
        #go to loading page
        return redirect(url_for("loading", task_id='startload', ref_fasta=ref_fasta, ref_shape=ref_shape, ref_db=ref_db, map_fasta=map_fasta, map_db=map_db, rgstart=(rg[0]+1), rgend=(rg[-1]+1), size=size, maxd=maxd, upload_id=upload_id, email=email))
    else:
        return render_template("inputs.html", form=form, error=error)


#display 404
@app.errorhandler(404)
def page_not_found(e):
    return render_template('404.html'), 404


#display 403
@app.errorhandler(403)
def page_not_found(e):
    return render_template('403.html'), 403


#display help
@app.route('/help')
def help():
    return render_template("help.html")


#display example
@app.route('/example')
def example():
    return render_template("example.html")


#display download
@app.route('/download')
def download():
    return render_template("download.html")


#display homepage
@app.route('/')
def homepage():
    return render_template("homepage.html")


##########################################################################################################
#CELERY TASK
##########################################################################################################
#check task
@app.route('/status', methods=('GET', 'POST'))
def status(celery_id=None):
    
    #restrict access
    if request.referrer is None:
        return render_template("403.html")
    if('/inputs' not in request.referrer and '/loading' not in request.referrer):
        return render_template("403.html")

    #POST request
    if request.method == 'POST':
        #initialize variable
        celery_id = request.args.get('celery_id')
        
        #check task
        task = runEnsembleRNA.AsyncResult(str(celery_id))
        if task.ready() == False:
            return 'False'
        else:
            return str(task.get())

    #GET request
    return render_template("403.html")




