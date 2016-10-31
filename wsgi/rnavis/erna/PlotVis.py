##########################################################################################################
#Functions for plotting structures
    #plotMap
    #plotRef
    #writePlotMDS
##########################################################################################################
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import pylab
import matplotlib.gridspec as gridspec
import numpy as np
import mpld3
from mpld3 import plugins, utils
import math
import subprocess
import csv
import matplotlib.image as mpimg
from sklearn.manifold import MDS
from sklearn.metrics import euclidean_distances
from scipy.spatial.distance import pdist
from scipy.spatial import distance
import warnings
from collections import Counter

##########################################################################################################
#Function plotting structures for map of conformational space
#Input: 2d matrix, frequency, pattern sequence, representative, outfile header
#Output: csv file, db file, pdf file, png file
##########################################################################################################
def plotMap(maparr, freq, nest, seqs, dbfile, map2d, outfile, plotm='T'):

    #mutli-dimensional scaling
    similarities = euclidean_distances(np.matrix(maparr))
    mds = MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=np.random.RandomState(seed=3), dissimilarity="precomputed", n_jobs=1)
    pos = mds.fit(similarities).embedding_

    #plot attributes
    N = len(pos)
    #size = [20*n for n in freq]
    size = 8000
    color = np.array(range(N))
    
    if str(plotm) == 'T':
    
        #plot MDS
        fig, ax = plt.subplots(figsize=(10,10))
        warnings.filterwarnings("ignore")
        scatter = ax.scatter(np.array(pos[:,0]), np.array(pos[:,1]), c=color, s=size, alpha=0.3, cmap=plt.cm.viridis, marker='s')
        plt.xlabel('Dimension 1', fontsize=20, labelpad=20)
        plt.ylabel('Dimension 2', fontsize=20, labelpad=20)
        #plt.axis([xmin, xmax, ymin, ymax])
        plt.tick_params(labelsize=15, length=14, direction='out', pad=15, top='off', right='off')

        #save figures
        fig.savefig(outfile + '.png', bbox_inches='tight', format='png')
        fig.savefig(outfile + '.pdf', bbox_inches='tight', format='pdf')
        plt.close(fig)
        warnings.resetwarnings()
        
        #write csv file
        writePlotMDS(freq, nest, seqs, dbfile, pos, maparr, map2d, outfile)

    return pos

##########################################################################################################
#Function plotting structures on bar plot
#Input: nest frequencies, nest pattern, nest patterns for each structure, dot bracket structures, maximum number of clusters, outfile header
#Output: png file, pdf file
##########################################################################################################
def plotRef(reffreq, refnest, refarr, mappos, mapnest, mapseqs, mapdb, maparr, map2d, outfile, refdb, refseqs, rg):
    
    #initialize variables
    freq = [0 for x in range(len(mapnest))]
    reffreq = list(reffreq)
    
    #reshape arrays
    diff = len(maparr[0])-len(refarr[0])
    if(diff>0):
        for k in range(diff):
            refarr = np.c_[refarr, np.zeros(len(refarr))]
            refseqs = [x+",0" for x in refseqs]
    if(diff<0):
        for k in range(abs(diff)):
            maparr = np.c_[maparr, np.zeros(len(maparr))]

    #get frequency
    maparr = np.array(maparr)
    for i in range(len(refarr)):
        maxdst = float("inf")
        r = list(refarr[i])
        if np.any(np.all((np.array(r)-np.array(maparr))==0, axis=1))==True:
            ind = int(list(np.where(np.all((np.array(r)-np.array(maparr))==0, axis=1))[0])[0])
            freq[ind] = reffreq[i]+freq[ind]
        else:
            for j in range(len(maparr)):
                m = list(maparr[j])
                dst = distance.euclidean(np.array(r),np.array(m))
                if(maxdst>dst):
                    mn = j
                    maxdst = dst
            freq[mn] = freq[mn]+reffreq[i]

    #plot attributes
    pos = mappos
    N = len(pos)
    size = [20*n for n in freq]
    color = np.array(range(N))

    #plot MDS
    fig, ax = plt.subplots(figsize=(10,10))
    warnings.filterwarnings("ignore")
    scatter = ax.scatter(np.array(pos[:,0]), np.array(pos[:,1]), c=color, s=size, alpha=0.3, cmap=plt.cm.viridis)
    plt.xlabel('Dimension 1', fontsize=20, labelpad=20)
    plt.ylabel('Dimension 2', fontsize=20, labelpad=20)
    #plt.axis([xmin, xmax, ymin, ymax])
    plt.tick_params(labelsize=15, length=14, direction='out', pad=15, top='off', right='off')
    warnings.resetwarnings()

    #save figures
    fig.savefig(outfile + '.png', bbox_inches='tight', format='png')
    fig.savefig(outfile + '.pdf', bbox_inches='tight', format='pdf')
    plt.close(fig)

    #write csv file
    ref = writePlotMDS(freq, mapnest, mapseqs, mapdb, mappos, maparr, map2d, outfile, refdb, refseqs, rg)

    #return output
    structure = ref['structs']
    diversity = ref['diversity']

    return{'structs':structure, 'diversity':diversity, 'freq':freq}

##########################################################################################################
#Function outputing csv file for scatter plot
#Input:  cluster number, vector representation, frequency, sequence, representative, outfile header
#Output: csv file
##########################################################################################################
def writePlotMDS(num, nest, seqs, dbfile, mappos, maparr, map2d, outfile, refdb=None, refseqs=None, rg=None):
    
    #initialize variables
    clusters = range(1,len(num)+1)
    frequency = list(num)
    
    #loop through clusters
    structure = [0 for i in range(len(num))]
    diversity = [[] for i in range(len(num))]
    for i in range(len(num)):
        indices = [j for j, x in enumerate(seqs) if x == nest[i]]
        db = [dbfile[j] for j in indices]
        
        #get cluster structure medoids
        structs = [j.replace('.','0') for j in db]
        structs = [j.replace('(','1') for j in structs]
        structs = [j.replace(')','1') for j in structs]
        structs = [[int(x) for x in list(j)] for j in structs]
        dst = pdist(1-np.matrix(structs),'jaccard')
        dst = np.sum(dst, axis=0)
        ind = np.argmin(dst)
        structure[i] = db[ind]

        #get diversity
        if refdb is not None:
            indices = [j for j, x in enumerate(refseqs) if x == nest[i]]
            db = [refdb[j] for j in indices]
            db = [x[rg[0]:rg[-1]] for x in db]
            structs = [j.replace('.','0') for j in db]
            structs = [j.replace('(','1') for j in structs]
            structs = [j.replace(')','1') for j in structs]
            structs = [[int(x) for x in list(j)] for j in structs]
            if not indices:
                diversity[i] = [0, 0]
            else:
                d = Counter(db)
                d = sorted(d.items())
                n = [x[0] for x in d] #unique structures
                m = [x[1] for x in d] #frequency
                divsz = 1
                if len(m) > 1:
                    if len(structs) < 2:
                        divsz = 1
                    if len(structs) < 3:
                        divsz = pdist(1-np.matrix(structs),'jaccard').tolist()[0]
                    else:
                        divsz = min(np.diag(np.matrix(pdist(1-np.matrix(structs),'jaccard')),k=1))
                divfreq = max(m)/len(db)
                diversity[i] = [divsz, divfreq]
       
    #write to file
    with open(outfile+'.csv', 'w') as csvfile:
        fieldnames = ['cluster', 'xy-coords', 'frequency', 'mediod-structure', 'vectorization']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for i in range(len(num)):
            basic_dict = {'cluster': clusters[i], 'xy-coords': np.array_str(mappos[i]), 'frequency': frequency[i], 'mediod-structure': structure[i], 'vectorization': maparr[i],}
            writer.writerow(basic_dict)

    return{'structs':structure, 'diversity':diversity}

#########################################
#Read dot-bracket to helix
#########################################
def _readDB(symbolString):
    #initialize variables
    s = []
    left = []
    right = []
    #get positions
    for index in range(len(symbolString)):
        symbol = symbolString[index]
        if symbol == "(":
            s.append(index)
        if symbol == ")":
            right.append(index)
            left.append(s[-1])
            s = s[:-1]
    #organize output
    output = [[left[x], right[x]] for x in range(len(left))]
    output = sorted(output, key=lambda x: x[0])
    return(output)

#########################################
#Get arc
#########################################
def _anglesArc(v, theta):
    #set theta.0X
    if v[1] >= 0:
        thetaOX = math.acos(v[0])
    else:
        thetaOX = 2 * math.pi - math.acos(v[0])
    #get angle
    angs = [thetaOX - theta, thetaOX + theta]
    return(angs)

#########################################
#Get arc
#########################################
def _arc(c, r, v, theta):
    #initialize variables
    angles = _anglesArc(v, theta)
    seqang = np.linspace(angles[0], angles[1], 100)
    xdata = [r*math.cos(x)+c[0]+1 for x in seqang]
    ydata = [r*math.sin(x)+c[1] for x in seqang]
    #plot arc
    return{'xdata':xdata, 'ydata':ydata}

#########################################
#Plot arc
#########################################
def _plotArc(i, j, y):
    #initalize variables
    center = [np.mean([i, j]), y]
    radius = abs(i-j)/2
    vector = [0, 1]
    theta = math.pi/2
    #get arc
    output = _arc(center, radius, vector, theta)
    return{'xdata':output['xdata'], 'ydata':output['ydata']}

#########################################
#Plot helix
#########################################
def _plotHelix(helix):
    #initialize variables
    i = [x for x,y in helix]
    j = [y for x,y in helix]
    xdata = [[] for i in range(len(helix))]
    ydata = [[] for i in range(len(helix))]
    #reversed so top helices appear on top of lower helices
    for n in range(len(helix))[::-1] :
        output = _plotArc(i[n], j[n], 0)
        xdata[n] = output['xdata']
        ydata[n] = output['ydata']
    return{'xdata':xdata, 'ydata':ydata}

#########################################
#Plot Interactive
#########################################
def plotInteractive(dir, header, pos, freq,  structs, diversity, vec, rg):
    class ClickInfo(plugins.PluginBase):
        """Plugin for getting info on click"""
        
        JAVASCRIPT = """
            mpld3.register_plugin("clickinfo", ClickInfo);
            ClickInfo.prototype = Object.create(mpld3.Plugin.prototype);
            ClickInfo.prototype.constructor = ClickInfo;
            ClickInfo.prototype.requiredProps = ["idpts", "idline", "data", "idline2", "data2"];
            function ClickInfo(fig, props){
            mpld3.Plugin.call(this, fig, props);
            };
            
            ClickInfo.prototype.draw = function(){
            var pts = mpld3.get_element(this.props.idpts);
            var line = mpld3.get_element(this.props.idline);
            var data = this.props.data;
            var line2 = mpld3.get_element(this.props.idline2);
            var data2 = this.props.data2;
            
            pts.elements().on("mousedown", function(d, i){
            line.data = data[i];
            c = line.data[0]
            for (j = 1; j < data[i].length; j++){
            c = c.concat(line.data[j])
            }
            line.elements().transition()
            .attr("d", line.datafunc(c))
            .style("stroke", this.style.fill);
            
            line.data2 = data2[i];
            line2.elements().transition()
            .attr("d", line2.datafunc(line.data2))
            .style("stroke", this.style.fill);
            });
            }
            """
        def __init__(self, points, line, linedata, line2, linedata2):
            self.dict_ = {"type": "clickinfo",
                "idpts": utils.get_id(points),
                "idline": utils.get_id(line),
                "data": linedata,
                "idline2": utils.get_id(line2),
                "data2": linedata2}
    
    #plot attributes
    N = len(pos)
    size = [20*n for n in freq]
    color = np.array(range(N))
    
    #plot figure
    #fig, ax = plt.subplots(3, figsize=(10,30))
    fig = plt.figure(figsize=(15,7.5))
    gs = gridspec.GridSpec(2, 2)
    ax = [plt.subplot(gs[:,0]), plt.subplot(gs[0,1]), plt.subplot(gs[1,1])]
    
    #plot reference scatter
    warnings.filterwarnings("ignore")
    scatter = ax[0].scatter(np.array(pos[:,0]), np.array(pos[:,1]), c=color, s=size, alpha=0.3, cmap=plt.cm.viridis)
    #ax[0].set_title('Reference Visualization', fontsize=20, y=1.2)
    ax[0].set_xlabel('Dimension 1', fontsize=17)
    ax[0].set_ylabel('Dimension 2', fontsize=17)
    ax[0].tick_params(labelsize=15, length=14, direction='out', pad=15, top='off', right='off')
    
    #order data
    mx = 0
    sx = np.array([0.05,0.1,0.1,0.1,0.1,0.05,0.05,0.05])
    sy = np.array([0.05,0.05,0.05,0.1,0.1,0.1,0.1,0.05])
    dd = [[] for x in range(len(diversity))]
    xdata = [[] for x in range(len(structs))]
    ydata = [[] for x in range(len(structs))]
    linedata = [[] for x in range(len(structs))]
    for i in range(len(diversity)):
        #diversity
        dx = sx+(diversity[i][0]-0.075)
        dy = sy+(diversity[i][1]-0.075)
        dd[i] = np.matrix.tolist(np.column_stack((dx, dy)))
        #medoid
        helix = _readDB(structs[i])
        output = _plotHelix(helix)
        xdata[i] = output['xdata']
        ydata[i] = output['ydata']
        if len(pos) == 1 and vec[0] == [0]:
            mx = None
        elif mx < max(sum(ydata[i],[])):
            mx = max(sum(ydata[i],[]))
        data = [[] for x in range(len(xdata[i]))]
        for j in range(len(xdata[i])):
            d = np.column_stack((np.array(xdata[i][j]), np.array(ydata[i][j]), np.array([i for x in xdata[i][j]])))
            d = np.matrix.tolist(d)
            data[j] = d
        linedata[i] = data
    linedata2 = dd #diversity
    
    #plot mediod
    ax[1].set_xlabel('Nucleotides', fontsize=17)
    ax[1].set_ylabel('', fontsize=17)
    #ax[1].set_title('Mediod Structure', fontsize=20, y=1.2)
    ax[1].axes.get_yaxis().set_ticks([])
    lines = ax[1].plot(np.linspace(0, 10, 100), np.linspace(0, 10, 100), '-w', lw=3, alpha=0.5)
    xdata = [[] for x in range(len(structs))]
    ydata = [[] for x in range(len(structs))]
    width = len(structs[0])
    rd = int(math.ceil(width/10.0))*10+1
    y1,y2 = ax[1].get_ylim()
    x1,x2 = ax[1].get_xlim()
    if mx is not None:
        ax[1].set_aspect("equal")
        ax[1].set_ylim(y1,mx+(width/50)+2)
        ax[1].set_xlim(x1-1.5,rd+0.5)
        ax[1].plot(np.array([rg[0]+1,rg[0]+1]), np.array([y1,mx+(width/50)+2]), color="black", lw=3, alpha=0.5)
        ax[1].plot(np.array([rg[-1]+1,rg[-1]+1]), np.array([y1,mx+(width/50)+2]), color="black", lw=3, alpha=0.5)
        ax[1].annotate(str((rg[0]+1)), xy=(1, 1), xytext=(rg[0]+2, mx+(width/50)-2), color="black", fontsize=15, alpha=0.5)
        ax[1].annotate(str((rg[-1]+1)), xy=(1, 1), xytext=(rg[-1]+2, mx+(width/50)-2), color="black", fontsize=15, alpha=0.5)
    else:
        ax[1].annotate("No Base Pairs Found", xy=(1, 1), xytext=(x2/4, y2/2), color="black", fontsize=40, alpha=0.3)
    
    #plot diversity
    ax[2].set_aspect(0.5)
    ax[2].set_ylim(-0.15,1.15)
    ax[2].set_xlim(-0.1,1.1)
    lines2 = ax[2].plot(np.linspace(0, 0, 100), np.linspace(0, 0, 100), '-w', lw=3, alpha=0.5)
    ax[2].plot(np.array([-1, 2]), np.array([0.5, 0.5]), '--', color="black")
    ax[2].plot(np.array([0.5, 0.5]), np.array([-1, 2]), '--', color="black")
    ax[2].annotate("Mid-H", xy=(1, 1), xytext=(0.14, 1), color="black", fontsize=15)
    ax[2].annotate("High", xy=(1, 1), xytext=(0.15, -0.1), color="black", fontsize=15)
    ax[2].annotate("Mid-L", xy=(1, 1), xytext=(0.74, -0.1), color="black", fontsize=15)
    ax[2].annotate("Low", xy=(1, 1), xytext=(0.75, 1), color="black", fontsize=15)
    ax[2].set_xlabel('Minimum Cluster Correlation', fontsize=17)
    ax[2].set_ylabel('Maximum Cluster Frequency', fontsize=17)
    ax[2].tick_params(labelsize=15, length=14, direction='out', pad=15, top='off', right='off')
    ax[2].set_title(' ', fontsize=10)
    
    #add click plugin
    plugins.connect(fig, ClickInfo(scatter, lines[0], linedata, lines2[0], linedata2))
    
    #add labels plugin
    labels = ['Rep: {0}'.format(x) for x in vec]
    tooltip = mpld3.plugins.PointLabelTooltip(scatter, labels=labels)
    mpld3.plugins.connect(fig, tooltip)
    
    #save figures
    plt.tight_layout()
    mpld3.save_html(fig, dir+'interactive.html')
    warnings.resetwarnings()



