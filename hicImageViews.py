import cooltools
import cooler

import numpy as np
import hicstraw

from skimage.transform import hough_line, hough_line_peaks
from skimage.feature import canny
from skimage.filters import threshold_otsu
from skimage.draw import line as draw_line
from skimage import data
import scipy

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import gridspec
from matplotlib import cm
import seaborn as sns

from CCCImageClasses import C3Image4, extractor


REDMAP = LinearSegmentedColormap.from_list("bright_red", [(1,1,1),(1,0,0)])
# helper function for plotting
def plot_hic_map(dense_matrix, maxcolor):
    plt.matshow(dense_matrix, cmap=REDMAP, vmin=0, vmax=maxcolor)
    plt.show()


def target_square_windowing(x1,x2,y1,y2,res):#,width):
    x1,x2,y1,y2,res=int(x1),int(x2),int(y1),int(y2),int(res)
    
    xline=x1+(x2-x1)//2
    yline=y1+(y2-y1)//2
    r1,r2=xline-res*10,xline+10*res
    r3,r4=yline-res*10,yline+10*res
    return r1,r2,r3,r4

def windowing(x1,x2,y1,y2,res,width):
    x1,x2,y1,y2,res=int(x1),int(x2),int(y1),int(y2),int(res)
    target = res*width
    #210,000
    if x2-x1>=target and y2-y1>=target:
        return x1,x2,y1,y2
    
    r1,r2,r3,r4 = x2-(x2-x1)//2-target//2, x2-(x2-x1)//2+target//2, y2-(y2-y1)//2-target//2, y2-(y2-y1)//2+target//2
    if r1<0 or r3<0:
        choose = min(r1,r3)
        r1-=choose
        r3-=choose
        r2+=choose
        r4+=choose
    return r1,r2,r3,r4

def choose_vmax(nimage):
    thresh = threshold_otsu(nimage)
    binary = nimage > thresh
    tested_angles = np.linspace(-np.pi / 2, np.pi / 2, 360, endpoint=False)
    h, theta, d = hough_line(binary, theta=tested_angles)
    
    c_angle,c_dist = 0,0
    chosen = 10
        
    for _,angle,dist in  zip(*hough_line_peaks(h, theta, d)):
        (x0, y0) = dist * np.array([np.cos(angle), np.sin(angle)])
        slope=np.tan(angle + np.pi / 2)
        
        if (abs(1-abs(slope))<abs(1-abs(chosen))):
            c_angle=angle
            c_dist=dist
            chosen=slope
    

    x0, y0 = c_dist * np.array([np.cos(c_angle), np.sin(c_angle)])
    slope = np.tan(c_angle + np.pi / 2)
    n=nimage.shape[0]

#     print(x0, y0, slope)
#     print(int(np.round(y0-x0)), n)
          
    if slope > 10:
          return nimage,np.max(nimage)

    mimage = nimage.copy()
    if int(np.round(y0-x0)) <= n//10:
        #the main diagonal is too close to the center of image.
        mimage = np.triu(mimage,k=n//10+2) 
#         np.median(mimage)+4*np.std(mimage)
        return mimage,np.max(mimage)

    else:
                  # For normalization for visualization, remove image below the expected focal point
        mimage = np.triu(mimage,k=-n//10) 
#         np.median(mimage)+4*np.std(mimage)
##################################################
#         # Old code takes the slope found in hough transform and blanks beneath it.
#         incrementor = 0
#         buffer=n//10+2
#         for coord in range(int(np.round(y0-x0))-buffer,n):
#             print(coord, np.round(y0-x0))
#             for i in range(0,min(incrementor+1,n)):
#     #         for i in range(0,min(incrementor+1,n)):
#                 mimage[coord,i]=0
#             incrementor+=1        
#         if np.max(mimage)>100:
#             return mimage, 30
        return mimage,np.max(mimage)
        # return mimage,np.median(mimage)+4*np.std(mimage)


def collect_numpy_matrices(hic_genericpath, bedpath, norm="NONE", res=10000, width=21, unrestricted=False):
    if hic_genericpath.endswith("hic"):
        print('1')
        return _hic_collect_numpy_matrices(hic_genericpath, bedpath, norm, res, width,unrestricted)
    elif hic_genericpath.endswith("cool"):
        print('2')
        return _cooler_collect_numpy_matrices(hic_genericpath, bedpath, norm, res, width,unrestricted)
    
    
    
def _hic_collect_numpy_matrices(hic_path, bedpath, norm="NONE", res=10000, width=21,unrestricted=False):
    hic = hicstraw.HiCFile(hic_path)
    hold_np = {}
    last_seen = ["chr0"]
    with open(bedpath, "r+") as f:
        f.readline()
        for no,line in enumerate(f):
#             print(line)
            c1,x1,x2,c2,y1,y2 = line.rstrip('\n').split('\t')[:6]
#             if not c1.startswith("chr"):
#                 c1="chr"+c1
#                 c2="chr"+c2
            if c1!=last_seen[0]:
                _matrix_object_ = hic.getMatrixZoomData(c1, c2, "observed", "NONE", "BP", res)
                last_seen.pop()
                last_seen+=[c1]

#             r1,r2,r3,r4 = target_square_windowing(x1,x2,y1,y2,res)                
            r1,r2,r3,r4 = windowing(x1,x2,y1,y2,res,width)
            _np_mat_ = _matrix_object_.getRecordsAsMatrix(r1,r2,r3,r4) 
            if not unrestricted and _np_mat_.shape==(1,1):
                continue
            
            hold_np[no] = {"original coordinates":[c1,x1,x2,c2,y1,y2], "window coordinates": [r1,r2,r3,r4], "numpy_window": _np_mat_, "mumpy_window": _np_mat_, "viewing_vmax":0}
            hold_np[no]["mumpy_window"],hold_np[no]["viewing_vmax"]=choose_vmax(_np_mat_)
        return hold_np   
    
    
    

def _cooler_collect_numpy_matrices(cooler_path, bedpath, norm="NONE", res=10000, width=21,unrestricted=False):
    #wip broken...
    hold_np = {}
    last_seen = ["chr0"]

    def coolload(filepath, resolution):
#         if filepath.endswith("mcool"):
        filepath+=f"::/resolutions/{resolution}"
        cool = cooler.Cooler(filepath)
        print(f"loaded {filepath}")
        return cool    
    
    cic = coolload(cooler_path,res)
    with open(bedpath, "r+") as f:
        f.readline()
        for no,line in enumerate(f):
            try:
                c1,x1,x2,c2,y1,y2 = line.rstrip('\n').split('\t')[:6]
            except ValueError:
                c1,x1,y1 = line.rstrip('\n').split('\t')[:3]
                x1=int(float(x1))
                y1=int(float(y1))
                x2=x1+res
                c2=c1
                y2=y1+res
#             print(c1,x1,x2,c2,y1,y2)
            if not c1.startswith("chr"):
                c1="chr"+c1
                c2="chr"+c2
            if c1!=last_seen[0]:
#                 print("here")
#                 _table_obj_ = cic.matrix(balance=False, as_pixels=True, join=True).fetch((c1,0,cic.chromsizes[c1]))
#                 print((c1,0,cic.chromsizes[c1]),(c2,0,cic.chromsizes[c2]))
                _mat_obj_ = cic.matrix().fetch(c1,c2)
#                     _mat_obj_ = cic.matrix(balance=False).fetch(f"chr{c1}",f"chr{c2}")
#                 _matrix_object_ = hic.getMatrixZoomData(c1, c2, "observed", "NONE", "BP", res)
                last_seen.pop()
                last_seen+=[c1]

#             r1,r2,r3,r4 = target_square_windowing(x1,x2,y1,y2,res)                        
            r1,r2,r3,r4 = windowing(x1,x2,y1,y2,res,width)
            _np_mat_ = _mat_obj_[r1//res:r2//res,r3//res:r4//res]
            if not unrestricted and _np_mat_.shape==(1,1):
                continue
            
            _np_mat_[np.isnan(_np_mat_)]=0
  
            hold_np[no] = {"original coordinates":[c1,x1,x2,c2,y1,y2], "window coordinates": [r1,r2,r3,r4], "numpy_window": _np_mat_, "mumpy_window": _np_mat_, "viewing_vmax":0}
#             vmax_flagged = np.max(_np_mat_)
#             if vmax_flagged<100:
#             if vmax_flagged<3:
#                 hold_np[no]["viewing_vmax"] = vmax_flagged
#             else:
            hold_np[no]["mumpy_window"],hold_np[no]["viewing_vmax"]=choose_vmax(_np_mat_)
        return hold_np    
    

def check_select(i,image_dict,squeeze):
    flags={0:False,1:False}
    nsize=image_dict[i]['numpy_window'].shape[0]
#     print(image_dict[i]['numpy_window'][nsize//2-squeeze:nsize//2+squeeze,nsize//2-squeeze:nsize//2+squeeze].shape)
    for oper in [0,1]:
        oper_peaks = scipy.signal.argrelextrema(np.sum(image_dict[i]['numpy_window'][nsize//2-squeeze:nsize//2+squeeze,nsize//2-squeeze:nsize//2+squeeze],axis=oper),np.greater)[0]
        if not oper_peaks.any():
            return False
#         print(oper_peaks)
        for peak in oper_peaks:
#             array_size = nsize//2+squeeze-nsize//2-squeeze
            if nsize//2-squeeze//2 <= peak+(nsize//2-squeeze)<=nsize//2+squeeze//2:
                flags[oper]=True
                continue
                
    return flags[0]*flags[1]


def subsample(i,image_dict,squeeze):
    nsize=image_dict[i]['numpy_window'].shape[0]
    return image_dict[i]['numpy_window'][nsize-squeeze:nsize+squeeze,nsize-squeeze:nsize+squeeze]

def only_populated_windows(numpy_hic_dictionary,window_scan_size,choose_scaler=4):
    mn = numpy_hic_dictionary[0]['numpy_window'].shape[0]
    choose_arr = []
    for i in list(numpy_hic_dictionary.keys()):
        if np.count_nonzero(numpy_hic_dictionary[i]['numpy_window'])>mn**2/choose_scaler: #and check_select(i,numpy_hic_dictionary,window_scan_size):
            choose_arr += [i]

    print(len(choose_arr))
    return choose_arr

def vis_block(npdictin, enum_index_array, offset, squeeze=40):
#     offset = 0

    plt.figure(0,figsize=(16,16))
    plots = []
    incr=0
    for i in range(4):
        for j in range(4):
            ax = plt.subplot2grid((4,4), (i,j))
            ax.set_title(f"Grand_canyon: {enum_index_array[incr+offset]}")
#             ax.matshow(subsample(enum_index_array[incr+offset],npdictin,40), cmap=REDMAP, vmin=0, vmax=npdictin[enum_index_array[incr+offset]]["viewing_vmax"])
            ax.matshow(subsample(enum_index_array[incr+offset],npdictin, squeeze), cmap=REDMAP, vmin=0, vmax=npdictin[enum_index_array[incr+offset]]["viewing_vmax"])
#             ax.matshow(npdictin[incr+offset]['numpy_window'], cmap=REDMAP, vmin=0, vmax=npdictin[enum_index_array[incr+offset]]["viewing_vmax"])
            incr+=1

    plt.subplots_adjust(hspace=.4)
    plt.show()
    
def vis_block_10(npdictin, enum_index_array, offset, squeeze=40):
#     offset = 0

    plt.figure(0,figsize=(16,16))
    plots = []
    incr=0
    for i in range(10):
        for j in range(10):
            ax = plt.subplot2grid((10,10), (i,j))
            ax.set_title(f"{enum_index_array[incr+offset]}")
#             ax.matshow(subsample(enum_index_array[incr+offset],npdictin,40), cmap=REDMAP, vmin=0, vmax=npdictin[enum_index_array[incr+offset]]["viewing_vmax"])
            ax.matshow(subsample(enum_index_array[incr+offset],npdictin, squeeze), cmap=REDMAP, vmin=0, vmax=npdictin[enum_index_array[incr+offset]]["viewing_vmax"])
#             ax.matshow(npdictin[incr+offset]['numpy_window'], cmap=REDMAP, vmin=0, vmax=npdictin[enum_index_array[incr+offset]]["viewing_vmax"])
            incr+=1
    plt.tight_layout()
    plt.subplots_adjust(hspace=.4)
    plt.show()
