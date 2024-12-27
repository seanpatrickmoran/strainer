#import json
#
#class C3Image3:
#    def __init__(self,name,coordinates,numpyarr,viewing_vmax,dimensions,resolution,hic_path,PUB_ID,**kwargs):
#        c1,x1,x2,c2,y1,y2=coordinates
#        _narr = numpyarr.astype('int32')
##         _narr = numpyarr.flatten().tolist()
#        if not c1.startswith("chr"):
#            c1="chr"+c1
#            c2="chr"+c2
#        self.coordinates=f"{c1},{x1},{x2},{c2},{y1},{y2}"
#        self.entity = {"name": name,
#                       "coordinates":self.coordinates,
#                       "numpyarr":_narr.flatten().tobytes(),
#                       "viewing_vmax": viewing_vmax,
#                       "dimensions": dimensions,
#                       "resolution": resolution,
#                       "hic_path": hic_path,
#                       "PUB_ID": PUB_ID,
#                       "meta":json.dumps(kwargs)
#                      }
#        
#    def __repr__(self):
##         return f"c3Image({self.c1}:{self.x1}-{self.x2}, {self.c2}:{self.y1}-{self.y2})"
#        val = self.entity["name"]
##         print('huh')
#        return f"{val}"

import json


class C3Image4:
    def __init__(self,name,coordinates,numpyarr,viewing_vmax,dimensions,resolution,hic_path,PUB_ID,**kwargs):
        c1,x1,x2,c2,y1,y2=coordinates
        _narr = numpyarr.astype('int32')
        if not c1.startswith("chr"):
            c1="chr"+c1
            c2="chr"+c2
        self.coordinates=f"{c1},{x1},{x2},{c2},{y1},{y2}"
        self.entity = {"name": name,
                       "dataset": kwargs.get("dataset",""),
                       "condition": kwargs.get("condition",""),
                       "coordinates":self.coordinates,
#                       maybe we convert the numpyarr into 4 channel images first?
                       "numpyarr":_narr.flatten().tobytes(),
                       "viewing_vmax": viewing_vmax,
                       "dimensions": dimensions,
                       "resolution": resolution,
                       "hic_path": hic_path,
                       "PUB_ID": PUB_ID,
                       "norm": kwargs.get("norm",""),
                       "meta":json.dumps(kwargs)
                      }
        
    def __repr__(self):
        val = self.entity["name"]
        return f"{val}"
        

 def extractor(dictionary, *args):
     return {value:dictionary[value] for value in args}

