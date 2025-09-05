#testfire

import sqlite3
from typing import List
import struct
import datetime
import numpy as np
import os
import pyBigWig



def serialize_f32(vector: List[float]) -> bytes:
    """serializes a list of floats into a compact "raw bytes" format"""
    return struct.pack("%sf" % len(vector), *vector)




def call(PATH,TIMEOUT):

    connection = sqlite3.connect(PATH, timeout=TIMEOUT)  # Set timeout to 10 seconds
    cursor = connection.cursor()
    return connection,cursor


def untouch(dbPATH,timeout):
    db = sqlite3.connect(dbPATH)
    try:
        db.execute("DROP TABLE imag")#, (kwargs['tablename'],))
        print(f"success")
    except Exception as e:
        print("connection failure")
        print(e)
    finally:
        db.commit()
        db.close()


def _createTable(dbPATH, timeout,**kwargs):

    connection,cursor=call(dbPATH,timeout)
    try:

        cursor = connection.execute("CREATE TABLE imag(key_id, name, dataset, condition, coordinates, numpyarr, viewing_vmax, true_max, hist_rel, hist_true, dimensions, hic_path, PUB_ID, resolution, norm, toolsource, featuretype, epigenomicFactors, motifDirection, meta)")
        print("table make; success")

    finally:
        cursor.close()
        connection.close()


    # sqlite_vec.load(db)
def _readSOURCE_writeVECTOR(dbPATH1, dbPATH2,timeout,**kwargs):
# def _readMatchAllTable(dbPATH,timeout,**kwargs):
    dmapped={
        "+":">",
        "-":"<"
    }
    def binary_search_motif(chromosome, target):
        cArray = kwargs["binary_search_CTCF"][chromosome]
        res=0
        p0=0
        p1=len(cArray)-1
        while p0<=p1:
            mid = p0 + (p1-p0)//2
            if abs(int(cArray[mid][2]) - target) < abs(int(cArray[res][2]) - target):
                res = mid
            elif abs(int(cArray[mid][2]) - target) == abs(int(cArray[res][2]) - target):
                if int(cArray[mid][2])>int(cArray[res][2]):
                    res=mid
                
            if int(cArray[mid][2]) == target:
                return mid
            elif int(cArray[mid][2]) < target:
                p0 = mid + 1
            elif int(cArray[mid][2]) > target:
                p1 = mid - 1
        
        return res


    def _readDB(offset, limit):
        connection_s,cursor_s=call(dbPATH1,timeout)
        connection_t,cursor_t=call(dbPATH2,timeout)

        try:
            cursor_s.row_factory = sqlite3.Row
            cursor_s.execute("SELECT * FROM imag LIMIT ? OFFSET ?", (limit,offset))
#             data = []
            for en in cursor_s.fetchall():
                key_id, name, dataset, condition, coordinates, numpyarr, viewing_vmax, true_max, hist_rel, hist_true, dimensions, hic_path, PUB_ID, resolution, norm, toolsource, featuretype, meta = en
                c1,_x1,_x2,c2,_y1,_y2=coordinates.split(",")
                x1=int(_x1)
                x2=int(_x2)
                y1=int(_y1)
                y2=int(_y2)
                
                
                factors = ["CTCF", 
                   "H3K27ac", 
                   "H3K27me3", 
                   "H3K9me3",
                   "H2AFZ",
                   "H3K36me3",
                   "H3K4me1",
                   "H3K9ac"]

                #use nBins
                nBins = 5
                epigenomicFactors = np.array([0]*nBins*2*len(factors), dtype=float)
                
                #get the file from inverseTableCell
                _cellName = kwargs["inverseTableCell"][dataset]                
                for i in range(len(factors)):
                    """ nBins = 5 example
                     <X-bins> | <Y-bins>
                    [0,0,0,0,0,0,0,0,0,0, <CTCF>
                     0,0,0,0,0,0,0,0,0,0, <H3K27ac>
                     0,0,0,0,0,0,0,0,0,0, <H3K27me3>
                     0,0,0,0,0,0,0,0,0,0, <H3K9me3>
                     0,0,0,0,0,0,0,0,0,0, <H2AFZ>
                     0,0,0,0,0,0,0,0,0,0, <H3K36me3>
                     0,0,0,0,0,0,0,0,0,0, <H3K4me1>
                     0,0,0,0,0,0,0,0,0,0] <H3K9ac>
                    """
                    if kwargs["registerPath"][_cellName][factors[i]]=="":
                        epigenomicFactors[i*2] = 0
                        epigenomicFactors[i*2+nBins] = 0
                        continue
                    bw = pyBigWig.open(kwargs["registerPath"][_cellName][factors[i]])
                    _midpoint = x1+(x2-x1)//2
                    # store1 = bw.stats(c1, x1, x2, nBins=nBins)
                    store1 = bw.stats(c1, _midpoint-int(resolution)//2-int(resolution)*2, _midpoint+int(resolution)//2+int(resolution)*2, nBins=nBins)
                    _midpoint = y1+(y2-y1)//2
                    # store2 = bw.stats(c2, y1, y2)
                    store2 = bw.stats(c2, _midpoint-int(resolution)//2-int(resolution)*2, _midpoint+int(resolution)//2+int(resolution)*2, nBins=nBins)
                    if store1 != []:
               #          for j in range(len(store1)):
                        for j in range(nBins):
                            epigenomicFactors[i*(2*nBins)+j] = store1[j]
                        # epigenomicFactors[i*2] = np.mean(store1)
                    if store2 != []:
                        for j in range(nBins):
                            epigenomicFactors[i*(2*nBins)+j+nBins] = store2[j]
                        # epigenomicFactors[i*2+1] = np.mean(store2)


                # key_id, name, dataset, condition, coordinates, numpyarr, viewing_vmax, true_max, hist_rel, hist_true, dimensions, hic_path, PUB_ID, resolution, norm, toolsource, featuretype, epigenomicFactors, meta = en

                xDirectional = ""
                yDirectional = ""

    #             b = bytearray()
    #             b.extend(map(ord, s))

                xPivot = binary_search_motif(c1, x1)
                yPivot = binary_search_motif(c2, y1)


                while xPivot<len(kwargs["binary_search_CTCF"][c1]) and int(kwargs["binary_search_CTCF"][c1][xPivot][2])<x2 :
                    while int(kwargs["binary_search_CTCF"][c1][xPivot][3])<x1 and xPivot<len(kwargs["binary_search_CTCF"][c1])-1:
                        xPivot+=1
    #                 b.extend(bytes(ord(dmapped[binary_search_CTCF[c1][xPivot][4]])))
                    if xPivot>=len(kwargs["binary_search_CTCF"][c1]):
                        break
                    xDirectional+=dmapped[kwargs["binary_search_CTCF"][c1][xPivot][4]]
                    xPivot+=1


                while yPivot<len(kwargs["binary_search_CTCF"][c2]) and int(kwargs["binary_search_CTCF"][c2][yPivot][2])<=y2:
                    while int(kwargs["binary_search_CTCF"][c2][yPivot][3])<y1 and yPivot<len(kwargs["binary_search_CTCF"][c2])-1:
                        yPivot+=1
    #                     b.extend(bytes(ord(dmapped[kwargs["binary_search_CTCF"][c2][yPivot][4]])))
                    if yPivot>=len(kwargs["binary_search_CTCF"][c2]):
                        break
                    yDirectional+=dmapped[kwargs["binary_search_CTCF"][c2][yPivot][4]]
                    yPivot+=1


                motifDirection = xDirectional+":"+yDirectional
                
                cursor_t.execute("INSERT INTO imag(key_id, name, dataset, condition, coordinates, numpyarr, viewing_vmax, true_max, hist_rel, hist_true, dimensions, hic_path, PUB_ID, resolution, norm, toolsource, featuretype, epigenomicFactors, motifDirection, meta) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", [key_id, name, dataset, condition, coordinates, numpyarr, viewing_vmax, true_max, hist_rel, hist_true, dimensions, hic_path, PUB_ID, resolution, norm, toolsource, featuretype, epigenomicFactors.tobytes(), motifDirection, meta],)

            print(f"success")
        except Exception as e:
            print(e)

        finally:
            connection_t.commit()
            cursor_t.close()
            cursor_s.close()
            connection_t.close()
            connection_s.close()            


    if not all(i in kwargs for i in ["limit","offset"]):
        raise Exception("need  \"limit\",\"offset\" in kwargs")
    return _readDB(kwargs['offset'],kwargs["limit"])


def mainProg():
    print('BOOTED')
    dbSOURCE = "/nfs/turbo/umms-drjieliu/proj/3C-FeatExt/012625_changeDBcalls/DB_DUMP/TEST_database_14_bin.db"
    dbTARGET = "/nfs/turbo/umms-drjieliu/proj/3C-FeatExt/012625_changeDBcalls/DB_DUMP/database_17_bin.db"


    factor_registry = ["CTCF", 
                       "H3K27ac", 
                       "H3K27me3", 
                       "H3K9me3",
                       "H2AFZ",
                       "H3K36me3",
                       "H3K4me1",
                       "H3K9ac"]


    registerPath = {
        "H1":{x:"" for x in factor_registry},
        "HFF":{x:"" for x in factor_registry},
        "GM12878":{x:"" for x in factor_registry},
        "K562":{x:"" for x in factor_registry}
    }
    print("registering bigWig paths")


    for fkey in factor_registry:
        registerPath["H1"][fkey]=f"/nfs/turbo/umms-drjieliu/proj/4dn/data/Epigenomic_Data/hESC/hESC_{fkey}_hg38.bigWig"
        registerPath["GM12878"][fkey]=f"/nfs/turbo/umms-drjieliu/proj/4dn/data/Epigenomic_Data/GM12878/GM12878_{fkey}_hg38.bigWig"
        registerPath["HFF"][fkey]=f"/nfs/turbo/umms-drjieliu/proj/4dn/data/Epigenomic_Data/human_tissues/foreskin_fibroblast/foreskin_fibroblast_{fkey}_hg38.bigWig"
        registerPath["K562"][fkey]=f"/nfs/turbo/umms-drjieliu/proj/4dn/data/Epigenomic_Data/human_tissues/K562/K562_{fkey}_hg38.bigWig"
        print(f".. init {fkey}")

    registerPath["K562"]["H3K36me3"] = "/nfs/turbo/umms-drjieliu/proj/3C-FeatExt/021725_with_epi/ENCFF633OZC.bigWig"
    registerPath["K562"]["H3K9ac"] = "/nfs/turbo/umms-drjieliu/proj/3C-FeatExt/021725_with_epi/ENCFF178QDA.bigWig"
    registerPath["K562"]["H2AFZ"] = "/nfs/turbo/umms-drjieliu/proj/3C-FeatExt/021725_with_epi/ENCFF202EVH.bigWig"
    registerPath["K562"]["H3K9me3"] = "/nfs/turbo/umms-drjieliu/proj/3C-FeatExt/021725_with_epi/ENCFF559MMQ.bigWig"
    registerPath["GM12878"]["H2AFZ"] = "/nfs/turbo/umms-drjieliu/proj/3C-FeatExt/021725_with_epi/ENCFF992GSC.bigWig"
    registerPath["HFF"]["H2AFZ"]=""
    registerPath["HFF"]["H3K9ac"]=""
    registerPath["HFF"]["H3K9me3"]="/nfs/turbo/umms-drjieliu/proj/3C-FeatExt/021725_with_epi/ENCFF542CZT.bigWig"
    registerPath["HFF"]["H3K36me3"]="/nfs/turbo/umms-drjieliu/proj/3C-FeatExt/021725_with_epi/ENCFF431XVV.bigWig"

    for key in registerPath.keys():
        print(f".Checking {key}")
        for fkey in factor_registry:
            print(f"... {fkey}")
            if registerPath[key][fkey]=="":
                continue
            assert os.path.isfile(registerPath[key][fkey]), "FileNotFoundError: path: " + registerPath[key][fkey] + "does not point to valid or readable bigWig file"
    print("checks passed",end="\n\n")

    inverseTableCell= {
        "HFFc6 (Tier 1)" : "HFF",
        "H1-hESC (Tier 1)": "H1",
        "GM12878": "GM12878",
        "K562": "K562"
    }

    binary_search_CTCF = {f"chr{x}":[] for x in [y for y in range(1,23)]+["X","Y"]}
    with open("/nfs/turbo/umms-drjieliu/proj/3C-FeatExt/021725_with_epi/matches.txt") as y:
        for line in y:
            name,chrom,_x1,_x2,direction=line.split()
            x1=int(_x1)
            x2=int(_x2)
            if name.split("_")[0]=="CTCF":
                binary_search_CTCF[chrom]+=[(name,chrom,x1,x2,direction)]

    for key in [y for y in range(1,23)]+["X"]:
        binary_search_CTCF[f"chr{key}"].sort(key=lambda x: x[2])
    print("binary search tree ready")



    try:
        _createTable(dbTARGET, 10)
    except sqlite3.OperationalError:
        untouch(dbTARGET,10)
        _createTable(dbTARGET, 10)

    hardLimiter = 99129;
    #check length of table

    insert_kwargs = {
        "limit": 8192,
        "offset": 0,
        "inverseTableCell": inverseTableCell,
        "binary_search_CTCF": binary_search_CTCF,
        "registerPath": registerPath
        }

    while insert_kwargs["offset"] < hardLimiter:
        tnow = datetime.datetime.now()
        _readSOURCE_writeVECTOR(dbSOURCE, dbTARGET, 10, **insert_kwargs)
        insert_kwargs["offset"] = insert_kwargs.get("offset", 0) + insert_kwargs.get("limit", 10)
        print(datetime.datetime.now() - tnow, datetime.datetime.now())
        print(insert_kwargs["offset"], insert_kwargs["limit"])

if __name__ == "__main__":
    now = datetime.datetime.now()
    mainProg();
    print(datetime.datetime.now() - now)



