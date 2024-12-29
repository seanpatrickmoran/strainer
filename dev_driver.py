import json
import os
import sqlite3
    
def call(PATH,TIMEOUT):
    connection = sqlite3.connect(PATH, timeout=TIMEOUT)  # Set timeout to 10 seconds
    cursor = connection.cursor()
    return connection,cursor

def checkSQL(databasePATH,timeout):
    assert os.path.isfile(databasePATH), "FileNotFoundError: databasePATH does not point to valid sqlite3 database"
    connection,cursor=call(databasePATH,timeout)
    promise = False
    try:
        cursor = connection.execute("SELECT name FROM sqlite_master WHERE type='table'")
        promise = True
        print(cursor.fetchall())
    except Exception as e:
        print("connection failure")
        print(e)
    finally:
        cursor.close()
        connection.close()
        return promise


from hicImageViews import collect_numpy_matrices
from touchSql import interact_table, _createTable, _writeManyToTable, _writeSingularToTable, _readFirstEntryTable, _readMatchAllTable, _readMatchHiCPathTable, _check_head, _check_tail, untouch
# from CCCImageClasses import C3Image4, extractor
def writeFunctionCalls(jsonListOfDicts,databasePATH,**kwargs):
    assert checkSQL(databasePATH,10), "No table exists at database, run 'interact_table(dbp,10,\"_createTable\")''"
    assert interact_table(dbp,10,"_check_head"), "No table exists at database, run 'interact_table(dbp,10,\"_createTable\")''"
    for jDict in jsonListOfDicts:
        hic_path    = jDict["hic_path"]
        featurePath = jDict["featurePath"]
        resolution  = jDict["resolution"]
        sqliteARGS  = jDict["sqliteARGS"]
        norm        = jDict["norm"]
        PUB_ID      = jDict["PUB_ID"]
        dataset     = jDict["dataset"]
        condition   = jDict["condition"]
        dimensions  = jDict["dimensions"]
        metadata    = jDict["meta-data"]
        nickname    = jDict["nickname"]
        
        print(f"writing {nickname} entries to {databasePATH}",end="")
        hicViewDict = collect_numpy_matrices(hic_path, featurePath, norm, int(resolution), int(dimensions))
        
        populousIndices = only_populated_windows(numpy_hic_dictionary,dimensions,choose_scaler=kwargs.get("choose_scaler",2))
        populousImgDict = {x:numpy_hic_dictionary[x] for x in populousIndices}
        
#        control_feat = only_populated_windows(hic_control_numpy,32,4)
#        populousDicts
        print(f".",end="")

#         hicViewDict[lineNumber] where lineNumber are lines from featurePath tabular file
#
#         hicViewDict[lineNumber] = {
#                     "original coordinates": [c1,x1,x2,c2,y1,y2], #c1 and c2 Str, else Int
#                     "window coordinates"  : [r1,r2,r3,r4],
#                     "numpy_window"        : 2D np.float32 Array, 
#                     "mumpy_window"        : 2D np.float32 Array, 
#                     "viewing_vmax"        : Int,
#
#
        insert_kwargs = {
            "image_np_dict": hicViewDict,
            "prefix_name"  : nickname,
            "resolution"   : resolution,
            "hic_path"     : hic_path,
            "PUB_ID"       : PUB_ID,
            "dataset"      : dataset,
            "condition"    : condition,
            "metadata"     : metadata}
    
        print(f".",end=" ")
        interact_table(databasePATH,int(sqliteARGS["SIG-timeout"]),sqliteARGS["touchSqlCmd"],**insert_kwargs)
        print(f"__@__@__@__@__@__")


if __name__ == "__main__":
    jsonPATH = "test.json"
    with open(jsonPATH) as f:
        d = json.load(f)

    
    dbp = "/home/spmoran/temp_smoran/Canyons/122724_test_mass_sqlwrites/databse1.db"
    try:
        if not interact_table(dbp,10,"_check_head"):
            interact_table(dbp,10,"_createTable")
    except:
        raise ("could not form table")
    
    finally:
        untouch(dbp,10)
    
    
    dbp = "/home/spmoran/temp_smoran/Canyons/122724_test_mass_sqlwrites/databse1.db"
    try:
        if not interact_table(dbp,10,"_check_head"):
            interact_table(dbp,10,"_createTable")
            writeFunctionCalls(d,dbp)
    except:
        raise ("could not form table")
    finally:
        print("all done!")
        interact_table(dbp,10,"_check_head")
