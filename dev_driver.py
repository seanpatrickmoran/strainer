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


from hicImageViews import collect_numpy_matrices, only_populated_windows
from touchSql import interact_table, _createTable, _writeManyToTable, _writeSingularToTable, _readFirstEntryTable, _readMatchAllTable, _readMatchHiCPathTable, _check_head, _check_tail, untouch
# from CCCImageClasses import C3Image4, extractor
def writeFunctionCalls(jsonListOfDicts,databasePATH,**kwargs):
    assert checkSQL(databasePATH,10), "No table exists at database, run 'interact_table(dbp,10,\"_createTable\")''"
    assert interact_table(dbp,10,"_check_head"), "No table exists at database, run 'interact_table(dbp,10,\"_createTable\")''"
    index_offset = 0;
    for jDict in jsonListOfDicts:

        print(jDict)
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
        toolsource  = jDict["tool-source"]
        featuretype = jDict["feat-type"]
        

# {'sqliteARGS': {'SIG-timeout': '10', 'touchSqlCmd': '_writeManyToTable'},
#  'hic_path': '4DNFIIMZB6Y9',
#  'resolution': '2000',
#  'norm': 'NONE',
#  'tool-source': 'mustache',
#  'feat-type': 'loop',
#  'dimensions': '64',
#  'meta-data': {'description': 'Hi-C experiments on H1-hESC crosslinked with both FA and DSG, digested with DpnII.',
#   'accession': '4DNESX75DD7R',
#   'Experiment Set Type': 'Replicate',
#   'Organism': 'H. sapiens',
#   'Sample Type': 'stem cells',
#   'Sample': 'H1-hESC (Tier 1)',
#   'Experiment Type(s)': 'in situ Hi-C',
#   'Modification Type': 'None',
#   'Treatment Type': 'None',
#   'Assay Details': 'Enzyme: DpnII',
#   'PMCID': 'PMC8446342'},
#  'author1': 'Akgol Oksuz B',
#  'nickname': '4DNFIIMZB6Y9_2000_mustache',
#  'featurePath': '/nfs/turbo/umms-drjieliu/proj/3C-FeatExt/020925_MassCalls/4DNFIIMZB6Y9_2kb.tsv',
#  'PUB_ID': 'PMC8446342',
#  'dataset': 'H1-hESC (Tier 1)',
#  'condition': 'in situ Hi-C',
#  'drawer-name': '/2021/Akgol/H1-hESC/mustache/4DNFIIMZB6Y9_2000_mustache.tsv'}





        
        print(f"writing {nickname} entries to {databasePATH}",end="")
        hicViewDict = collect_numpy_matrices(hic_path, featurePath, norm, int(resolution), int(dimensions))
        
        try:
            populousIndices = only_populated_windows(hicViewDict,dimensions,choose_scaler=kwargs.get("choose_scaler",2))
            populousImgDict = {x:hicViewDict[x] for x in populousIndices}
        
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
                "metadata"     : metadata,
                "toolsource"   : toolsource,
                "featuretype"  : featuretype,
                "key_id"       : index_offset}
        
            print(f".",end=" ")
            _,resOffset = interact_table(databasePATH,int(sqliteARGS["SIG-timeout"]),sqliteARGS["touchSqlCmd"],**insert_kwargs)
            index_offset = resOffset+1 
            print(f"__@__@__@__@__@__")

        except Exception as e:
            print(e, index_offset)
            continue


if __name__ == "__main__":
    #jsonPATH = "/nfs/turbo/umms-drjieliu/usr/temp_smoran/Canyons/122824_loops_2_json_write/122924_test_mustache_GM12878.json"
    jsonPATH = "/nfs/turbo/umms-drjieliu/proj/3C-FeatExt/020925_MassCalls/021325_for_extended_DB.json"
    with open(jsonPATH) as f:
        d = json.load(f)


    #dbp = "/home/spmoran/temp_smoran/Canyons/122724_test_mass_sqlwrites/databse6_binary.db"
#    dbp = "/nfs/turbo/umms-drjieliu/proj/3C-FeatExt/012625_changeDBcalls/DB_DUMP/database_8_bin.db"
    dbp = "/nfs/turbo/umms-drjieliu/proj/3C-FeatExt/012625_changeDBcalls/DB_DUMP/database_11_bin.db"
    #vecdb = "/nfs/turbo/umms-drjieliu/proj/3C-FeatExt/012625_changeDBcalls/DB_DUMP/sql_vec.db"
 
    try:
        if not interact_table(dbp,10,"_check_head"):
            interact_table(dbp,10,"_createTable")
    except:
        raise ("could not form table")
    
    finally:
        untouch(dbp,10)
    
    
    try:
        if not interact_table(dbp,10,"_check_head"):
            interact_table(dbp,10,"_createTable")
            writeFunctionCalls(d,dbp)
    except:
        raise ("could not form table")
    finally:
        print("all done!")
        interact_table(dbp,10,"_check_head")
