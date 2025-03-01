import sqlite3
from CCCImageClasses import C3Image4, extractor
#    

def call(PATH,TIMEOUT):

    connection = sqlite3.connect(PATH, timeout=TIMEOUT)  # Set timeout to 10 seconds
    cursor = connection.cursor()
    return connection,cursor


def _createTable(dbPATH, timeout,**kwargs):
    connection,cursor=call(dbPATH,timeout)
    try:

#         cursor = connection.execute("CREATE TABLE imag(rowid, name, coordinates, numpyarr, viewing_vmax, dimensions, hic_path, PUB_ID, resolution, meta)")
        cursor = connection.execute("CREATE TABLE imag(key_id, name, dataset, condition, coordinates, numpyarr, viewing_vmax, true_max, hist_rel, hist_true, dimensions, hic_path, PUB_ID, resolution, norm, toolsource, featuretype, meta)")
        print("table make; success")
    finally:
        cursor.close()
        connection.close()

def _writeManyToTable(dbPATH,timeout,**kwargs):
    def _write_db(data):
        connection,cursor=call(dbPATH,timeout)
        try:
            data = tuple(x for x in data)
            cursor.executemany("INSERT INTO imag VALUES(:key_id, :name, :dataset, :condition, :coordinates, :numpyarr, :viewing_vmax, :true_max, :hist_rel, :hist_true, :dimensions, :hic_path, :PUB_ID, :resolution, :norm, :toolsource, :featuretype, :meta)", data)
            print(f"success")

        except Exception as e:
            template = "An exception of type {0} occurred. Arguments:\n{1!r}"
            message = template.format(type(e).__name__, e.args)
            print(message)
            print(e)

        finally:
            connection.commit()
            cursor.close()
            connection.close()
    
    for val in ["dataset", "condition","norm"]:
        kwargs[val] = kwargs.get(val, "")
    if not all(i in kwargs for i in ["image_np_dict","prefix_name","hic_path","PUB_ID","resolution","norm"]):
        raise Exception("\"image_np_dict\",\"prefix_name\",\"hic_path\",\"PUB_ID\", \"resolution\", or \"norm\" not in kwargs")
        return
    
    data = []
    key_id        = kwargs["key_id"]
    image_np_dict = kwargs["image_np_dict"]
    prefix_name   = kwargs["prefix_name"]
    resolution    = kwargs["resolution"]
    hic_path      = kwargs["hic_path"]
    PUB_ID        = kwargs["PUB_ID"]
    dataset       = kwargs["dataset"]
    condition     = kwargs["condition"]
    norm          = kwargs["norm"] 
    toolsource    = kwargs["toolsource"]
    featuretype   = kwargs["featuretype"]
    for i in image_np_dict.keys():
        if len(data)>999:
            _write_db(data)
            data = []
        focus = image_np_dict[i]
        fields = extractor(focus, 'original coordinates', 'numpy_window', 'viewing_vmax', 'true_max', 'hist_rel', 'hist_true')
        # print(C3Image4(key_id, f"{prefix_name}_#{i}", *[fields[k] for k in fields.keys()],focus['numpy_window'].shape[0],resolution,hic_path,PUB_ID,genome="hg38", dataset=dataset,condition=condition,norm=norm).entity)
        data += [C3Image4(key_id, f"{prefix_name}_#{i}", *[fields[k] for k in fields.keys()],focus['numpy_window'].shape[0],resolution,hic_path,PUB_ID,genome="hg38", dataset=dataset,condition=condition,norm=norm, toolsource=toolsource, featuretype=featuretype).entity]
        key_id += 1
    _write_db(data)
    return "", key_id

def _writeSingularToTable(dbPATH,timeout,**kwargs):
    data = (
            kwargs["C3Image4"].entity
            )
    def _write_db(data):
        connection,cursor=call(dbPATH,timeout)
        try:
            data = tuple(x for x in data)
    #         print(data)
            cursor.executemany("INSERT INTO imag VALUES(:key_id, :name, :dataset, :condition, :coordinates, :numpyarr, :viewing_vmax, :true_max, :hist_rel, :hist_true, :dimensions, :hic_path, :PUB_ID, :resolution, :norm, :toolsource, :featuretype, :meta)", data)
            print(f"success")
        except Exception as e:
            print(e)

        finally:
            connection.commit()
            cursor.close()
            connection.close()
    if "C3Image4" not in kwargs:
        raise Exception("need \"C3Image4\": C3Image4 object in kwargs")
    _write_db(kwargs["C3Image4"])


def _readFirstEntryTable(dbPATH,timeout,**kwargs):
    def _readDB(name):
        connection,cursor=call(dbPATH,timeout)
        try:
            cursor.row_factory = sqlite3.Row
            params = (name,)
            cursor.execute("SELECT * FROM imag WHERE name = ? ", params)
            print(f"success")
            row = cursor.fetchone()
        except Exception as e:
            print(e)
        finally:
            cursor.close()
            connection.close()
        return row
    if "name" not in kwargs:
        raise Exception("need query \'name\' in kwargs")
    return _readDB(kwargs['name'])
    

def _readMatchAllTable(dbPATH,timeout,**kwargs):
    def _readDB(name):
        connection,cursor=call(dbPATH,timeout)
        try:
            cursor.row_factory = sqlite3.Row
            params = (name,)
            cursor.execute("SELECT * FROM imag WHERE name = ? ", params)
            print(f"success")
            reply = []
            for en in cursor.fetchall():
                reply+=[en]
        except Exception as e:
            template = "An exception of type {0} occurred. Arguments:\n{1!r}"
            message = template.format(type(e).__name__, e.args)
            print(message)
            print(e)

        finally:
            cursor.close()
            connection.close()
        return reply
    if "name" not in kwargs:
        raise Exception("need query \'name\' in kwargs")
    return _readDB(kwargs['name'])


def _readMatchHiCPathTable(dbPATH,timeout,**kwargs):
    def _readDB(name):
        connection,cursor=call(dbPATH,timeout)
        try:
            cursor.row_factory = sqlite3.Row
            params = (name,)
            cursor.execute("SELECT * FROM imag WHERE hic_path = ? ", params)
            print(f"success")
            reply = []
            for en in cursor.fetchall():
                reply+=[en]
        except Exception as e:
            print(e)
        finally:
            cursor.close()
            connection.close()
        return reply
    if "hic_path" not in kwargs:
        raise Exception("need query \'hic_path\' in hic_path")
    return _readDB(kwargs['hic_path'])

def _check_head(dbPATH,timeout,**kwargs):
    connection,cursor=call(dbPATH,timeout)
    promise = False
    try:
        cursor.execute("SELECT * FROM imag ORDER BY ROWID ASC LIMIT 10")#, (kwargs['tablename'],))
        print(cursor.fetchall())
        print(f"success")
        promise = True
    except Exception as e:
        print("connection failure")
        print(e)
    finally:
        cursor.close()
        connection.close()
        print(promise)
        return promise

        
def _check_tail(dbPATH,timeout,**kwargs):
    connection,cursor=call(dbPATH,timeout)
    promise = False
    try:
        cursor.execute("SELECT * FROM imag ORDER BY ROWID DESC LIMIT 10")#, (kwargs['tablename'],))
        print(cursor.fetchall())
        print(f"success")
        promise = True
    except Exception as e:
        print("connection failure")
        print(e)
    finally:
        cursor.close()
        connection.close()
        return promise

def interact_table(databasePath,timeout,f_name,**kwargs):
    fdict = {"_createTable":_createTable,
             "_writeManyToTable":_writeManyToTable,
             "_readFirstEntryTable":_readFirstEntryTable,
             "_readMatchHiCPathTable":_readMatchHiCPathTable,
             "_check_head":_check_head,
             "_check_tail":_check_tail,
             "_writeSingularToTable":_writeSingularToTable}
    return fdict[f_name](databasePath,timeout,**kwargs)
    

def untouch(dbPATH,timeout):
    connection,cursor=call(dbPATH,timeout)
    try:
        cursor.execute("DROP TABLE imag")#, (kwargs['tablename'],))
#         print(cursor.fetchall())
        print(f"success")
    except Exception as e:
        print("connection failure")
        print(e)
    finally:
        connection.commit()
        cursor.close()
        connection.close()


#databasePATH = '/nfs/turbo/umms-drjieliu/proj/GrandCanyons/image_db/images2.db'
#insert_kwargs = {"hic_path":'/nfs/turbo/umms-drjieliu/proj/4dn/data/bulkHiC/haley_HSPC_hic/6527.hic'}
#rows = interact_table(databasePATH,10,"_readMatchHiCPathTable",**insert_kwargs)
#
#interact_table(databasePATH,10,"_check_tail")
#


