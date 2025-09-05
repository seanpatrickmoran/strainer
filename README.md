# strainer
Takes in 3C data (hi-c, cooler), and tabular 3C architecural annotations, writes their information and **their images** to sqlite3 database.

```
python strainer.py mustache-loops-2B-readpairs.json output-table.db
```

### Prepare JSON or YAML for strainer to find and process features from Hi-C.
Preparation of data for strainer.py requires an input datafile to format the properties of the queries features.

*2Bit genome files are required, simply list all needed.
*"datasets" a list of json objects that require
  1. name
  2. hic_path: path to hic or cooler file for the 3C data.
  3. feature_path: path to the feature caller outputs, i.e. mustache outputs.
  4. genome: such as hg19, needs to be used such that it is the prefix for filepaths at top of json i.e. "HG19_PATH"
  5. resolutions: \[resolution\]. Hi-C resolution, should be availible in Hi-C file and the resolution used for the Feature Caller.

```json
{
  "HG38_PATH": "/2bitPath/hg38.analysisSet.2bit"
  "HG19_PATH": "/2bitPath/hg19.2bit"
  "datasets": [
    {
      "name": "GM12878_4DNFIFLJLIS5",
      "hic_path": "/data/path/4DNFIFLJLIS5.hic",
      "feature_path": "/mustache/resultsPath/4DNFIFLJLIS5_2kb.tsv",
      "dataset": "HFFc6 (Tier 1)",
      "condition": "in situ Hi-C",
      "genome": "hg38",
      "resolutions": [2000],
      "options": {
        "norm": "NONE",
        "toolsource": "mustache",
        "featuretype": "loop"
      }
    },
    {
      "name": "GM12878_4DNFIFLJLIS5",
      "hic_path": "/data/path/GM12878_4DNFIFLJLIS5.hic",
      "feature_path": "/mustache/resultsPath/4DNFIFLJLIS5_5kb.tsv",
      "dataset": "HFFc6 (Tier 1)",
      "condition": "in situ Hi-C",
      "genome": "hg38",
      "resolutions": [5000],
      "options": {
        "norm": "NONE",
        "toolsource": "mustache",
        "featuretype": "loop"
      }
    }
  ]
}
```


Alternatively, YAML can be used instead. Simply follow the template below.


```yaml
HG38_PATH: "/2bitPath/hg38.analysisSet.2bit"
HG19_PATH: "/2bitPath/hg19.2bit"

datasets:
  - name: GM12878_4DNFIFLJLIS5
    hic_path: /data/path/4DNFIFLJLIS5.hic
    feature_path: /mustache/resultsPath/4DNFIFLJLIS5_2kb.tsv
    dataset: HFFc6 (Tier 1)
    condition: in situ Hi-C
    genome: hg38
    resolutions: 
      - 2000
    options:
      norm: NONE
      toolsource: mustache
      featuretype: loop

  - name: GM12878_4DNFIFLJLIS5
    hic_path: /data/path/GM12878_4DNFIFLJLIS5.hic
    feature_path: /mustache/resultsPath/4DNFIFLJLIS5_5kb.tsv
    dataset: HFFc6 (Tier 1)
    condition: in situ Hi-C
    genome: hg38
    resolutions: 
      - 5000
    options:
      norm: NONE
      toolsource: mustache
      featuretype: loop
```


### Munge the written sql table into JSON for inference
Once the SQL table is formed, its necessary to process the SQL table into a JSON parsable for the inference model. Use sql2json.py 

```
# munge everything.
python data2json.py database.db \
  --output data.json \
  --image-dir ./images


# select 5000 only from db1, all db2
python data2json.py db1.db:5000 db2.db \
  --output data.json \
  --image-dir ./images \


# choose only resolutions 2000, 5000
python data2json.py database.db \
  --output data.json \
  --image-dir ./images \
  --resolution 2000 5000
```

For training the model, we can choose to only select labelled data.
```
python data2json.py db1.db \
  --output data.json \
  --image-dir ./images \
  --train

```




