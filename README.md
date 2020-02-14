# xpore

Dillinger is a collection of tools for Nanopore data analysis.

### Installation

xpore requires [Python3](https://www.python.org) to run.

To install our xpore package and its dependencies, run

```sh
$ python setup.py install
$ pip install -r requirements.txt 
$ pyensembl install --release 91 --species homo_sapiens
```

### Detection of differential modification
Downlaod the dataset from our S3 bucket.

```sh
$ wget http://s3.ap-southeast-1.amazonaws.com/all-public-data.store.genome.sg/xpore/diffmod_test_data.tar.gz
```

After extraction, you will find 
```
.
|-- db
    |-- model_kmer.csv
    |-- mapping
        |-- ENSG00000123989.csv
        |-- ENSG00000168496.csv
        |-- ENSG00000170144.csv
        |-- ENSG00000204388.csv
|-- tests
    |-- config
        |-- diffmod.ini 
    |-- data 
        |-- GohGIS_Hek293T_directRNA_Rep2 
        |-- GohGIS_Hek293T-METTL3-KO_directRNA_Rep2_Run1
```

Each dataset under the `tests/data` directory contains the following directories
```
bamtx 
fast5
fastq
nanopolish
dataprep
```

First, we need to preprocess the data for each data set using `xpore-dataprep`. Within each dataset directory, run
```sh
$ xpore-dataprep --eventalign nanopolish/eventalign.txt \
--summary nanopolish/summary.txt \
--bamtx bamtx/aligned.bam \
--mapping ../../../db/mapping \
--out_dir dataprep --n_processes 4
```
Output files: `eventalign.hdf5`, `eventalign.log`, `data.json`, `data.index` , `data.log`

Now the data are ready for estimating differential modification using `xpore-diffmod`. To call the differential modification between Hek293T_directRNA and Hek293T-METTL3-KO_directRNA, we can use the example confgiuration file 
```sh
$ xpore-diffmod --config_filepath ./tests/config/diffmod.ini --n_processes 4 --save_table
```
Output files: `out/diffmod.ini/diffmod.table` `out/diffmod.ini/diffmod.log`
### License
MIT

