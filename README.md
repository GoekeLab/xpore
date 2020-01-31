# xpore

Dillinger is a collection of tools for Nanopore data analysis.

### Installation

xpore requires [Python3](https://www.python.org) to run.

To install, run

```sh
$ python setup.py install
```

### Detection of differential modification
Download the [datasets](s3://all-public-data.store.genome.sg/xpore/tests) and [annotations](s3://all-public-data.store.genome.sg/xpore/db) from our s3 bucket.

You will find 
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
    |-- out
        |--
```
Each dataset contains the following directories
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

Now the data are ready for estimating differential modification using `xpore-diffmod`. To , Hek293T_directRNA and Hek293T-METTL3-KO_directRNA, we can use the example confgiuration file 
```sh
$ xpore-diffmod --config_filepath ./tests/config/diffmod.ini --n_processes 4
```
### License
MIT

