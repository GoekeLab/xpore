# xPore

### Installation

xpore requires [Python3](https://www.python.org) to run.

To install our xpore package and its dependencies, run

```sh
$ python setup.py install
$ pyensembl install --release 91 --species homo_sapiens
```

### Detection of differential modification
Download the demo dataset from our S3 bucket.

```sh
$ wget http://s3.ap-southeast-1.amazonaws.com/all-public-data.store.genome.sg/xpore/diffmod_test_data.tar.gz
```

After extraction, you will find 
```
|-- diffmod.ini 
|-- data 
    |-- GohGIS_Hek293T_directRNA_Rep2 
    |-- GohGIS_Hek293T-METTL3-KO_directRNA_Rep2_Run1
|-- db
	|-- model_kmer.csv
	|-- Homo_sapiens.GRCh38.cdna.mmi
```

Each dataset under the `data` directory contains the following directories
```
fast5 
fastq
bamtx
nanopolish
```

(To run the demo, you can skip the first two steps, because `bamtx` and `nanopolish` are already there.)

1. First, align the sequence to transcriptome using `minimap2`. Within each dataset directory, run
```sh
$mkdir bamtx
$minimap2 -ax map-ont -uf -t 8 --secondary=no \
db/Homo_sapiens.GRCh38.cdna.mmi \
fastq/basecalled.fastq.gz > bamtx/aligned.sam 2>> bamtx/aligned.sam.log
$samtools view -Sb bamtx/aligned.sam | samtools sort -o bamtx/aligned.bam - &>> bamtx/aligned.bam.log
$samtools index bamtx/aligned.bam &>> bamtx/aligned.bam.index.log
```

2. Second, align the raw signal to the reference sequence using `nanopolish`. Within each dataset directory, run
```sh
$mkdir nanopolish
$nanopolish index -d fast5 fastq/basecalled.fastq
$nanopolish eventalign --reads fastq/basecalled.fastq \
--bam bamtx/aligned.bam \
--genome $proj_dir/db/Homo_sapiens.GRCh38.cdna.ncrna_wtChrIs_modified.fa \
--scale-events \
--signal-index \
--summary nanopolish/summary.txt \
--threads 4 > nanopolish/eventalign.txt
```

3. Third, we need to preprocess the data for each data set using `xpore-dataprep`. Within each dataset directory, run
```sh
$xpore-dataprep --eventalign nanopolish/eventalign.txt \
        --summary nanopolish/summary.txt \
        --read_count_min 10 \
        --read_count_max 500 \
        --out_dir dataprep \
        --genome \
        --n_processes 2
```

Output files: `eventalign.hdf5`, `eventalign.log`, `data.json`, `data.index` , `data.readcount`, `data.log`

4. Finally, now that the data are ready for estimating differential modification using `xpore-diffmod`. To call the differential modification between Hek293T_directRNA and Hek293T-METTL3-KO_directRNA, we can use the example confgiuration file `diffmod.yaml` by running

```sh
$ xpore-diffmod --config_filepath ./tests/config/diffmod.yaml --n_processes 2 --save_table
```

Output files: `out/diffmod.yaml/diffmod.table` `out/diffmod.ini/diffmod.log`.


We can rank the significanly differentially modified sites based on `p_w_mod_HEK293T-KO_vs_HEK293T-WT`. The results are shown below.

```
idx position    kmer    coverage_HEK293T-METTL3-KO-rep1 coverage_HEK293T-WT-rep1    mu_unmod    mu_mod  sigma2_unmod    sigma2_mod  conf_mu_unmod   conf_mu_mod mod_assignment  w_mod_HEK293T-METTL3-KO-rep1    w_mod_HEK293T-WT-rep1   p_w_mod_HEK293T-KO_vs_HEK293T-WT    w_mod_mean_diff_HEK293T-KO_vs_HEK293T-WT    z_score_HEK293T-KO_vs_HEK293T-WT
ENSG00000114125 141745412   GGACT   167.00000000000009  77.99999999999997   123.64198305264463  117.62845573389104  5.925237677872507   18.048686652338954  0.9686894976263544  0.19542869203353666 lower   0.122081280515318   0.9453989811254184  4.241373321581284e-115  -0.8233177006101003 -22.803411286539568
ENSG00000159111 47824212    GGACT   115.0   56.99999999999999   126.04060818513784  120.32286061375729  2.6865489759165357  13.820088773078876  0.6444364495129247  0.4640590683780786  lower   0.12675220252612124 0.9547753654686716  1.1037896604310229e-88  -0.8280231629425505 -19.965292828395782
ENSG00000159111 47824138    GGGAC   105.0   55.0    118.8431338231134   115.38819698904041  3.965195468986447   9.877299131873366   0.8614802593826912  0.35998415978405274 lower   0.2420911154423771  0.9999818188429512  1.8981606007746968e-73  -0.7578907034005742 -18.128515052229204
ENSG00000159111 47824213    CGGAC   111.0   58.0    120.0023181524575   125.19517965940052  16.09490223670403   2.517386156153043   0.7770385571640749  0.1754346779458279  higher  0.6714153939678753  1.7240784800524122e-05  3.0229603394241693e-51  0.6713981531830748  15.058784020930725
ENSG00000114125 141745411   GGGAC   166.00000000000003  85.00000000000003   117.03987287411272  120.2784827935068   8.177643930183974   2.8216439842252683  0.6933138912876065  0.5304746373270921  higher  0.7056088802507199  0.12806065000998446 4.010247723322406e-30   0.5775482302407354  11.403633554535956
```

### License
MIT

