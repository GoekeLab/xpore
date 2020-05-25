.. _quickstart:

Quickstart - Detection of differential RNA modifications
=========================================================

Download the demo dataset from our S3 bucket.::
    wget http://s3.ap-southeast-1.amazonaws.com/all-public-data.store.genome.sg/xpore/diffmod_test_data.tar.gz

After extraction, you will find::
    |-- diffmod.yaml
    |-- data
        |-- Hek293T-WT-rep1
        |-- Hek293T-KO-rep1
    |-- db
        |-- model_kmer.csv
        |-- Homo_sapiens.GRCh38.cdna.mmi

Each dataset under the ``data`` directory contains the following directories:
``fast5`` : Raw signal FAST5 files.
``fastq`` : Basecalled reads.
``bamtx`` : Transcriptome-aligned sequence.
``nanopolish``: Eventalign files using `nanopolish <https://nanopolish.readthedocs.io/en/latest>`


