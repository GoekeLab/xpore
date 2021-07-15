.. _quickstart:

Quickstart - Detection of differential RNA modifications
=========================================================

Download and extract the demo dataset from our `zenodo <https://zenodo.org/record/5103099/files/demo.tar.gz>`_::

    wget http://s3.ap-southeast-1.amazonaws.com/all-public-data.store.genome.sg/xpore/demo.tar.gz
    tar -xvf demo.tar.gz

After extraction, you will find::
    
    demo
    |-- Hek293T_config.yml  # configuration file
    |-- data
        |-- HEK293T-METTL3-KO-rep1  # dataset dir
        |-- HEK293T-WT-rep1 # dataset dir
    |-- demo.gtf # GTF (general transfer format) file required for mapping transcriptomic to genomic coordinates  
    |-- demo.fa # transcriptome reference FASTA file required for mapping transcriptomic to genomic coordinates

Each dataset under the ``data`` directory contains the following directories:

* ``fast5`` : Raw signal, FAST5 files
* ``fastq`` : Basecalled reads, FASTQ file
* ``bamtx`` : Transcriptome-aligned sequence, BAM file
* ``nanopolish``: Eventalign files obtained from `nanopolish eventalign <https://nanopolish.readthedocs.io/en/latest/quickstart_eventalign.html>`_

Note that the FAST5, FASTQ and BAM files are required to obtain the eventalign file with Nanopolish, xPore only requires the eventalign file. See our :ref:`Data preparation page <preparation>` for details to obtain the eventalign file from raw reads.

1. Preprocess the data for each data set using ``xpore dataprep``. Note that the ``--gtf_path_or_url`` and ``--transcript_fasta_paths_or_urls`` arguments are required to map transcriptomic to genomic coordinates when the ``--genome`` option is chosen, so that xPore can run based on genome coordinates. (This step will take approximately 5h for 1 million reads)::

    # Within each dataset directory i.e. demo/data/HEK293T-METTL3-KO-rep1 and demo/data/HEK293T-WT-rep1, run
    xpore-dataprep \
    --eventalign nanopolish/eventalign.txt \
    --gtf_path_or_url demo.gtf \
    --transcript_fasta_paths_or_urls demo.fa \
    --out_dir dataprep \
    --genome  

The output files are stored under ``dataprep`` in each  dataset directory:

* ``eventalign.index`` : Index file to access ``eventalign.txt``, the output from nanopolish eventalign
* ``data.json`` : Preprocessed data for ``xpore-diffmod``
* ``data.index`` : File index of ``data.json`` for random access per gene
* ``data.readcount`` : Summary of readcounts per gene
* ``data.log`` : Log file

Run ``xpore dataprep -h`` or visit our :ref:`Command line arguments <cmd>` to explore the full usage description. 

2. Prepare a ``.yml`` configuration file. With this YAML file, you can specify the information of your design experiment, the data directories, the output directory, and the method options.
In the demo directory, there is an example configuration file ``Hek293T_config.yaml`` available that you can use as a starting template.
Below is how it looks like::

    notes: Pairwise comparison without replicates with default parameter setting.

    data:
        KO:
            rep1: ./data/HEK293T-METTL3-KO-rep1/dataprep 
        WT:
            rep1: ./data/HEK293T-WT-rep1/dataprep

    out: ./out # output dir


See the :ref:`Configuration file page <configuration>` for more details.

3. Now that we have the data and the configuration file ready for modelling differential modifications using ``xpore-diffmod``. 

::

    # At the demo directory where the configuration file is, run.
    xpore diffmod --config Hek293T_config.yml

The output files are generated within the ``out`` directory:

* ``diffmod.table`` : Result table of differential RNA modification across all tested positions
* ``diffmod.log`` : Log file

Run ``xpore diffmod -h`` or visit our :ref:`Command line arguments <cmd>` to explore the full usage description.

We can rank the significantly differentially modified sites based on ``pval_HEK293T-KO_vs_HEK293T-WT``. The results are shown below.::

    id                position   kmer  diff_mod_rate_KO_vs_WT  pval_KO_vs_WT  z_score_KO_vs_WT  ...  sigma2_unmod  sigma2_mod  conf_mu_unmod  conf_mu_mod  mod_assignment        t-test
    ENSG00000114125  141745412  GGACT               -0.823318  4.241373e-115        -22.803411  ...      5.925238   18.048687       0.968689     0.195429           lower  1.768910e-19
    ENSG00000159111   47824212  GGACT               -0.828023   1.103790e-88        -19.965293  ...      2.686549   13.820089       0.644436     0.464059           lower  5.803242e-18
    ENSG00000159111   47824138  GGGAC               -0.757891   1.898161e-73        -18.128515  ...      3.965195    9.877299       0.861480     0.359984           lower  9.708552e-08
    ENSG00000159111   47824137  GGACA               -0.604056   7.614675e-24        -10.068479  ...      7.164075    4.257725       0.553929     0.353160           lower  2.294337e-10
    ENSG00000114125  141745249  GGACT               -0.514980   2.779122e-19         -8.977134  ...      5.215243   20.598471       0.954968     0.347174           lower  1.304111e-06

4. (Optional) We can consider only one modification type per k-mer by finding the majority ``mod_assignment`` of each k-mer. 
For example, the majority of the modification means of ``GGACT`` (``mu_mod``) is lower than the non-modification counterpart (``mu_unmod``). 
We can filter out those positions whose ``mod_assigment`` values are not in line with those of the majority in order to restrict ourselves with one modification type per kmer in the analysis.
This can be done by running ``xpore postprocessing``.

::

    xpore postprocessing --diffmod_dir out

With this command, we will get the final file in which only kmers with their ``mod_assignment`` different from the majority assigment of the corresponding kmer are removed. The output file ``majority_direction_kmer_diffmod.table`` is generated in the ``out`` directtory. You can find more details in our paper.

Run ``xpore postprocessing -h`` or visit our :ref:`Command line arguments <cmd>` to explore the full usage description.
