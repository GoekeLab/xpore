.. _quickstart:

Quickstart - Detection of differential RNA modifications
=========================================================

.. note::
   **Updates in xPore v2.2:** xPore is now compatible with genome alignments and RNA004 data — see the table below and the steps for more information.

This page walks through the full pipeline, from raw nanopore signal (POD5, or FAST5 for older runs) through basecalling, alignment, resquiggling, ``xpore dataprep`` and ``xpore diffmod``, to a ranked table of differentially modified sites. If you would rather try xPore on a ready-made dataset first, see :ref:`Running demo data <demo>`.

xPore is compatible with both transcriptome- and genome-aligned data. See below for the minimal commands to run xPore on each:

.. list-table::
   :header-rows: 1
   :widths: 22 48 30

   * -
     - **Dataprep**
     - **Diffmod**
   * - Transcriptome alignment (output transcriptome coordinates)
     - | ``xpore dataprep``
       | ``--eventalign <eventalign.txt>``
       | ``--out_dir <out_dir>``
     - | ``xpore diffmod``
       | ``--config <config.yml>``
   * - Transcriptome alignment (output genome coordinates)
     - | ``xpore dataprep``
       | ``--eventalign <eventalign.txt>``
       | ``--out_dir <out_dir>``
       | ``--genome``
       | ``--transcript_fasta <transcript.fa>``
       | ``--gtf_or_gff <annotation.gtf>``
     - | ``xpore diffmod``
       | ``--config <config.yml>``
   * - Genome alignment
     - | ``xpore dataprep``
       | ``--eventalign <eventalign.txt>``
       | ``--out_dir <out_dir>``
       | ``--kmer_source model_kmer``
     - | ``xpore diffmod``
       | ``--config <config.yml>``

1. Basecalling
--------------

Typically MinKNOW runs `Dorado <https://github.com/nanoporetech/dorado>`_ and basecalls live during acquisition.

If you need to (re)basecall offline — e.g. a different model, a newer Dorado version, or basecalling wasn't enabled during the run — run Dorado directly on the raw signal files (POD5, or FAST5 for older runs). Dorado replaces the older Guppy/Albacore basecallers and supports both RNA002 and RNA004 chemistries — select the model that matches your chemistry. You can find more detail about basecalling at `Oxford Nanopore Technologies <https://nanoporetech.com>`_::

    dorado basecaller <MODEL> </PATH/TO/POD5_DIR> --emit-fastq > <PATH/TO/FASTQ>

For RNA004 data use an ``rna004`` model (e.g. ``rna004_130bps_sup@v5.1.0``); for RNA002 data use an ``rna002`` model. See the `Dorado documentation <https://github.com/nanoporetech/dorado>`_ for the available models and options.

2. Alignment
------------

Align the basecalled reads with `minimap2 <https://github.com/lh3/minimap2>`_. xPore supports both **transcriptome** and **genome** alignments — choose one depending on which coordinate system you want in the output.

**Transcriptome alignment** (align to a transcriptome reference). xPore reports transcriptomic coordinates, or genomic coordinates if you also pass ``--genome`` together with ``--gtf_or_gff`` and ``--transcript_fasta`` to ``xpore dataprep``::

    minimap2 -ax map-ont -uf -t 3 --secondary=no <TRANSCRIPTOME.MMI> <PATH/TO/FASTQ.GZ> > <PATH/TO/SAM> 2>> <PATH/TO/SAM_LOG>
    samtools view -Sb <PATH/TO/SAM> | samtools sort -o <PATH/TO/BAM> - &>> <PATH/TO/BAM_LOG>
    samtools index <PATH/TO/BAM> &>> <PATH/TO/BAM_INDEX_LOG>

**Genome alignment** (align directly to a genome reference; use spliced alignment so that reads spanning introns map correctly). Genome alignments contain reverse-strand reads, so run ``xpore dataprep`` with ``--kmer_source model_kmer`` for these (see :ref:`Command line arguments <cmd>`)::

    minimap2 -ax splice -uf -k14 -t 3 --secondary=no <GENOME.MMI> <PATH/TO/FASTQ.GZ> > <PATH/TO/SAM> 2>> <PATH/TO/SAM_LOG>
    samtools view -Sb <PATH/TO/SAM> | samtools sort -o <PATH/TO/BAM> - &>> <PATH/TO/BAM_LOG>
    samtools index <PATH/TO/BAM> &>> <PATH/TO/BAM_INDEX_LOG>

3. Resquiggling
----------------

Resquiggle (align the raw signal to the reference) to produce the eventalign file. We recommend `f5c <https://github.com/hasindu2008/f5c>`_, an optimised, CPU/GPU-accelerated re-implementation of ``nanopolish eventalign`` that produces equivalent output much faster on large datasets::

    # index the raw signal: use -d <FAST5_DIR> for FAST5, or --slow5 <FILE.blow5> for SLOW5/BLOW5
    f5c index -d <PATH/TO/FAST5_DIR> <PATH/TO/FASTQ_FILE>
    f5c eventalign --reads <PATH/TO/FASTQ_FILE> \
    --bam <PATH/TO/BAM_FILE> \
    --genome <PATH/TO/REFERENCE_FASTA> \
    --rna \
    --signal-index \
    --scale-events \
    --threads 32 > <PATH/TO/eventalign.txt>

For **RNA004** data, add ``--kmer-model <PATH/TO/5-mer-model>``: recent versions of f5c auto-select the 9-mer model for RNA004, which xPore cannot use — xPore requires the 5-mer model (see `xPore issue #215 <https://github.com/GoekeLab/xpore/issues/215>`_).

``nanopolish eventalign`` can be used instead with the same arguments. Note that the ``--genome`` argument here refers to the **alignment reference** (the transcriptome or genome FASTA used in step 2).

4. Preprocess with ``xpore dataprep``
--------------------------------------

Preprocess the eventalign file with ``xpore dataprep``, using the command below that matches your alignment mode from step 2. This step will take approximately 5h for 1 million reads.

.. list-table::
   :header-rows: 1
   :widths: 33 34 33

   * - **Transcriptome alignment (output transcriptome coordinates)**
     - **Transcriptome alignment (output genome coordinates)**
     - **Genome alignment (output genome coordinates)**
   * - | ``xpore dataprep``
       | ``--eventalign <eventalign.txt>``
       | ``--out_dir <OUT_DIR>``
     - | ``xpore dataprep``
       | ``--eventalign <eventalign.txt>``
       | ``--out_dir <OUT_DIR>``
       | ``--genome``
       | ``--transcript_fasta <transcript.fa>``
       | ``--gtf_or_gff <annotation.gtf>``
     - | ``xpore dataprep``
       | ``--eventalign <eventalign.txt>``
       | ``--out_dir <OUT_DIR>``
       | ``--kmer_source model_kmer``

The ``--gtf_or_gff`` and ``--transcript_fasta`` arguments map transcriptomic to genomic coordinates; GTF is the recommended option — if GFF is the only file available, note that it works with GENCODE or ENSEMBL FASTA files, but not UCSC FASTA files. We plan to remove the requirement of FASTA files in a future release.

The output files are stored under ``<OUT_DIR>``:

* ``eventalign.index`` : Index file to access ``eventalign.txt``, the output from f5c/nanopolish eventalign
* ``data.json`` : Preprocessed data for ``xpore diffmod``
* ``data.index`` : File index of ``data.json`` for random access per gene
* ``data.readcount`` : Summary of readcounts per gene
* ``data.log`` : Log file

Run this for each of your samples (e.g. each condition/replicate). Run ``xpore dataprep -h`` or visit our :ref:`Command line arguments <cmd>` to explore the full usage description.

5. Configuration file
-----------------------

Prepare a ``.yml`` configuration file. With this YAML file, you can specify the information of your design experiment, the data directories, the output directory, and the method options::

    data:
        <CONDITION_NAME_1>:
            <REP1>: <PATH/TO/OUT_DIR>
        <CONDITION_NAME_2>:
            <REP1>: <PATH/TO/OUT_DIR>

    out: <PATH/TO/OUT_DIR> # output dir

    # Since xPore v2.2 the default unmodified-signal prior is the RNA004 model.
    # For RNA002 data, set prior to the bundled RNA002 model instead:
    # prior: /path/to/xpore/diffmod/RNA002_5mer_model.csv

See the :ref:`Configuration file page <configuration>` for more details.

6. Run ``xpore diffmod``
--------------------------

Now that we have the data and the configuration file ready, model differential modifications using ``xpore diffmod``::

    xpore diffmod --config <PATH/TO/config.yml>

The output files are generated within the output directory specified in the config:

* ``diffmod.table`` : Result table of differential RNA modification across all tested positions
* ``diffmod.log`` : Log file

Run ``xpore diffmod -h`` or visit our :ref:`Command line arguments <cmd>` to explore the full usage description.

We can rank the significantly differentially modified sites based on the p-value column, e.g. ``pval_<CONDITION_1>_vs_<CONDITION_2>``. The results look like::

    id                position   kmer  diff_mod_rate_KO_vs_WT  pval_KO_vs_WT  z_score_KO_vs_WT  ...  sigma2_unmod  sigma2_mod  conf_mu_unmod  conf_mu_mod  mod_assignment        t-test
    ENSG00000114125  141745412  GGACT               -0.823318  4.241373e-115        -22.803411  ...      5.925238   18.048687       0.968689     0.195429           lower  1.768910e-19
    ENSG00000159111   47824212  GGACT               -0.828023   1.103790e-88        -19.965293  ...      2.686549   13.820089       0.644436     0.464059           lower  5.803242e-18
    ENSG00000159111   47824138  GGGAC               -0.757891   1.898161e-73        -18.128515  ...      3.965195    9.877299       0.861480     0.359984           lower  9.708552e-08
    ENSG00000159111   47824137  GGACA               -0.604056   7.614675e-24        -10.068479  ...      7.164075    4.257725       0.553929     0.353160           lower  2.294337e-10
    ENSG00000114125  141745249  GGACT               -0.514980   2.779122e-19         -8.977134  ...      5.215243   20.598471       0.954968     0.347174           lower  1.304111e-06

7. (Optional) Postprocessing
------------------------------

We can consider only one modification type per k-mer by finding the majority ``mod_assignment`` of each k-mer.
For example, the majority of the modification means of ``GGACT`` (``mu_mod``) is lower than the non-modification counterpart (``mu_unmod``).
We can filter out those positions whose ``mod_assigment`` values are not in line with those of the majority in order to restrict ourselves with one modification type per kmer in the analysis.
This can be done by running ``xpore postprocessing``::

    xpore postprocessing --diffmod_dir <PATH/TO/OUT_DIR>

With this command, we will get the final file in which only kmers with their ``mod_assignment`` different from the majority assigment of the corresponding kmer are removed. The output file ``majority_direction_kmer_diffmod.table`` is generated in the output directory. You can find more details in our paper.

Run ``xpore postprocessing -h`` or visit our :ref:`Command line arguments <cmd>` to explore the full usage description.
