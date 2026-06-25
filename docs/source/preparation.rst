.. _preparation:

Data preparation from raw reads
===================================

1. After obtaining the raw signal files (POD5, or FAST5 for older runs), the first step is to basecall them. Below is an example using `Dorado <https://github.com/nanoporetech/dorado>`_, Oxford Nanopore's current basecaller (it replaces the older Guppy/Albacore basecallers and supports both RNA002 and RNA004 chemistries — select the model that matches your chemistry). You can find more detail about basecalling at `Oxford Nanopore Technologies <https://nanoporetech.com>`_::

    dorado basecaller <MODEL> </PATH/TO/POD5_DIR> --emit-fastq > <PATH/TO/FASTQ>

   For RNA004 data use an ``rna004`` model (e.g. ``rna004_130bps_sup@v5.1.0``); for RNA002 data use an ``rna002`` model. See the `Dorado documentation <https://github.com/nanoporetech/dorado>`_ for the available models and options.

2. Align the basecalled reads with `minimap2 <https://github.com/lh3/minimap2>`_. xPore supports both **transcriptome** and **genome** alignments — choose one depending on which coordinate system you want in the output.

   **Transcriptome alignment** (align to a transcriptome reference). xPore reports transcriptomic coordinates, or genomic coordinates if you also pass ``--genome`` together with ``--gtf_or_gff`` and ``--transcript_fasta`` to ``xpore dataprep``::

       minimap2 -ax map-ont -uf -t 3 --secondary=no <TRANSCRIPTOME.MMI> <PATH/TO/FASTQ.GZ> > <PATH/TO/SAM> 2>> <PATH/TO/SAM_LOG>
       samtools view -Sb <PATH/TO/SAM> | samtools sort -o <PATH/TO/BAM> - &>> <PATH/TO/BAM_LOG>
       samtools index <PATH/TO/BAM> &>> <PATH/TO/BAM_INDEX_LOG>

   **Genome alignment** (align directly to a genome reference; use spliced alignment so that reads spanning introns map correctly). Genome alignments contain reverse-strand reads, so run ``xpore dataprep`` with ``--kmer_source model_kmer`` for these (see :ref:`Command line arguments <cmd>`)::

       minimap2 -ax splice -uf -k14 -t 3 --secondary=no <GENOME.MMI> <PATH/TO/FASTQ.GZ> > <PATH/TO/SAM> 2>> <PATH/TO/SAM_LOG>
       samtools view -Sb <PATH/TO/SAM> | samtools sort -o <PATH/TO/BAM> - &>> <PATH/TO/BAM_LOG>
       samtools index <PATH/TO/BAM> &>> <PATH/TO/BAM_INDEX_LOG>

3. Resquiggle (align the raw signal to the reference) to produce the eventalign file. We recommend `f5c <https://github.com/hasindu2008/f5c>`_, an optimised, CPU/GPU-accelerated re-implementation of ``nanopolish eventalign`` that produces equivalent output much faster on large datasets::

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

   ``nanopolish eventalign`` can be used instead with the same arguments. Note that the ``--genome`` argument here refers to the **alignment reference** (the transcriptome or genome FASTA used in step 2), not xPore's ``--genome`` flag.
