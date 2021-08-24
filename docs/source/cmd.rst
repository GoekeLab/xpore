.. _cmd:

Command line arguments
=======================

We provide 2 main scripts to run the analysis of differential RNA modifications as the following.

``xpore-dataprep``
********************

* Input

Output files from ``nanopolish eventalgin``. Please refer to :ref:`Data preparation <preparation>` for the full Nanopolish command.

=================================   ==========  ===================  ============================================================================================================
Argument name                       Required    Default value         Description
=================================   ==========  ===================  ============================================================================================================
--eventalign=FILE                   Yes         NA                    Eventalign filepath, the output from nanopolish.         
--out_dir=DIR                       Yes         NA                    Output directory.
--gtf_path_or_url                   No          NA                    GTF file path or url used for mapping transcriptomic to genomic coordinates.
--transcript_fasta_paths_or_urls    No          NA                    Transcript FASTA paths or urls used for mapping transcriptomic to genomic coordinates.
--skip_eventalign_indexing          No          False                 To skip indexing the eventalign nanopolish output.
--genome                            No          False                 To run on Genomic coordinates. Without this argument, the program will run on transcriptomic coordinates.
--n_processes=NUM                   No          1                     Number of processes to run.
--readcount_max=NUM                 No          1000                  Maximum read counts per gene.
--readcount_min=NUM                 No          1                     Minimum read counts per gene.
--resume                            No          False                 With this argument, the program will resume from the previous run.
=================================   ==========  ===================  ============================================================================================================

* Output

======================  ==============  ===============================================================================================================================================================
File name               File type       Description
======================  ==============  ===============================================================================================================================================================
eventalign.index        csv             File index indicating the position in the `eventalign.txt` file (the output of nanopolish eventalign) where the segmentation information of each read index is stored, allowing a random access.
data.json               json            Intensity level mean for each position.
data.index              csv             File index indicating the position in the `data.json` file where the intensity level means across positions of each gene is stored, allowing a random access.
data.log                txt             Gene ids being processed.
data.readcount          csv             Summary of readcounts per gene.
======================  ==============  ===============================================================================================================================================================

``xpore-diffmod``
******************

* Input

Output files from ``xpore-dataprep``.

===================  ==========  ===============      ==============================================================================
Argument name         Required    Default value       Description
===================  ==========  ===============      ==============================================================================
--config=FILE           Yes         NA                YAML configurtaion filepath.
--n_processes=NUM       No          1                 Number of processes to run.
--save_models           No          False             With this argument, the program will save the model parameters for each id.
--resume                No          False             With this argument, the program will resume from the previous run.
--ids=LIST              No          []                Gene / Transcript ids to model.
===================  ==========  ===============      ==============================================================================

* Output

======================  ===============     =================================================================================================================================================
File name                File type           Description
======================  ===============     =================================================================================================================================================
diffmod.table            csv                 Output table information of differential modification rates. Please refer to :ref:`Output table description <outputtable>` for the full description.   
diffmod.log              txt                 Gene/Transcript ids being processed.
======================  ===============     =================================================================================================================================================

``xpore-postprocessing``
**************************

* Input

The ``diffmod.table`` file  from ``xpore-diffmod``.

======================  ===============     =======================================================================
Argument name            Required           Description
======================  ===============     =======================================================================
--diffmod_dir            Yes                Path of the directory containing ``diffmod.table``.
======================  ===============     =======================================================================

