.. _scripts:

Scripts
==========

We provide 2 main scripts to run the analysis of differential RNA modifications as the following.

``xpore-dataprep``
********************

* Input

Output files from `nanopolish eventalgin`. Please refer to :ref:`Data preparation <preparation>` for the full Nanopolish command.

* Usage example

=================================   ==========  ===================  ============================================================================================================
Argument name(s)                    Required    Default value         Description
=================================   ==========  ===================  ============================================================================================================
--eventalign=FILE                   Yes         NA                    Eventalign filepath, the output from nanopolish.         
--summary=FILE                      Yes         NA                    Eventalign summary filepath, the output from nanopolish.
--out_dir=DIR                       Yes         NA                    Output directory.
--ensembl=NUM                       No          91                    Ensembl version for gene-transcript mapping.
--species=STR                       No          homo_sapiens          Species for ensembl gene-transcript mapping.
--customised_genome                 No          False                 If customised genome provided.
--reference_name                    No          NA                    fasta reference name.
--annotation_name                   No          NA                    gtf annotation name.
--gtf_path_or_url                   No          NA                    gtf file path or url.
--transcript_fasta_paths_or_urls    No          NA                    Transcript fasta paths or urls.
--skip_eventalign_indexing          No          False                 To skip indexing the eventalign nanopolish output.
--genome                            No          False                 To run on Genomic coordinates. Without this argument, the program will run on transcriptomic coordinates.
--n_processes=NUM                   No          1                     Number of processes to run.
--readcount_max=NUM                 No          1000                  Maximum read counts per gene.
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

Output files from `xpore-dataprep`.

* Usage example

===================  ==========  ===============      ==============================================================================
Argument name(s)      Required    Default value       Description
===================  ==========  ===============      ==============================================================================
--config=FILE           Yes         NA                Yaml configuraion filepath.
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
   
