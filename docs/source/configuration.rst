.. _configuration:

Configuration file
==================

The format of configuration file which is one of the inputs for ``xpore-diffmod`` is YAML.

::
    
    data:
        <DATASET_NAME_1>:
            condition_name: <CONDITION_NAME>
            dirpath: <DIR_PATH_TO_DATA_JSON>

        <DATASET_NMAE_2>:
            condition_name: <CONDITION_NAME>
            dirpath: <DIR_PATH_TO_DATA_JSON>

    paths:
        model_kmer: <PATH_TO_MODEL_KMER>
        out_dir: <DIR_PATH_FOR_OUTPUTS>

    criteria:
        readcount_min: <MINIMUM_READCOUNT>
        readcount_max: <MAXIMUM_READCOUNT>

    method:
        prefiltering:
            method: t-test
            threshold: <P_VALUE_THRESHOLD>


