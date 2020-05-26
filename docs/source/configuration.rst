.. _configuration:

Configuration file
==================

The format of configuration file which is one of the inputs for ``xpore-diffmod`` is YAML.

::
    
    data:
        <CONDITION_NAME_1>:
            <REP1>: <DIR_PATH_TO_DATA_JSON>
            ...

        <CONDITION_NMAE_2>:
            <REP1>: <DIR_PATH_TO_DATA_JSON>
            ...

        ...

    paths:
        model_kmer: <PATH_TO_MODEL_KMER>
        out_dir: <DIR_PATH_FOR_OUTPUTS>

    criteria:
        readcount_min: <MINIMUM_READCOUNT>
        readcount_max: <MAXIMUM_READCOUNT>
    
    (OPTIONAL)
    method:
        prefiltering:
            method: t-test
            threshold: <P_VALUE_THRESHOLD>


