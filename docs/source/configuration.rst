.. _configuration:

Configuration file
==================

The format of configuration file which is one of the inputs for ``xpore-diffmod`` is YAML.

Only the ``data`` and ``out`` sections are required, other sections are optional. Below is the detail for each section.

::
    
    data:
        <CONDITION_NAME_1>:
            <REP1>: <DIR_PATH_TO_DATA_JSON>
            ...

        <CONDITION_NAME_2>:
            <REP1>: <DIR_PATH_TO_DATA_JSON>
            ...

        ...

    out: <DIR_PATH_FOR_OUTPUTS>
    
    criteria:
        readcount_min: <15>
        readcount_max: <1000>
        
    method:
        # To speed up xpore-diffmod, you can use a statistical test (currently only t-test is implemented) can be used 
        # to remove positions that are unlikely to be differentially modified. So, xpore-diffmod will model only 
        # those significant positions by the statistical test -- usually the P_VALUE_THRESHOLD very high e.g. 0.1. 
        # If you want xPore to test every genomic/transcriptomic position, please remove this prefiltering section.
        prefiltering:
            method: t-test
            threshold: <P_VALUE_THRESHOLD>
        
        # Here are the parameters for Bayesian inference. The default values shown in <> are used, if not specified. 
        max_iters: <500>
        stopping_criteria: <0.00001>
        



