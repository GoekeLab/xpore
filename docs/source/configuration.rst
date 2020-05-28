.. _configuration:

Configuration file
==================

The format of configuration file which is one of the inputs for ``xpore-diffmod`` is YAML.

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
    
    (OPTIONAL)
    method:
        prefiltering:
            method: t-test
            threshold: <P_VALUE_THRESHOLD>


