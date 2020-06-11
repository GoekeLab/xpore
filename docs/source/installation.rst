.. _installation:

Installation
=======================

xPore requires `Python3 <https://www.python.org>`_ to run.

PyPI installation (recommended)
---------------------------------
::

    pip install xpore
    pyensembl install --release 91 --species homo_sapiens  # please specify the compatible Ensembl release with your data when you install it.

Installation from our GitHub repository
---------------------------------------
::

    git clone https://github.com/GoekeLab/xpore.git
    cd xpore
    python setup.py install
    pyensembl install --release 91 --species homo_sapiens  # please specify the compatible Ensembl release with your data when you install it.
