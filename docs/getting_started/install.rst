
::::::::::::
Installation
::::::::::::

McClintock is written in Python3 leveraging the `SnakeMake <https://snakemake.readthedocs.io/en/stable/>`_ workflow system and is designed to run on linux operating systems. Installation of software dependencies for McClintock is automated by `Conda <https://docs.conda.io/en/latest/>`_, thus a working installation of Conda is required to install McClintock. Conda can be installed via the `Miniconda installer <https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh>`_.


Installing Miniconda (Python 3.X)
"""""""""""""""""""""""""""""""""

.. code:: bash

   wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O $HOME//miniconda.sh
   bash ~/miniconda.sh -b -p $HOME/miniconda # silent mode
   echo "export PATH=\$PATH:\$HOME/miniconda/bin" >> $HOME/.bashrc # add to .bashrc
   source $HOME/.bashrc
   conda init

.. note:: :code:`conda init` requires you to close and open a new terminal before it takes effect

**Update Conda**

.. code:: bash

   conda update conda

.. seealso:: 
   For more information on Conda, see the `Conda User Guide <https://docs.conda.io/projects/conda/en/latest/index.html>`_

Installing McClintock 
"""""""""""""""""""""

After installing and updating Conda, you should now clone the repository and create the base mcclintock conda environment.

**Clone McClintock Repository**

.. code:: bash

   git clone git@github.com:bergmanlab/mcclintock.git
   cd mcclintock

**Create McClintock Conda Environment**

.. code:: bash

   conda env create -f install/envs/mcclintock.yml --name mcclintock

* This installs the base dependencies (Snakemake, Python3, BioPython) needed to run the main McClintock script into the McClintock Conda environment

**Activate McClintock Conda Environment**

.. code:: bash

   conda activate mcclintock

* This adds the dependencies installed in the McClintock conda environment to the environment PATH so that they can be used by the McClintock scripts.

.. warning:: 
   
   This environment must always be activated prior to running any of the McClintock scripts

.. note:: 
   
   Sometimes activating conda environments does not work via :code:`conda activate myenv` when run through a script submitted to a queueing system, this can be fixed by activating the environment in the script as shown below

   .. code:: bash

      CONDA_BASE=$(conda info --base)
      source ${CONDA_BASE}/etc/profile.d/conda.sh
      conda activate mcclintock


Installing component methods
""""""""""""""""""""""""""""
McClintock utilizes various TE detection tools, each with their own set of dependencies, that are not included in the McClintock repository. These must be installed separately. 

To install all of the component methods and create their associated conda environments, use the following command:

.. code:: bash

   python3 mcclintock.py --install

If you only want to install specific methods to save space and time, you can specify method(s) using the :code:`-m` flag

.. code:: bash

   python3 mcclintock.py --install -m <method1>,<method2>

If you want to install additional methods to an already existing mcclintock installation, you can use the :code:`--resume` flag

.. code:: bash

   python3 mcclintock.py --install -m <method1> --resume

.. note::

   If you do not use the :code:`--resume` flag when installing specific method(s) with `-m`, the installation script will remove all existing method installations and only install what you have specified. For example:

   .. code:: bash

      python3 mcclintock.py --install -m relocate
      python3 mcclintock.py --install -m ngs_te_mapper

   * This command will only leave you with :code:`ngs_te_mapper` installed as :code:`relocate` will be removed during the second installation step. To install :code:`ngs_te_mapper` in addition to :code:`relocate` use:

   .. code:: bash

      python3 mcclintock.py --install -m relocate
      python3 mcclintock.py --install -m ngs_te_mapper --resume


