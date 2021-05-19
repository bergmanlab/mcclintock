
###################
Method config files
###################

Each component method has associated config file(s) that enable the modification of run parameters and filtering parameters. TE detection components have separate run and post processing config files. These files are located in the `config directory <https://github.com/bergmanlab/mcclintock/tree/master/config>`_ of mcclintock :code:`/path/to/mcclintock/config`.

***********************
Example run config file
***********************
.. code:: python

    PARAMS = {
        '-l' : 10,
        '-m' : 0.0,
        '-bm' : 10,
        '-bt' : 7,
        '-f' : 100
    }

Each config file contains a python :code:`dict` named :code:`PARAMS` which contains :code:`key : value` pairs with the :code:`key` being the flag or name of the parameter, and the :code:`value` being the value associated with the parameter. Most of the :code:`values` are set to the default parameters for each tool. The :code:`values` can be modified to better fit the data being used.

********
Coverage
********
:code:`/path/to/mcclintock/config/coverage/coverage.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/coverage/coverage.py>`_

.. code:: python

    PARAMS = {
        "omit_edges": True,
        "omit_edges_read_length" : True,
        "omit_edges_length" : 300
    }

This config file contains the parameters that can be modified for the :code:`coverage` component method.

omit_edges
  * Omits the edges of the TE sequence coverage when calculating average depth across the element.

omit_edges_read_length
  * If :code:`omit_edges: True` and :code:`omit_edges_read_length True`, then the read length will be used as the length of the edges to omit

omit_edges_length
  * If :code:`omit_edges: True` and :code:`omit_edges_read_length False`, the value of :code:`omit_edges_length` will be used as the length of the edges to omit

*************
ngs_te_mapper
*************

run config
==========
:code:`/path/to/mcclintock/config/ngs_te_mapper/ngs_te_mapper_run.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/ngs_te_mapper/ngs_te_mapper_run.py>`_

.. code:: python

    PARAMS = {
        "tsd=" : 20
    }

This config file contains the parameters that can be modified to influence the running of the :code:`ngs_te_mapper` component method.

tsd=
  * The max length of junction read overlap to consider a target site duplication

postprocessing config
=====================
:code:`/path/to/mcclintock/config/ngs_te_mapper/ngs_te_mapper_post.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/ngs_te_mapper/ngs_te_mapper_post.py>`_

.. code:: python

    PARAMS = {
        "min_read_support" : 0
    }

This config file contains the parameters that influence the post processing of the :code:`ngs_te_mapper` predictions.

min_read_support
  * the minimum number of reads supporting an insertion required to be considered a valid prediction