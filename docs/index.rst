.. mcclintock documentation master file, created by
   sphinx-quickstart on Tue Jan  5 10:08:16 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==========
McClintock
==========

Many methods have been developed to detect transposable element (TE) insertions from whole genome shotgun next-generation sequencing (NGS) data, each of which has different dependencies, run interfaces, and output formats. Here, we have developed a meta-pipeline to install and run multiple methods for detecting TE insertions in NGS data, which generates output in the UCSC Browser extensible data (BED) format. A detailed description of the original McClintock pipeline and evaluation of the original six McClintock component methods on the yeast genome can be found in Nelson, Linheiro and Bergman (2017) G3 7:2763-2778.

---------------------------------
TE Dectection Software Components
---------------------------------

.. csv-table::
   :header: "Method", "Citation"
   :widths: 10, 15

   `ngs_te_mapper <https://github.com/bergmanlab/ngs_te_mapper>`_, `Linheiro and Bergman (2012) <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030008>`_
   `RelocaTE <https://github.com/srobb1/RelocaTE>`_, `Robb et al. (2013) <http://www.g3journal.org/content/3/6/949.long>`_
   `RelocaTE2 <https://github.com/stajichlab/RelocaTE2>`_, `Chen et al. (2017) <https://peerj.com/articles/2942/>`_
   `TEMP <https://github.com/JialiUMassWengLab/TEMP>`_, `Zhuang et al. (2014) <http://nar.oxfordjournals.org/content/42/11/6826.full>`_
   `RetroSeq <https://github.com/tk2/RetroSeq>`_, `Keane et al. (2012) <http://bioinformatics.oxfordjournals.org/content/29/3/389.long>`_
   `PoPoolationTE <https://sourceforge.net/projects/popoolationte/>`_, `Kofler et al. (2012) <https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002487>`_
   `PoPoolationTE2 <https://sourceforge.net/p/popoolation-te2/wiki/Home>`_, `Kofler et al. (2016) <https://academic.oup.com/mbe/article/33/10/2759/2925581>`_
   `TE-locate <https://sourceforge.net/projects/te-locate/>`_, `Platzer et al. (2012) <http://www.mdpi.com/2079-7737/1/2/395>`_
   `TEFLoN <https://github.com/jradrion/TEFLoN>`_, `Adrion et al. (2017) <https://academic.oup.com/gbe/article/9/5/1329/3064433>`_


.. toctree::
   :maxdepth: 1
   :caption: Getting Started
   :hidden:

   getting_started/install


