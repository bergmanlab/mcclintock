
================================
McClintock-formatted predictions
================================

For each TE-prediction component method, McClintock generates a BED file and VCF of the reference and non-reference predictions. This ensures a consistant format for the output which makes it easier to compare the performance of each method. Results in these files are filtered by parameters defined in the :code:`config` files for each method (see :doc:`../running/config` for more details). The McClintock-formatted prediction files are located in :code:`<output>/<sample>/results/<method>/`.

ngs_te_mapper
-------------

:code:`<reference>_ngs_te_mapper_nonredundant.bed`

* BED file containing 0-based intervals corresponding to TSDs for non-reference predictions and 0-based intervals corresponding to the reference TEs. This file contains the same predictions from :code:`unfiltered/<reference>_insertions.bed` with the bed line name adjusted to match the standard McClintock naming convention. By default, no filtering is performed on the raw ngs_te_mapper predictions aside from removing redundant predictions. However, the config file: (:code:`/path/to/mcclintock/config/ngs_te_mapper/ngs_te_mapper_post.py`) can be modified to increase the minimum read support threshold if desired.

RelocaTE
--------

:code:`<reference>_relocate_nonredundant.bed`

* BED file containing predictions from :code:`unfiltered/combined.gff` converted into 0-based intervals with bed line names matching the standard McClintock naming convention. By default, no filtering is performed on the raw predictions aside from removing redundant predictions. However, the config file: (:code:`/path/to/mcclintock/config/relocate/relocate_post.py`) can be modified to increase the minimum left and right prediction support thresholds for both reference and non-reference predictions.

RelocaTE2
---------

:code:`<reference>_relocate2_nonredundant.bed`

* BED file containing all reference and non-reference predictions from :code:`ALL.all_ref_insert.gff` and :code:`ALL.all_nonref_insert.gff`. Coordinates are adjusted to be 0-based. By default, no filtering is performed on split-read and split-pair evidence. However, the config file: (:code:`/path/to/mcclintock/config/relocate2/relocate2_post.py`) can be modified to increase the default threshold for these values.


TEMP
----

:code:`<reference>_temp_nonredundant.bed`

* BED file containing all reference TEs not reported as absent by TEMP in the :code:`unfiltered/<reference>.absence.refined.bp.summary` file. Also contains non-reference TE predictions :code:`unfiltered/<reference>.insertion.refined.bp.summary` formatted as a bed line using the McClintock naming convention. Positions for both reference and non-reference predictions are 0-based. Non-reference predictions from :code:`unfiltered/<reference>.insertion.refined.bp.summary` are only added to this file if the prediction has read support on both ends ("1p1") and has a sample frequency of > 10%. These filtering restrictions can be modified in the config file: (:code:`/path/to/mcclintock/config/TEMP/temp_post.py`). Non-reference TEs with split-read support at both ends are marked in the bed line name with "sr" and the Junction1 and Junction2 columns from :code:`unfiltered/<reference>.insertion.refined.bp.summary` are used as the start and end positions of the TSD in this file (converted to 0-based positions). If the non-reference TE prediction does not have split-read support on both ends of the insertions, the designation "rp" is added to the bed line name and the Start and End columns from :code:`unfiltered/<reference>.insertion.refined.bp.summary` are used as the start and end positions of the TSD in this file (converted to 0-based). Note: TEMP reference insertions are labeled nonab in the bed line name since they are inferred by no evidence of absence to contrast them from reference insertions detected by other components that are inferred from evidence of their presence.

RetroSeq
--------
:code:`<reference>_retroseq_nonredundant.bed`

* BED file containing non-reference TE predictions from :code:`unfiltered/<reference>.call.PE.vcf` with a `Breakpoint confidence threshold <https://github.com/tk2/RetroSeq/wiki/RetroSeq-Tutorial#interpreting-the-output>`_ of >6 are retained in this file. This filtering threshold can be changed by modifying the config file: (:code:`/path/to/mcclintock/config/retroseq/retroseq_post.py`). The position interval reported in the INFO column of :code:`unfiltered/<reference>.call.PE.vcf` are converted to 0-based positions and used as the start and end positions in the bed lines in this file.

PopoolationTE
-------------

:code:`<reference>_popoolationte_nonredundant.bed`

* BED file containing only TE predictions with read support on both ends ("FR") and with percent read support >0.1 for both ends were retained in this file. The entire interval between the inner coordinates of the of the two terminal spans (not the midpoint) was converted to 0-based coordinates. Filtering parameters can be modified in the config file: (:code:`/path/to/mcclintock/config/popoolationte/popoolationte_post.py`)

PopoolationTE2
--------------

:code:`<reference>_popoolationte2_nonredundant.bed`

* BED file containing all of the predictions from :code:`unfiltered/teinsertions.txt` that have support on both ends ("FR") and have a frequency >0.1. The filtering criteria can be modified in the config file: (:code:`/path/to/mcclintock/config/popoolationte2/popoolationte2_post.py`). If predictions overlap a TE in the reference genome, that reference TE is reported in this file using the positions of the reference TE annotation (not the position reported by PoPoolationTE2). If the prediction does not overlap a reference TE, it is designated a non-reference TE insertion :code:`|non-reference|`. The coordinates for all predictions are adjusted to be 0-based.

TE-Locate 
---------

:code:`<reference>_telocate_nonredundant.bed`

* BED file containing all reference and non-reference predictions from :code:`unfiltered/te-locate-raw.info`. Coordinates for both reference and non-reference TE predictions are converted to a 0-based interval. The reference TE end position is extended by the len column in :code:`unfiltered/te-locate-raw.info`. Non-reference TE predictions are a single position as TE-Locate does not predict the TSD size.

TEFLoN
------

:code:`<reference>_teflon_nonredundant.bed`

* BED file containing all reference and non-reference predictions from :code:`unfiltered/genotypes/sample.genotypes.txt`. Reference predictions use the coordinates for the TE with the reference ID from column 7. By default, only non-reference predictions with both breakpoints (C2 and C3) are kept in this file. Non-reference predictions must also have at least 3 presence reads (C10) and an allele frequency greater than 0.1 (C13). These filtering restrictions can be changed by modifying the TEFLoN config file: :code:`/path/to/mcclintock/config/teflon/teflon_post.py`