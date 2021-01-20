
======================
Raw method predictions
======================

Each component method produces files that contain the TE predictions and the level of support for each TE prediction. These files are used directly to create the standardized BED and VCF output (see See :doc:`formatted` for more details).
The raw output from component methods can be found in :code:`<output>/<sample>/results/<method>/unfiltered/`

ngs_te_mapper
-------------

:code:`<reference>_insertions.bed`

* BED file containing raw 0-based intervals corresponding to TSDs for non-reference predictions and 0-based intervals corresponding to the reference TEs. Reference TE intervals are inferred from the data, not from the reference TE annotations. Strand information is present for both non-reference and reference TEs.

RelocaTE
--------

:code:`combined.gff`

* GFF containing 1-based TSDs for non-reference predictions and 1-based intervals for reference TEs. The reference intervals are based on the reference TE annotations.

RelocaTE2
---------

:code:`repeat/results/ALL.all_ref_insert.gff`

* GFF file containing reference TE predictions with 1-based coordinates. The final column also contains read counts supporting the junction (split-read) and read counts supporting the insertion (read pair).

:code:`unfiltered/repeat/results/ALL.all_nonref_insert.gff`

* GFF file containing non-reference TE predictions with 1-based coordinates. The final column also contains read counts supporting the junction (split-read) and read counts supporting the insertion (read pair).

TEMP
----

:code:`<reference>.absence.refined.bp.summary`

* Tab-delimited table containing reference TEs that are predicted to be absent from the short read data. Position intervals are 1-based.

:code:`<reference>.insertion.refined.bp.summary`

* Tab-delimited table containing non-reference TE predictions. Position intervals are 1-based.

RetroSeq
--------

:code:`<reference>.call.PE.vcf`

* VCF file containing non-reference TE predictions. Non-reference TEs are annotated as 1-based intervals in the POS column and two consecutive coordinates in the INFO field. No predictions are made for reference TEs. Strand information is not provided.

PopoolationTE
-------------

:code:`te-poly-filtered.txt`

* Tab-delimited table with non-reference and reference TE predictions and support values. Predictions are annotated as 1-based intervals on either end of the predicted insertion, and also as a midpoint between the inner coordinates of the two terminal spans (which can lead to half-base midpoint coordinates)

PopoolationTE2
--------------

:code:`teinsertions.txt`

* Tab-delimited table with TE predictions and TE frequency values (ratio of physical coverage supporting a TE insertion to the total physical coverage). PoPoolationTE2 does not indicate which predictions are reference and non-reference TEs. Also, only a single position is reported for each prediction, so the TSD is not predicted. Predictions may only have support from one side of the junction ("F" or "R") or both sides ("FR"). Prediction coordinates are 1-based.

TE-Locate 
---------

:code:`te-locate-raw.info`

* A tab-delimited table containing reference ("old") and non-reference ("new") predictions using 1-based positions. TSD intervals are not predicted for non-reference TEs, instead a single position is reported.

TEFLoN
------

:code:`genotypes/sample.genotypes.txt`

* A tab-delimited table containing all of the breakpoints and support information for insertion predictions. Predictions are treated as reference predictions if they contain a TE ID in column 7.

.. code:: text

    # from: https://github.com/jradrion/TEFLoN
    C1: chromosome
    C2: 5' breakpoint estimate ("-" if estimate not available)
    C3: 3' breakpoint estimate ("-" if estimate not available)
    C4: search level id (Usually TE family)
    C5: cluster level id (Usually TE order or class)
    C6: strand ("." if strand could not be detected)
    C7: reference TE ID ("-" if novel insertion)
    C8: 5' breakpoint is supported by soft-clipped reads (if TRUE "+" else "-")
    C9: 3' breakpoint is supported by soft-clipped reads (if TRUE "+" else "-")
    C10: read count for "presence reads"
    C11: read count for "absence reads"
    C12: read count for "ambiguous reads"
    C13: genotype for every TE (allele frequency for pooled data, present/absent for haploid, present/absent/heterozygous for diploid) #Note: haploid/diploid caller is under construction, use "pooled" for presence/absence read counts
    C14: numbered identifier for each TE in the population