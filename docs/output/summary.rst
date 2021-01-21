
==================
Summary text files
==================

McClintock provides a basic summary of the results of the component methods in a fer text files. These files show the number of predictions for each method, number of predictions for each TE family, as well as a summmary of the coverage for each consensus TE or TE provided in coverage fasta.

:code:`<output>/<sample>/results/summary/data/run/summary_report.txt`

* Summary Report of McClintock run. Contains information on the McClintock command used, when and where the script was run, details about the mapped reads, and table that shows the number of TE predictions produced from each method.

:code:`<output>/<sample>/results/summary/data/run/te_prediction_summary.txt`

* A comma-delimited table showing reference and non-reference predictions for each component method

:code:`<output>/<sample>/results/summary/data/families/family_prediction_summary.txt`

* a comma-delimited table showing TE predictions (all, reference, non-reference) from each method for each TE family

:code:`<output>/<sample>/results/summary/data/coverage/te_depth.txt`

* (Only produced if coverage module is run) a comma-delimited table showing normalized depth for each consensus TE or TE provided in coverage fasta All tables and plots contain a link to the raw data so that users can manually filter or visualize it with other programs.
