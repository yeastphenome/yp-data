YeastPhenome.org Dataset README File
====================================

Warning
-------

THESE DATA MAY CONTAIN ERRORS. YeastPhenome.org is running in beta version. The data are available for download, but we cannot guarantee the absence of errors or code bugs introduced during processing. This warning will be removed after all cross-checks and validations have been completed. In the meantime, please, be careful when using the data.


Citations
---------

If you use yeast phenome data in a talk or a manuscript, please acknowledge the data source by citing the original publication and the yeast phenome database.
The database should be cited as follows:

> Database of Systematic Phenotypic Screens in Yeast, www.yeastphenome.org, accessed on XXXX-XX-XX.


The Basics
----------

This folder contains:

1. the filtered, edited & re-formatted data from this publication (`.txt`): a tab-delimited matrix in TXT format;
2. the filtered, edited & re-formatted data from this publication (`.mat`): a Matlab structure;
3. the raw data from which #1 and #2 were derived (`/raw_data/*`): in most cases, 1 or more text or Excel files;
4. the code used to transform #1 into #2 (`code.m`): a Matlab function.

Running `code.m` will regenerate the `.mat` and the `.txt` file from the `raw_data/*` files and output the list of files that were used in the process:

    filenames = code();


The Details
-----------

The raw data files (#1) are identical to what was provided by the authors in the original publication or in their personal communication with the curators of Yeastphenome.org. No edits to the content or even the name of the file were performed, with the exception of:

* PDF -> text file conversions
* XLS -> XLSX conversions (ensures a more stable and consistent loading process)

The raw data files are used as input to the `code.m` function which performs a number of standard checks and converts the data into a common format. 

Some of the checks include:

* naming conversions (common names -> systematic ORF names)
* naming errors (typos, spaces or capitalization errors in the ORF names)
* unexpected values ("N/A", "blank", "empty", etc. -> "NaN")

The list of conversions includes:

* normalization to wild-type or control, whenever the values for experiment and control are provided separately in the raw data
* reversing the order of the values (by taking the opposite or the reciprocal values) such that "higher" and "lower" values correspond to the chosen description of the phenotype. For example, if the data were reported as "sensitivity to drug X", we describe the phenotype to "growth in the presence of drug X" and convert the values such that higher sensitivity corresponds to lower growth.

Help
----

Please direct questions & comments to Anastasia Baryshnikova (<abaryshnikova@calicolabs.com>) or check in regularly at <http://www.yeastphenome.org>.
