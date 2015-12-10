YeastPhenome.org Dataset README File
====================================

The Basics
----------

The current folder contains:

1. the raw data corresponding to this publication (`/raw_data/*`): usually, text or Excel files;
2. the filtered, edited & re-formatted data (`*.mat`): Matlab structure;
3. the code used to transform #1 to #2 (`code.m`): Matlab function.

Running `code.m` will regenerate the `*.mat` file from the `raw_data/` and output the list of files that were used in the process:

    filenames = code();


The Details
-----------

The raw data files (#1) are identical to the ones provided by the authors in the original publication or in their personal communication with the curators of Yeastphenome.org. No edits to the content or even the name of the file were performed. The only exception is when the raw data was derived from a PDF (main text or supplement): in that case, the data were extracted from PDF and saved as a text file.

The raw data files are used as input to the `code.m` function which performs a number of standard checks and converts the data into a common format. 

The list of checks includes:

* naming conversions (common names -> systematic ORF names)
* naming errors (typos in the ORF names)
* unexpected values ("N/A", "blank", "empty", etc. -> "NaN")

The list of conversions includes:

* normalization to wild-type or control, whenever the values for experiment and control are provided separately in the raw data
* reversing the order of the values (by taking the opposite or the reciprocal values) such that "higher" and "lower" values correspond to the chosen description of the phenotype. For example, if the data were reported as "sensitivity to drug X", we describe the phenotype to "growth in the presence of drug X" and convert the values such that higher sensitivity corresponds to lower growth.

Help
----

Please direct questions & comments to Anastasia Baryshnikova (<abarysh@princeton.edu>) or check in regularly at <http://www.yeastphenome.org>.

