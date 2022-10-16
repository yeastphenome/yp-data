YeastPhenome-data README File
====================================

Quick intro
------------

This repository contains the code used to clean, transform and normalize the data included in the YeastPhenome.org project. The data and code are organized into folders (1 folder per publication, identified by its Pubmed ID). Each folder contains the raw data (i.e., data obtained from the publication or any other sources), a Jupyter notebook containing code for processing, and the transformed data.

For more detailed information about the sources of data and the description of processing/transformation steps, see www.yeastphenome.org and the corresponding publication (coming soon).


Installation
-------------

To run the processing code (`code.ipynb`), create the following environment:

    conda create -n yp-data python=3.7.7 --file requirements.txt
    conda activate yp-data
    ipython kernel install --user --name=yp-data


To clone the repository:

    git clone git@github.com:yeastphenome/yp-data.git

To download only the code/data for a specific publication:

    mkdir yp-data
    cd yp-data
    mkdir Utils
    mkdir Datasets
    cd Datasets
    mkdir <PMID>

Download all contents of the Github `Utils/` folder into the local `Utils/` and all contents of the Github `<PMID>/` folder into the local `<PMID>/` (<PMID> is the Pubmed ID of the publication of interest).


Help
----

Please direct questions & comments to Anastasia Baryshnikova (<abaryshnikova@calicolabs.com>).

