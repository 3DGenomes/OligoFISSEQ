# OligoFISSEQ

## Analysis Walkthrough

### Requirements

<p>The walkthrough is based on [Jupyter notebooks](https://jupyter.org/install). The rest of the requirements are basically the ones required by the [ImageJ tutorials](https://github.com/imagej/tutorials) plus the [Integrative Modeling Platform](https://integrativemodeling.org/) and the python implementation of the [Constrained K-means](https://github.com/Behrouz-Babaki/COP-Kmeans) algorithm are listed in the environment.yml file.</p>
<p>The easiest to install them in an isolated environment in to use [Miniconda](https://conda.io/miniconda.html). A way to install it could be:</p>

    wget --quiet --no-check-certificate https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh \
         -O ~/miniconda.sh && \
         /bin/bash ~/miniconda.sh -b && \
         rm ~/miniconda.sh

<p>Then open follow this steps:</p>

1. Clone this repository.
2. Open a console and `cd` to your cloned working copy.
3. `conda env create -f environment.yml` to create a conda environment with the
   dependencies the notebook need.
4. `conda activate scijava` to activate the environment.
5. `jupyter notebook` to launch Jupyter Notebook in a web browser.
6. In the browser, open the `OligoFISSEQ pipeline.ipynb`.

### The OligoFISSEQ pipeline notebook

<p>In the notebook there will be an example of the analysis workflow of one image produced by OligoFISSEQ. All the files needed are contained in this repository</p>   

## Scripts

### Preprocessing

**create_hyperstack_from_rounds.ijm**
<p>ImageJ script to compile each round of OligoFISSEQ into hyperstacks composed by 5 channels, a series of z-slices and one frame per round. If an image where all the punctas are labelled, like in toto image, is available it is included as a new additional frame. The hyperstacks are aligned using Fiji plugin “Correct 3D Drift” (Parslow, Cardona, and Bryson-Richardson 2014). Images of DAPI stained nuclei are used to perform threshold segmentation and extract each individual cell from the initial image as a separate region of interest (roi).</p>

### Tier 1 detection of barcodes

**detect_barcode_patches.py**
<p>ImageJ jython script to detect barcodes in the hyperstacks created in the preprocessing. For each hyperstack this script produces:</p>

* A tiff file with the maximum intensity projection over all rounds over all channels
* A tiff file where pixels yielding a full barcode are colored by barcode
* A tiff file where each slice contains the signal of each barcode individually
* A tsv file with information of every detected patch with full and subsampled barcoding

### Tier 2 Chromosome tracing

**trace_OligoFISSEQ.py**
<p>Python script to trace chromosomes for the 36plex datasets using the tsv file with the detected barcodes. The output of this script is:</p>

* A tsv file with the positions of the traced loci

**trace_OligoFISSEQ_chrX.py**
<p>Python script to trace chromosomes for the chrX-46plex datasets using the tsv file with the detected barcodes. The output of this script is:</p>

* A tsv file with the positions of the traced loci


