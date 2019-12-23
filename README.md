# OligoFISSEQ

## Tier 1 detection

### Preprocessing

create_hyperstack_from_rounds.ijm : ImageJ script to compile each round of OligoFISSEQ into hyperstacks composed by 5 channels, a series of z-slices and one frame per round. If an image where all the punctas are labelled, like in toto image, is
available it is included as a new additional frame. The hyperstacks are aligned using Fiji plugin “Correct 3D Drift” (Parslow, Cardona, and Bryson-Richardson 2014). Images of DAPI stained nuclei are used to perform threshold segmentation and extract each individual cell from the initial image as a separate region of interest (roi).

### Detection of barcodes

detect_barcode_patches.py: ImageJ jython script to detect barcodes in the hyperstacks created in the preprocessing. For each hyperstack this script produces:

* A tiff file with the maximum intensity projection over all rounds over all channels
* A tiff file where pixels yielding a full barcode are colored by barcode
* A tiff file where each slice contains the signal of each barcode individually
* A tsv file with information of every detected patch with full and subsampled barcoding

## Tier 2

### Chromosome tracing

domino_oligo_cluster.py: Python script to trace chromosomes for the 36plex datasets using the tsv file with the detected barcodes. The output of this script is:

* A tsv file with the positions of the traced loci

domino_oligo_cluster_chrX.py: Python script to trace chromosomes for the chrX-46plex datasets using the tsv file with the detected barcodes. The output of this script is:

* A tsv file with the positions of the traced loci


