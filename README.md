Non-disruptive methylation monitoring of cellular states with cell-free DNA 
===========

Author: Anja Hess

Affiliation: Max Planck Institute for Molecular Genetics, Berlin, Germany

Date: 2025-JULY-21

## 1. Description
Scripts to reproduce figures from the manuscript:
Hess, Anja et al. 2025: Non-disruptive methylation monitoring of cellular states with cell-free DNA

## 2. Before you start

Please make sure you have installed **python => 3.12** and **R =>4.4.2**. Next, download the required packages:

For Python-based scripts (majority):

    pip3 install numpy pandas seaborn matplotlib 

For R-based scripts:

    R
    install.packages(c("ggplot2","RColorBrewer", "pheatmap", "circlize", "vioplot"))


The majority of the scripts were run on a 32 GB memory local PC on Ubuntu 24.04.2 LTS. Scripts using larger sequencing data 
as an input (e.g. BAM data such as in Fig. 1E) were run in a **Unix** environment on a remote server (required memory **~ 20-40 GB**).


## 3. Package architecture:

### 3.1 Main (./)
    The main folder contains 
    - this README
    - scripts for each figure, named Figure_FigNr_SubfigureNr.py or Figure_FigNr_SubfigureNr.R
    - R_plots.R, which is a utility to store simple plotting functions frequently used across figures
    - plots.py, serving a similar function as R_plots.R, but for python-based scripts

### 3.2 Source data (./sourcedata)

This folder contains subdirectories for all main figures and the connected panels. Whenever compatible with Github's
file size limitations, the source data is provided to ease reproducing the results. Importantly, resulting plots from execution of the figure scripts will be saved **into their figure's folder**.

    │── FIG_1
    │   ├── 1B
    │   ├── 1CD
    │   ├── 1E
    │   ├── 1F
    │   └── 1G
    ├── FIG_1EXT
    │   ├── E1A
    │   ├── E1B
    │   ├── E1C
    │   ├── E1D
    │   ├── E1E
    │   ├── E1F
    │   └── E1G
    ├── FIG_2
    │   ├── 2B
    │   ├── 2CD
    │   └── 2E
    └── FIG_2EXT
        ├── E2A
        ├── E2C
        ├── E2D
        └── E2E-F


### 3.3 A word on large source data 

There are some instances where the scripts take sourcedata as input that are **too large** to be stored within this Github repository. In such cases, we have uploaded the underlying files at the Gene Expression Omnibus (GEO) 
under the accession GSE293866. You will find detailed instructions on how to download and process the files in the **README** stored in the figure's subdirectory.

    ├── 1F
    │   ├── README

In the example above, the .h5ad file underlying Figure 1F can be accessed following the instructions in its folder's README.

## 4. Usage example

First, make sure you follow the **README** to download data not provided with this repository. Make sure you put them into the **correct folder** matching the figure number.

    cd sourcedata/FIG_1/1F/
    wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE293866&format=file&file=WHATEVER_FILE_ID
    ll
    3115144544 Dez  8  2024 WGBS_PCA_DNMT1i_mESC.h5ad

The file is now in the correct location, so we can continue. 
To reproduce a figure of choice simply navigate to the main directory and execute the corresponding script. The example below shows you how to reproduce the principal component analysis (PCA) plot for Figure 1F.

    cd liquidbio/
    python3 Figure_1_F.py 

This will result in the following output.

    WARNING: saving figure to file sourcedata/FIG_1/1F/pcamean-meth.pdf
    UserWarning: FigureCanvasAgg is non-interactive, and thus cannot be shown
      plt.show()

Now navigate to liquidbio/sourcedata/FIG_1/1F/ to find the newly created plot pcamean-meth.pdf

![Test Image 5](./sourcedata/FIG_1/1F/pcamean-meth.png)

## 5. Data access

Raw and processed sequencing data generated in this study are publicly available on GEO under the accession GSE293866.

## 6. Citation

**Anja Hess, Alexander Kovacsovics, Fabian Bachinger, Helene Kretzmer, Ludovic Vallier and Alexander Meissner: 
Non-disruptive methylation monitoring of cellular states with cell-free DNA (2025)**

*Fabian Bachinger and  Alexander Kovacsovics contributed equally as second authors to this work.*