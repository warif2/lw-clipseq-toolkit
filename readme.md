# lw-clipseq-toolkit
A collection of command line python scripts to annotate CLIP
peaks. Additional scripts are provided for further analysis
such as overlapping with rMATS data and simple data summary. 

## Hardware and Software Requirements
  * 64 bit Linux, Mac OS X, or WSL2 on Windows
  * Python 3.8 or higher
  * bedtools suite in PATH
  
## Setting up a virtual environment (Optional)
It is recommended to create a virtual environment before setting
up the script to avoid errors with dependencies. To install Anaconda,
follow the instructions available here, https://docs.anaconda.com/anaconda/install/.

```bash
# Creating venv with the name 'clipseq-tlk' using anaconda
conda create -n clipseq-tlk python=3.8

# Activate venv
conda activate clipseq-tlk
```

## Install required dependencies

```bash
# Download the latest release bedtools using conda
conda install -c bioconda bedtools 

# Install bedops using conda
conda install -c bioconda bedops
```

## Download and Installation

```bash
# Download the latest release of toolkit using git
git clone https://github.com/warif2/lw-clipseq-toolkit.git

# Install python packages dependencies using
pip install -r /lw-clipseq-toolkit/requirements.txt
```

Finally, run setup.py to obtain the license key needed to run the tools.
```bash
# Run setup.py
cd lw-clipseq-toolkit/
python setup.py
```
*Note: enter the 'user' and 'pin' information provided. If you want to use the script and need a 'user' and 'pin', please contact me on waqar.github@gmail.com. 

## Quick Command Line Usage

The following section describes the intended steps in the analysis. The scripts include gff2bed.py, 
annotate_peaks.py, peak_summary.py, peaks2rmats.py, and rmats_filter.py. Proper usage of 
scripts can be understood using the --help option.

Step 1: Convert a gff file to bed format.
```bash
# Example on how to run command gff2bed.py
python gff2bed.py -gff <path_to_gff_file>  -out <desired output name of bed>
```

Once the gff file is converted to a bed file, peak annotation can be performed.
  
Step 2: Run annotate_peaks.py on peak bed file.
```bash
# Example on how to run command annotate_peaks.py
python annotate_peaks.py -pk <path_to_peak_file> -gbed <path_to_gff_bed_file>
```

Refer to the clip_analysis.ipynb python notebook for additional details on script usage and input and output file formats.

## Installing Jupyter and Using Python Notebook
To use the provided clip_analysis.ipynb notebook, jupyter needs to be installed.
```bash
# Install jupyter within the conda environment
# Make sure to have the clipseq-tlk environment activated. 
conda install -c anaconda jupyter

# Launch jupyter notebook
jupyter notebook
```

Open up a web browser and navigate to 'http://localhost:8888/' to view the jupyter console. 
Navigate to the lw-clipseq-toolkit directory and open the clip_analysis.ipynb notebook.