# Phage growth curves analyses

## Installation

Create virtual environment. If Conda is not yet present, install [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html).

```bash
conda create -n phageGC python=3.8.17
conda activate phageGC
```

Install dependencies

```bash
pip install -r requirements.txt
```

## Calculate AUC of growth curves with given cut-off timepoints

**Input:** an Excel .xlsx file with at least two sheets:

1. Sheet with data describing layouts of the culture plates and cut-off timepoints for AUC calculation. Sheet name: `Schema`.

2. Sheet(s) with data of OD measurements. Data of different plates (replicates) should be in separate sheets. Sheet names should correspond to values in the `Replicate` column in the `Schema` sheet. 

**Output:** a .CSV file with AUC data for all sample wells defined in the `Schema` sheet of the input file.

- (default) AUC values for each sample well
- (with option `-a`, or `--average`) Averaged AUC values across replicates (culture plates) for each sample type

**Basic command**:

```bash
conda activate phageGC

python phage_auc.py -f <input_Excel_file> -o <output_filename.csv>
```

**Full description:**

```
usage: phage_auc.py [-h] -f DATA_FILE [-o OUT_FILE] [-a] [-s SCHEMA_SHEETNAME]
                    [-ts T_START] [-te T_END] [-ti T_INTERVAL] [-i IGNORE_BACTERIA]

Calculate AUC of growth curves with given cut-off timepoints.

optional arguments:
  -h, --help            show this help message and exit
  -f DATA_FILE, --data_file DATA_FILE
                        Excel input file [required]
  -o OUT_FILE, --out_file OUT_FILE
                        CSV output file [required]
  -a, --average         Enable averaging AUC values across replicates if specified.
  -s SCHEMA_SHEETNAME, --schema_sheetname SCHEMA_SHEETNAME
                        Name of the Excel sheet for the schema of culture plates.
                        [default: 'Schema']
  -ts T_START, --t_start T_START
                        Start time of experiments, in hour. [default: 0.0]
  -te T_END, --t_end T_END
                        End time of experiments, in hour. [default: 24.0]
  -ti T_INTERVAL, --t_interval T_INTERVAL
                        Time interval between every 2 consecutive OD measurements, in
                        hour. [default: 1/12 (5 min)]
  -i IGNORE_BACTERIA, --ignore_bacteria IGNORE_BACTERIA
                        ID of bacteria (sp1_ID) to ignore in calculations. Use multiple
                        times when there're >1 bacteria to ignore. [default: 'KP85']
```
