## Quality Control analysis

#### These scripts generate a report of the quality of your data (currently only single-cell).

All the information goes into the configuration file (YAML format). There is an example (config.yaml) with comments regarding the files' format.

The files you need to prepare are:
	- config.yaml

### Install
Clone this repository (your ~/bin folder is a good place).
```
git clone https://github.com/vijaybioinfo/quality_control.git
cd quality_control
```

### Run the the single-cell QC.
After you've added the necessary information to the YAML file you can call this step like so:
```
Rscript /path/to/quality_control/single_cell.R -y /path/to/project/config.yaml
```
