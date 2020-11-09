## Quality Control analysis

#### These scripts generate a report of the quality of your data (currently only single-cell).

All the information goes into the configuration file (YAML format). There is an example (config.yaml) with comments regarding the files' format.

The files you need to prepare are:
- config.yaml
- metadata (optional): if you wish to add information describing each library in a multiple library project. Check example in the [example](./data/metadata_library.csv) in the data folder.
- aggregation (optional): if you're dealing with several libraries you haven't aggregated or don't know yet if you will. Check example in the [example](./data/aggregation.csv) in the data folder.

### Install
Clone this repository (your ~/bin folder is a good place).
```
git clone https://github.com/vijaybioinfo/quality_control.git
```

### Run the the single-cell QC.
After you've added the necessary information to the YAML file you can call this step like so:
```
Rscript /path/to/quality_control/single_cell.R -y /path/to/project/config.yaml
```
