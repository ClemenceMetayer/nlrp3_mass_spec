### Install necessary packages

[CosinorPy](https://github.com/mmoskon/CosinorPy)
[RAIN](https://www.bioconductor.org/packages/release/bioc/html/rain.html)

### Run `create_dataset.py`

This will create the data structure I use from the MS data.

### Run `run_cosinor_list_period.py`

This will save files with the cosinor results and best models founds.
This file uses the updated version of CosinorPy (cosinor_float_period.py) which allows decimal periods. 
The function get_best_fits has also been modified to determine the best period with a F-test.

### Have a look at the notebook `cosinor_analysis_list_period.ipynb`

This notebook presents the results of the Cosinor analysis for the range of period (20h-30h). 
