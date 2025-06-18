# Compute CPMs from Gasperini et al. 2019

**Author:** Eugenio Mattei  
**Date:** June 12, 2025

## Overview

Compute counts per million (CPM) values for genes from the Gasperini et al. 2019 dataset.
This analysis processes the K562 single-cell raw counts data to generate normalized CPM values for each gene.


## How to Run

1. **Install dependencies**
    ```R
    # Install the 'renv' R package if not already installed
    if (!requireNamespace("renv", quietly = TRUE)) {
      install.packages("renv")
    }
    # Restore the project environment
    renv::restore()
    ```

2. **Run the analysis**  
   ```sh
   Rscript -e "rmarkdown::render('compute_cpms_from_gasperini.Rmd', output_file='compute_cpms_from_gasperini.html')"
   ```

## Outputs

These files are considered final and will be used in subsequent analyses. All the intermediates files are deleted.

### HTML summary
The R markdown notebook knitted into an HTML file.
[compute_cpms_from_gasperini.html]("./compute_cpms_from_gasperini.html")

### Per gene CPM values
The final CPM values for each gene can be found here: [outputs/gasperini_cpms.csv](outputs/gasperini_cpms.csv).

*The file `outputs/gasperini_cpms.csv` is comma-separated and includes a header row:*
```
gencode_v26lift37_id,cpm
ENSG00000123,1050.4
ENSG00000456,987.2
ENSG00000789,876.5
ENSG00000234,765.1
ENSG00000567,654.3
```



## Dependencies

- R >= 4.0

## Notes

- The `monocle` R package is being installed otherwise the RDS object cannot be read successfully.

## Contact

For questions, please open an issue on GitHub.