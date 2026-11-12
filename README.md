# Spatial mediation analysis

Reproduction files for *Disentangling Social and Political Drivers of COVID-19 Mortality Through Spillover Effects in Spatial Mediation Analysis*.

- [Simulations](Simulations/) contains the Gaussian-conditioned simulation, the optional Slurm launcher, the saved results for 500 replications at each of $N=200$ and $N=500$, and the script that reconstructs the two manuscript tables.
- [Data_Analysis](Data_Analysis/) contains the county data, shared distance matrix, Census county boundaries, fitting and processing scripts, saved fits, and figures. The national and four regional fits cover the main analysis, alternative exposure contrast, binary adjacency with distance decay, and omission of mobility.

Run the commands below from the repository root. The simulation and table-processing scripts use base R. For fitting and plotting the application, install the R packages used by these scripts.

    install.packages(c("data.table", "sf", "spdep", "Matrix", "ggplot2", "patchwork"))

## Recreate the reported results from the saved fits

    Rscript Simulations/Process_Results.R
    Rscript Data_Analysis/Process_Results.R
    Rscript Data_Analysis/plot_revision.R

The processing scripts write the rounded values used in the paper to each folder's `Results/Reported` directory. The unrounded simulation metrics, pathway estimates, coefficients, fit diagnostics, and county inclusion records are in `Results`. The plotting script recreates the four effect plots, the national map, and the Q-Q plots in `Data_Analysis/Figures`. That folder also contains the supplied conceptual figures.

## Refit the models

    Rscript Simulations/Spatial_Mediation_MLE_Gaussian_Primary.R
    Rscript Data_Analysis/Main.R

Both scripts use the bundled `Data_Analysis/Data/dist.rds`. By default they write to their respective `Results` directories. To preserve the saved results while rerunning the fits, pass an output directory.

    Rscript Simulations/Spatial_Mediation_MLE_Gaussian_Primary.R Data_Analysis/Data/dist.rds rerun/simulations
    Rscript Data_Analysis/Main.R Data_Analysis/Data Data_Analysis/Data/dist.rds rerun/application

The simulation uses fixed seeds and recreates the two graphs used in the paper. The application can take an optional fourth argument of `USA`, `Midwest`, `Northeast`, `Southeast`, or `West` for a single region. The original application session information is in `Data_Analysis/Results/R_session_info.txt`.

## Application inputs

| File | Use |
| --- | --- |
| `VaccineHesitancyData.csv.gz` | Full source CSV containing county vaccine-hesitancy estimates, SVI, and racial composition. The fitting script reads this losslessly compressed copy. |
| `VaccineHesitancyData_v2.csv` | County Biden vote shares used for the exposure. |
| `trip_mortality_county_weekly_data.csv` | Supplied county-level weekly mobility, population, and mortality input. The mortality analysis selects the records dated June 1, 2021, through February 1, 2022, inclusively. |
| `dist.rds` | County distance matrix used by both analyses. |
| `county_boundaries/` | 2020 Census county cartographic boundary files used to construct adjacency and draw the map. |
| `urban_rural.csv`, `unemployment_wages.csv` | Additional supplied county data files. The reported fits use SVI, racial composition, and mobility for adjustment. |

The vaccine-hesitancy estimates are documented by [CDC and ASPE](https://data.cdc.gov/Vaccinations/Vaccine-Hesitancy-for-COVID-19-County-and-local-es/q9mh-h2tw/about_data). The supplied weekly file is the direct mortality input to the reproduction scripts. Its upstream mortality provenance and aggregation are outside the included source files.
