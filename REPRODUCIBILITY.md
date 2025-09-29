# How to Reproduce the Analyses

Welcome! In here we collect all the documentation related to the reproducibility of the analyses.

The use everything more conveniently and to leverage the project folder structure (on which the scripts rely on), it is recommended to download this repository. To do so, head to the [Repo's Homepage]() and click on the `Code` button.

If you don't use git or if you are unfamiliar with it, select `Download ZIP`.
If you use git you can either copy the _https_ and _ssh_ code and use `git clone` or pasting this directly in your terminal

```
git clone https://github.com/
```

## Requirements

R must be installed on your system. The scripts use a variety of packages, which will be installed by {renv} or by the scripts (see below).
RStudio or any "R-friendly" IDE will make everything easier (the analyses have been originally run in [Positron](https://positron.posit.co/)).


## Reproducibility environment

This project leverages `{renv}` to ensure exact reproducibility. If you are not familiar with `{renv}`, take a look at the main [documentation](https://rstudio.github.io/renv/articles/renv.html). Put simply, renv (through the renv.lock file) records the packages and their versions used during the analyses, and it will try to reinstall these versions on your machine (inside the project folder). In this way the analyses will run under the same development environment.

> [!NOTE]
> you don't have to use {renv}, you can continue using your own version of packages, we will make sure to update and check the scripts regularly to ensure compatibility with newer versions.

In the repository you will find a `use_renv.txt` file.
If this file is present, {renv} will be initialized (installed if not present in your library) and will ask you if you want to download or update the packages listed in the renv.lock file.

> [!TIP]
> In case you don't want to use {renv}, remove the `use_renv.txt` file BEFORE opening the project, or rename it.

Given the wide variety of operating systems and environments we cannot guarantee that all the packages will install correctly (some might need additional system dependencies). In case you activated {renv} but run into troubles, follow the instructions to [deactivate or uninstall it](https://rstudio.github.io/renv/articles/renv.html#uninstalling-renv), or more simply run in the R console:

```R
renv::deactivate(clean = TRUE) # this deletes the renv folder and renv.lock

# If necessary
root <- renv::paths$root()
unlink(root, recursive = TRUE)
```

The scripts are also coded to install the necessary packages when not using renv, so not using renv will not create troubles.

## Reproducibility steps

> [!IMPORTANT]
> File size limitations on github prevented us to upload the .RData output of the `01_radiocarbon.R` script and the DEM to use in the `05_density_map.R`. You will need to download the DEM to reproduce the analyses and the RData file if you don't want to rerun the whole script. Inside the scripts a function will open a download page if the files are not found locally.

> [!WARNING]
> _These files need to be placed at the exact location specified below so that the scripts can find them_.

1. Clone or download the repository.
2. Download the DEM from [Figshare](https://figshare.com/ndownloader/files/58300900) and placed it in the directory `data/raw/raster` inside the repository folder.
3. (optional) Download the radiocarbon_data.RData file from [Figshare](https://figshare.com/ndownloader/files/58302040?private_link=61729c5d2acfb1f2d6a8) and place it in the directory `output/rda` inside the repository folder.
4. Run the `01_radiocarbon.R` script, which will generate the necessary data for the other scripts and also Figures 3, 4, and 5.
5. Run the `02_settlements.R` script, which will generate archaeological proxy data for the `04` script and also Figure 6.
6. Run the `03_paleoclimate.R`, which will generate the paleoclimate data for the `04` script and also Figure 7.
7. Run the `04_paleoclimate_correlation.R` script, which will generate Figure 8.
8. Run the `05_density_map.R` script, which will also generate Figure 9.

If you run into errors or troubles, please open an issue here on github or contact one of the maintainers.

