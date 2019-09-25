# ‘Lentimonas’ use an array of exo-enzymes to degrade sulfated fucans in the ocean
## Sichert et al., 2019

This repository contains the analyses for Figure 6 ***Verrucomicrobia* are abundant and specialized polysaccharide degraders**

To begin this analysis follow the instructions below:

1. Clone this repository
```
# Clone
git clone https://github.com/mschecht/lentimonas_env_analysis.git

# Move into repo
cd lentimonas_env_analysis
```

2. Download the data
```
# Download
wget https://ndownloader.figshare.com/articles/9904793/versions/1 -O mitags.zip

# Unzip
unzip mitags.zip

# Move data to correct directories
mv tara_mitag_tax.tsv.gz data/raw/tara/ && mv osd2014_mitag_tax.tsv.gz data/raw/osd/
```

3. Open the `plotting_Verrmucomicrobia.Rmd` file in the `notebooks/` directory in RStudio and being the analyses
