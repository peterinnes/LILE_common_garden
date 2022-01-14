## This repository contains all the R scripts to reproduce results from Innes et al 2022, "Assessment of biogeographic variation in traits of Lewis flax (Linum lewisii) for use in restoration and agriculture"

First, run the code within 'Ephraim_EDA.R.' This contains an exploratory analysis of the Ephraim garden trait data, and also applies filters to the data (some Lewis flax accessions are excluded here because they were found to be contaminated with plants from the 'Appar' accession). This script also contains code for the calculation of the composite traits estimated seed yield and fecundity. Importantly, the .csv files generated within thisscript are used for downstream data analysis.  

Next, analysis of Ephraim non-oil traits, Ephraim oil traits, and Millville traits is carried out separately in three different scrips: 'Ephraim_Analysis.R' 'Oil_Analysis.R' and 'Millville_Analysis.R'. Generation of some of the combined figures (e.g. Figure 4) relies on results from all three of these scripts.

'misc_Env_Analysis.R' contains code for the transfer distance test and the environmental PCA.

'Maps.R' contains code for creation of Figure 1.

All data is stored in the accompanying Dryad repository: https://doi.org/10.5061/dryad.tdz08kq1n