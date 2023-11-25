# Proteomics analysis

Pathway activities in 7 tissues are analyzed using proteomics data obtained from Nie, Xiu, et al. "Multi-organ proteomic landscape of COVID-19 autopsies." Cell 184.3 (2021): 775-791.
Author: Nhung Pham April 2022
The work is published in Pham, Nhung, et al. "Tissue-specific pathway activities: A retrospective analysis in COVID-19 patients." Frontiers in immunology 13 (2022): 963357. 
# Code structure:

    Data folder contains proteomics data, wikipathways and COVID19 Disease Map GMT files that are needed for the analysis
    "MultiOrganProteomicFunction.R" contains custom-made functions that needed for the analysis
    "WilcoxonTissueAnalysis.R" is the master script for the analysis

# How to run:

    Download the code and the Data folder
    Follow the steps in the "WilcoxonTissueAnalysis.R" to recreate the result in the paper
