# rf_hvtn505
This document aims to instruct how to reproduce results in the manuscript "Improving Random Forest Predictions in Small Datasets from Two-phase Sampling Designs". The source codes can be generally split into two main R files "myhelper.R" and "experiments". The detailed explanations are as follows.

## myhelper.R
This R file contains several helper functions that help to easly get results.

### Function list
- *screen_lasso*: implements lasso variable screening.
- *screen.dat.index*: generates a dataset and a screened variable index. Users can specify a candidate marker set and whether lasso variable screening is applied.
- *get_nms_group_all_antigens*: generates variable index by an assay type.
- *get.rf.cvauc*: fits different types of random forests(standard random forest (RF), RF with under-sampling (RF_under), RF with over-sampling (RF_over), and tuned RF (tRF)). Users can specify inver sampling probability weighting (ipw).
- *get.glm.cvauc*: fits generalized linear models. Users can specify inver sampling probability weighting (ipw).
- *get.st.cvauc*: fits stacking that is based on random forests and generalized linear models. 

## experiments.R
This R file containes codes for all experiments in the manuscript.

Step 1) Importing HVTN 505 dataset
- i) First download the HVTN 505 dataset at https://atlas.scharp.org/cpas/project/HVTN%20Public%20Data/HVTN%20505/begin.view
- ii) Install R package HVTN505 using the code "devtools::install_local("HVTN505_2019-4-25.tar.gz")".

Step 2) Conducting Experiments
- Codes in "2-1. Random forest (RF)" generate all random forest-based results in Table 1, 2, and 3 in the manuscript. You can customize whether clinical covariates are included or whether lasso variable screening is applied.
- Codes in "2-2. Generalized linear models (GLM)" generate all generalized linear models-based results in Table 3 and Table A.1 (supplemantary material). You can customize whether clinical covariates are included or whether lasso variable screening is applied.
- Codes in "2-3. Stacking (RF + GLM)" generate all stacking results in Table 4. To fit stacking, you should first import four R files, which are "helper_functions.R", "caretList.R", "caretEnsemble.R", and "caretStack.R", in "caretEmsemblerR" folder. These four files are copied from the *caretEnsemble* R package (Deane-Mayer and Knowles, 2016), and changes were made so that candidate learners in a stacking can be fit on different training datasets.
