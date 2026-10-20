The 2026 North Carolina blue crab (Calinectes sapidus) benchmark stock assessment requires the following:

1) The Assessment model and figure code, which contains the ADMB .tpl file, compiled executable (.exe),
   other supporting .R files used to connect ADMB model output with R, plots from the original size structured
   model (Cao et al. 2017), and biological reference point code
2) The main R files for running the model: 01_RunModel.R, 01_SetData.R, 01_DataGen.R, and 01_Functions_old.R, which set up the model
   configuration and run the model from R
3) The data files: BC_catch_2seas.csv, which contains the seasonal catch and fishery-dependent size composition data
   during 1982-2024; and BC_survey_2seas.csv, which contains the fishery-independent index of relative abundance
   and combined size composition data

How to run the stock assessment:

1) Create a folder on your desktop, e.g., "2026 NC blue crab stock assessment model"
2) Download the contents of this repository to that folder
3) Use R Studio to open the R Project file: "size_structured_assessment_model-master"
4) Within that R Project, open the R file: "01_RunModel.R"
5) Go to line 64 of 01_RunModel.R and give the scenario variable a unique name for model results and plotting code to go to
6) Run the whole 01_RunModel.R script

Additional details: 

You can change which model parameters are estimated directly in the 01_RunModel.R file, but changes to other fishery and life history 
characteristics have to be changed manually in the 01_SetData.R file. 

The paper documenting the model by Cao et al. (2017) and associated user guide are located within the Size structured assessment model paper and guide folder.

To access the original model code, visit Dr. Jie Cao's GitHub: https://github.com/SOUKETU > size_structured_assessment_model

If you encounter any issues running the code, or have questions about the model and its outputs, 
contact Dr. Matt Damiano at matt.damiano@deq.nc.gov







