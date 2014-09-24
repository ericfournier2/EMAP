# INSTRUCTIONS
#
# 1. For each experiment, gather all of your raw intensity files into a single folder. You
#    can have multiple folders if you have multiple experiments.
# 2. Using a text editor, create target files for your epigenetic and/or transcriptomic experiments.
#    The created file should be in the same directory as your raw intensity data for the experiment.
#    This file should follow the format of the example file (SAM.target).
#    That is to say, it should be tab-separated and contain at least 3 columns:
#     - Filename, the name of text file the containing the raw scanned intensities,
#     - Cy3, the name of the experimental condition for the green channel,
#     - Cy5, the name of the experimental condition for the red channel.
#    Since this analysis pipeline only supports contrasts between two conditions, you
#    will need to create separate target files for each pair of conditions you may have.
#
# 3. Substitute the values of the following two variables with the name of the target files
#    you have just created. Make sure the names are correct and that the file extension is present.
#    If you cannot see the file extension, it might be ".txt".
#    If you do not have a transcriptomic experiment, simply write "".
epigenetic_Target <- "HighEpi.target"
transcriptomic_Target <- "HighTrans.target"

# 4. Substitute the name of the folders containing your target files and raw data. 
#    If your target file are directly in the analysis folder, you can leave this empty.
epigenetic_Folder <- "Raw/Denise"
transcriptomic_Folder <- "Raw/Denise"

# 5. Substitute the path to the directory containing the analysis pipeline in the following
#    expression. If the path contains slashes, substitute all backslashes with forwars slashes.
setwd("C:/Dev/Projects/Epigenetics/cow\EMAP")

# 6. Give names to your epigenetic and transcriptomic experiments.
epigenetic_Name <- "Denise/High/Epigenetic"
transcriptomic_Name <- "Denise/High/Transcriptomic"

# 7. Set the "reference" condition for your experiment. If you are analysing a combined 
#    transcriptomic and epigenetic analysis, you MUST specify a reference condition, and
#    it must exist in both. Otherwise, if the reference condition has no importante, set 
#    this to "" and the script will pick one for you.
reference_Condition <- "Zero"

# 8 Set a name for the results of the "Combined" analysis.
combined_Name <- " Denise/High/Combined"

# 9. Set thresholds for significance.
epi_foldchange_Threshold <- log2(1.5)
epi_pvalue_Threshold <- 0.05

trans_foldchange_Threshold <- log2(1.5)
trans_pvalue_Threshold <- 0.05


# 10. If you have been using the beta version  of the chip (with probe names GT_HQ...), change
#     this to "v1". Otherwise (if the probe names are EDMA_MET...), keep this set as "v2".
VERSION <- "v2"

# 11. Open R, select "Source R code" from the "File" menu, and load this analysis script.
#     The script will automatically analyze your data, and output the results in the
#     "Results/Experiment_Name" directory. This process can take several minutes, or even hours
#     on low-end machines. Once it is completed, you can consult the EmbryoGENE Epigenetics Analysis
#     Pipeline Manual for a detailed explanation of what each output file contains.

###################################################################################################
#                                       END OF INSTRUCTIONS                                       #
###################################################################################################

source("Scripts/main.R")