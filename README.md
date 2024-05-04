# Computational-neuroanatomy-v1
The code developed for building the computational neuroanatomy atlas and biclusters as described in Kaczanowska et al, "Molecular archaeology of human cognitive traits"

Copyright (C) 2022 VRVis.
All rights reserved.
Contact: VRVis Forschungs-GmbH (office@vrvis.at)

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. All advertising materials mentioning features or use of this software
   must display the following acknowledgement:
   This product includes software developed by the VRVis Forschungs-GmbH.
4. Neither the name of the VRVis Forschungs-GmbH nor the
   names of its contributors may be used to endorse or promote products
   derived from this software without specific prior written permission.

Author Florian Ganglberger (ganglberger@vrvis.at)

# General information
1) Before you can start the analysis, you need to download gene expression data from
the Allen Human Brain Atlas. You can do this with the MATLAB
scripts run_download_AHBA.m provided in MATLAB_DATA_DOWNLOAD. Before you can run this file,
you need to download  parse_json.m from https://www.mathworks.com/matlabcentral/fileexchange/20565-json-parser 
and copy it to this folder MATLAB_DATA_DOWNLOAD. Then you just need to run run_download_AMBA_AHBA.m , which can
take several hours or up to days.
We also provide the output of this analysis at storage/all_genes_expression_Human.mat, storage/atlasRegionsHuman.mat and storage/ontologyHuman.json

2) For the analysis, also taskFMRI network data from the Human Connectome Project is needed. The task data is available
at the 900 Subjects Data Release  https://www.humanconnectome.org/study/hcp-young-adult/document/900-subjects-data-release).
The task fMRI data HCP_S900_787_tfMRI_ALLTASKS_level3_zstat1_hp200_s2_MSMAll.dscalar.nii needs to be copied to the storage folder.
It can be converted to biopsy-site level (raw data is on cortical surface level) with generateTaskFMRIDataOnBiopsyLevel.m 
in MATLAB_DATA_DOWNLOAD. Note that this script  requires the human connectome project workbench 
from http://www.humanconnectome.org/software/connectome-workbench.html
We also provide the output of this analysis at storage/taskfMRI.mat

3) Analysis can be performed with the files in R_MAIN_CODE. It contains a R notebook (COMPUTATIONAL_NEUROANATOMY.Rmd) that
generates the figures from the paper step-by-step provided with comments (with R version 4.0.3). The analysis requires the prediction of functional
neuroanatomical maps with updated code from a previous publication Ganglberger et al [1]. It can be found also in this folder
(predict_human_expression_synergy.R and predict_mouse_expression_synergy) and will be automatically called in the R notebook.
After performing the analysis, subspace patternmining with biclustering (GABi) can be performed with GABI_BICLUSTERING.R. 
Please note that for this pattern mining can take several days to weeks without an HPC environment.

[1] Ganglberger, F., Kaczanowska, J., Penninger, J. M., Hess, A., Bühler, K., & Haubensak, W. (2017). 
Predicting functional neuroanatomical maps from fusing brain networks with genetic information. 
NeuroImage. https://doi.org/10.1016/j.neuroimage.2017.08.070


# File overview

1. MATLAB_DATA_DOWNLOAD
- run_download_AHBA.m 
-- script to download gene-expression and connectivity data
- get_gene_expression_human.m
- generateTaskFMRIDataOnBiopsyLevel.m
-- generates task FMRI data on biopsy level for further analysis
	
2. R_MAIN_CODE	
- COMPUTATIONAL_NEUROANATOMY.Rmd
-- central R notebook file to perform the analysis and generate paper figures
	COMPUTATIONAL_NEUROANATOMY.html
		html file which shows R notebook with all figures/tables
	GABI_BICLUSTERING.R
		file for subspace pattern mining (biclustering) after the analysis in COMPUTATIONAL_NEUROANATOMY.Rmd
	GABi_fixed.R
		original GABi code from (https://cran.r-project.org/web/packages/GABi) with bugfixes
	GABi_parallel.R
		parallelized and optimized version of GABi_fixed.R
	predict_human_expression_synergy.R
		code for prediction of functional neuroanatomical maps, first order (gene expression synergy only) for human
	
3. storage
		folder that will contain intermediate files, downloaded data and results for the analysis
	dnds_data_table
		contains DNDS data table used for the analysis
	taskfMRI.mat
		output of generateTaskFMRIDataOnBiopsyLevel.mat
	all_genes_expression_Human.mat
		brainwide gene expression at 3702 biopsy sites--> output of get_gene_expression_human
	atlasRegionsHuman.mat
		region ID association of 3702 biopsy sites --> output of get_gene_expression_human
	ontologyHuman.json
		region ID association to brain regions -> ooutput of get_gene_expression_human
	results_functionalmaps_human
		results generated by predict_human_expression_synergy.R
	biclustering_results
		results generated by the biclustering

		
