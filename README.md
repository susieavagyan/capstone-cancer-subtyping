# capstone-cancer-subtyping
# Molecular Subtyping of Colorectal Cancer: In the Frontiers of Personalized Diagnostics and Treatment

### Abstract
Being the third most diagnosed and second most deadly cancer worldwide, colorectal cancer is a
highly complex, multigenic disease that has very high inter-patient variability in terms of the genetics of the
tumor. This raises the need for developing a personalized treatment for CRC patients for better efficacy and
reduced toxicity. Molecular subtyping of the disease is a way to define biological subgroups for which targeted
treatment can be optimized. Our research aimed to test different Machine Learning and Deep Learning models
that combined theoretically and practically tested state-of-the-art concepts to obtain biologically meaningful clusters from somatic mutations and copy number alterations as CRC patient subgroups. Four different methods with different types of inputs were tested, Spectrum, xGeneModel,Kmeans clustering, and Deep Embedded Clustering. The obtained clusters were treated as labels to build classifiers as a predictive tool for incoming patient records. Survival Analysis was also performed to analyse obtained clusters.

This repository contains all the scripts used for the project (R and Jupyter Notebook files), as well as the datasets used, obtained figures and results, and the package and library versions which have been in use at the time of writing and running the scripts.

### Data, Preprocessing, EDA
Obtaining, preprocessing and initial analysis of the data is done using the script **data_prep.R**.
The data was dumped directly from cBioportal's public datahub repository in Github. The link to that repository is **https://github.com/cBioPortal/datahub/tree/master/public/msk_met_2021**, the link to the original cBioportal study page is **https://www.cbioportal.org/study/summary?id=msk_met_2021**, also referenced in the paper. The raw and preprocessed mutational data files are saved in the **data** folder (*mut_raw.csv, cnv_raw.csv, mut_cnv_onehot.csv*). We also obtain clinical data with this script, which is saved in the **data** folder as well (*clinical_patients.csv, clinical_primary_raw.csv*). The plots for EDA are saved in the **figures** folder.

### Data Enrichment
Enrichment/stratification step was performed using the script **som_analysis.R**. This code runs the full oposSOM pipeline on the data and outputs large files as results, which are not saved in this directory (the pipeline takes around 15-20 minutes to run). The information needed is the propagated data matrix, which is obtained explicitly in the script and saved in the **results\som** folder (*propagated_mut.csv*)

### Clustering Models
For each clustering algorithm a separate file is in the directory. 

1. Spectrum model is ran by the script **clustering_models_Spectrum.R** passing the non propagated data **mut_cnv_onehot.csv** as input, and the results are saved in **results\spectrum** (cluster assignments  *spec_cluster_assignments_5.csv* and cluster statistics *cluster_stats_plot.pdf*, *cluster_stats.csv* )
2. xGene model is ran by the script **clustering_models_xGene.R** passing the non propagated data **mut_cnv_onehot.csv**, functional similarity matrix **FSM.txt**, and weighted similarity matrix **WM.txt** as input, and the results are saved in **results\xgene** (cluster assignments  *cluster_labels.csv*, intermediary distnace matrix *distance.txt*, clustering evaluation results *clustering_and_evaluation.pdf*, *result_overview.txt* )
3. KMeans model is ran by the script **kmeans_stratified.ipynb** passing the propagated data **propagated_mut.csv**, functional similarity matrix **FSM.txt**, and weighted similarity matrix **WM.txt** as input, and the results are saved in **results\kmeans** (cluster assignments  *kmeans_stratified_cluster_assignments.csv*) 
4. DEC model is ran by the script **DEC.ipynb** passing the propagated data **propagated_mut.csv** as input, and the results are saved in **results\DEC** (cluster assignments  *dec_cluster_assignments_5.csv*). We also obtained PCA plots of this clustering, which are saved in the **figures** folder (*5_clusters_dec_pca2.png*, 5_clusters_dec_pca3.png*)

### Classification Models

Classification models are in the script **classification.ipynb**. It takes 2 inputs: cluster assignments (of a specific clustering method) and data(propagated or non-propagated depending on from which data were the cluster assignments obtained from). 10 classifiers (9 traditional ML classifiers from sklearn library and 1 custom DL classifier - MLP) were constructed in this script. The best parameters for each classifier were obtained using GridSearchCV function. The performance evaluation metrics were measured and barcharts of those performance measures for all trained classifiers for each clustering were saved in the **figures** folders. The best performing classifier information was also obtained for each clustering and saved in the respective subfolders of the **results** folder (*best_*_pm.csv*)

### Biological Cluster Analysis

Further biologic analysis of clusters was performed via the script **Cluster_Analysis.ipynb**. This script takes as input the cluster assignments file (of a specific clustering method) and the clinical data of samples (*clinical_primary_raw*) and patients (*clinical_patients.csv*). This script performs Batch Effect Counfounding Analysis, and the obtained figures of counfounder distributions per cluster are saved in the **figures** folder (*batch_\*.png*). Kaplan Meier Survival Analysis is also contained in this script. The resulting KM survival curves are saved in the  **figures** folder (*KM_curves\*.png*).

### *Used packages and libraries

The list of R packages and Python libraries with respective versions used in this project are saved in the files **session_info_R.txt** for R (using the script *session_info.R*) and **requirements_PY.txt** for Python (using pip freeze command from terminal)






