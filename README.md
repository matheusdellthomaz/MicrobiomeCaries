Here you can find the R scripts, RDS files, csv and txt datasets used in Oral bacteria makers of caries: and in silico machine learning analysis

MicrobiomeCombinations.R is the script used for testing different ML algorithms, ASVs combinations and combination sizes developed for predicting Carie STATUS from 5 ASVs. In this script the algorithms available are xgboost, glm, rf, and svm.

Microbiome.R contains the feature selection and importance analysis, datasets management and metrics, ML, correlation, GLM, ANOVA, PCA and PLS-DA analysis and plotting.

The RDS files named "RF765_393_1445_273_100X.rds" and "GLMNet669_2294_273_765_652.rds" are, respectively, the Random Forest and GLMNet models developed in this study, and the "EERF765_393_1445_273_100.rds" and "EEGLMNet669_2294_273_765_652.rds" are their respectively explainer files. 
The RDS files can be loaded directly to get predictions from other datasets, run the script and get the STATUS predictions based on the ASVs used.

The CSV file named "subset_db_rf.csv" is the centered log-ratio transformed compositional dataset used for training the ML algorithms and subsequent analysis.

The txt file named "Taxonomy.txt" contains each ASV taxonomy.
The txt file named "Metadata.txt" contains participant's metadata.
The txt file named "Metadata.txt" contains participant's metadata.
The txt file named "ASV Table.txt" contains the untransformed ASV table for all studies participants.
