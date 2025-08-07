z# SALVE: prediction of interorgan communication with transcriptome latent space representation 

## Abstract
Massive transcriptomics data allow gene relationships to be discovered from their correlated expression. We describe SALVE, a method to infer the associations between secretome-encoding transcripts and gene modules in a distal organ from RNA sequencing data. This method builds upon similar bioinformatics approaches by introducing transcriptome latent space representations and transfer learning to simultaneously increase discovery power and predict downstream functional associations. Applied to GTEx v8 data, we show the method readily recapitulates canonical endocrine relationships, including insulin and adiponectin signaling, while inferring new candidate organokines and their signaling modality.  We also explore its utility for generating new hypotheses on cardiokine candidates and finding distal factors that may affect cardiac protein synthesis and metabolism. The predictions suggest a potential role of circulating galectin-3 (LGALS3) in regulating cardiac protein synthesis and homeostasis, which can be recapitulated in part in human induced pluripotent stem cell (hiPSC)-derived cardiomyocytes. This method may aid in ongoing efforts to delineate interorgan communications and endocrine networks in various areas of study.

## Data
- GTEx V8 data: https://www.gtexportal.org


## Code
### Processing GTEx v8 data
- 01_ReadGTExData.R: Read in GTEx v8, all samples, save as DDS, then batch correct with ComBAT
- 02_CompareGtex8BatchEffect.R: Assess batch effects of the GTEx v8 data prior and post ComBAT corrections
- 03_MakePLIERGTExVSTComBatAgeCorrNoSecretome.R: Applies the PLIER via wrapper function from MultiPLIER

### Analysis
- 10_CalculateTissueCorrelation_V2.R: Loop through all tissue pairs to generate correlation matrices.
- 14_CaclulateAllSScores.R: Calculate Ssec score between tissue pairs






