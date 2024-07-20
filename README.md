# COVID 19 Variant Classifer (98.9% accuracy)

## Intro
SARS-COV2 (COVID-19) impacted modern society in an unprecedented way, changing
government policies, altering societal norms, driving vaccine advancements, and fostering
drastic economic changes (Hosseinzadeh et al., 2022). One major factor driving COVID-19
virulence and persistence is their high viral mutation rate, giving rise to a multitude of variants,
including major strains like Alpha, Beta, Delta, and Omicron (Markov et al., 2023). Although
each strain differs in the range of mutations present, a major indicator differentiating between
strains is mutations to the surface glycoproteins (spike proteins) responsible for entry and
invasion of the host cell (Markov et al., 2023). This project attempts to create a high-accuracy
identifier of COVID-19 major variants based on spike protein sequence differences to aid the
classification and investigation of COVID-19 variants.

## Methods
We derive COVID-19 sequence data from the NCBI nucleotide database, isolate spike protein sequences, and explore different distance measures, parameters, and machine learning algorithms and their impact on the accuracy of the classifier.

## Results
Overall, the best performing model reached an accuracy of 98.9% running on a Random Forest framework with k-mer 3. 

## Figures
![image](https://github.com/user-attachments/assets/13f736e9-c464-42fc-9963-a9680be2a622)

## Citation
Please cite this work if used for research/industrial purposes.
