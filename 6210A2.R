
#### LOAD LIBRARIES & SEARCH DATABASE----
# install.packages("Tidyverse", "ggpmisc") #install packages if needed

#load libraries
library("tidyverse")
library("rentrez")
library("Biostrings")
library("stringr")

#check databases & search terms available
entrez_dbs()
entrez_db_searchable(db = "nuccore")

#search for entries of each variant and save the search history
#search terms: ensure it is SARS-COV2, contains gene for S protein, RNA-based only (no DNA entries) for consistancy, and variant as indicated below
#VARIENTS: Alpha = B.1.1.7 (1000), Beta = B.1.351 (100), Delta = B.1.617.2 (500), Omicron = BA.1 (only original omicron lineage taken for simplicity) (according to https://cov-lineages.org/ PANGOLIN COVID lineage)
Alpha <- entrez_search(db = "nuccore", term = "Severe acute respiratory syndrome coronavirus 2[ORGN] AND surface glycoprotein[PROT] NOT DNA[ALL] AND B.1.1.7[ALL]", retmax = 1, use_history=TRUE)
Beta <- entrez_search(db = "nuccore", term = "Severe acute respiratory syndrome coronavirus 2[ORGN] AND surface glycoprotein[PROT] NOT DNA[ALL] AND B.1.351[ALL]", retmax = 1, use_history=TRUE)
Delta <- entrez_search(db = "nuccore", term = "Severe acute respiratory syndrome coronavirus 2[ORGN] AND surface glycoprotein[PROT] NOT DNA[ALL] AND B.1.617.2[ALL]", retmax = 1, use_history=TRUE)
Omicron <- entrez_search(db = "nuccore", term = "Severe acute respiratory syndrome coronavirus 2[ORGN] AND surface glycoprotein[PROT] NOT DNA[ALL] AND BA.1[ALL]", retmax = 1, use_history=TRUE)

#put variant names and search results into a list (iteratable later)
variant_results <- list(Alpha, Beta, Delta, Omicron)
variant_names <- list("Alpha", "Beta", "Delta", "Omicron")

#NOTE: number of entries varies as seen below. Max of 1000 entries taken from each group. (lack of entries may cause bias in machine learning model and cause more uncertain labeling)
#check total number of entries found for SARS-COV2 Alpha, Beta, Delta, Omicron respectively:
print(c(Alpha$count, Beta$count, Delta$count, Omicron$count))

#NOTE: there is clear imbalance in the number of entries found. For the purposes of this classifier there are too few beta variant entries for complete undersampling and no extra data available from the database in our scope, thus we will accept the potential class imbalance bias and only undersample to 1500 entries max. One option we an try is to oversample by removing a random character from each Beta variant entry to artificially generate more data, however I ran out of time for this project. 

#### DOWNLOAD & CLEAN DATA----


#NOTE: many steps below are piped/looped to save memory/variable use, but can (and have been) be easily taken apart to test that they work as intended.


#initialize variables/counters
dfCOVFull=data.frame() #set up overall dataframe
i <- 1 #counter

#loop to fetch data for each variant (this section in the novel code component for data extraction)
for (variant_name in variant_names){
  
  #fetch (using browser history because too many files) all CDS from a set number of entries 
  entrez_fetch(db = "nuccore", web_history=variant_results[[i]]$web_history, rettype = "fasta_cds_na", retmax=1500, use_history = TRUE) %>% 
    #isolate only sequences (CDS) coding surface glycoprotein (spike protein) using regex
    str_extract_all(">[^\\]]+\\]\\s\\[protein\\=surface\\sglycoprotein[^>]+", simplify = TRUE) %>% 
    #write to temporary file
    write("COV_S_proteins_temp.fasta", sep = "\n")
  
  #read as DNAStringSet
  stringSet <- readDNAStringSet("COV_S_proteins_temp.fasta")
  
  #convert our data to a dataframe with the correct variant type, remove empty sequences (due to regex searches) and append new set of data to the current set
  dfCOVFull <- data.frame(Variant=variant_name, id=names(stringSet), S_Sequence=paste(stringSet)) %>% 
    #NOTE: above you can check id column to double check that regex extraction worked (should all be gene=S and protein=surface glycoprotein) 
    filter(str_length(S_Sequence)>"100") %>% 
    #NOTE: you can rehome the above step to outside the for loop to check number of sequences removed (but it is placed here to save memory for my computer)
    rbind(dfCOVFull) #append to dataframe
  
  i<- i+1 #increment counter
  cat("Processing data set",i-1,"...") #note for user
}

cat("Dataframe creation DONE") #note for user

#CHECK DATA

#some tests to ensure code working as expected
# View(dfCOVFull) #Optional view dataframe to check
# dfCOVFull[1,2] == dfCOVFull[11,2] #Check that values are not being replicated (on simplified model with retmax=10)


dfCOVFull %>% group_by(Variant) %>% count() #check the number of available data after simple processing
sum(is.na(dfCOVFull$S_Sequence)) #should be no NAs in data (All false (0), therefore sum = 0)

#check length of sequences to ensure they are more or less the same protein...
str_length(dfCOVFull$S_Sequence) %>% hist(xlab = "Sequence Length (NAs)", ylab = "Frequency", main = "Frequency Histogram of Spike Protein Sequence Lengths", ) 
str_length(dfCOVFull$S_Sequence) %>% table()
str_length(dfCOVFull$S_Sequence) %>% summary()
#NOTE: very good distribution, except clear outlier at 2562. Remove these datapoints.
dfCOVFull <- dfCOVFull %>% filter(str_count(S_Sequence) > 3000)

#check again, looks good!
str_length(dfCOVFull$S_Sequence) %>% hist(xlab = "Sequence Length (NAs)", ylab = "Frequency", main = "Frequency Histogram of Spike Protein Sequence Lengths", ) 
str_length(dfCOVFull$S_Sequence) %>% table()
str_length(dfCOVFull$S_Sequence) %>% summary()



#### RANDOM FOREST DATA PROCESSING & DISTANCE MEASURES----

#process data for RandomForest (remove beginning and trailing Ns, gap indicators, and limit number of Ns to <10%)
dfCOVFull <- dfCOVFull %>%
  mutate(S_Sequence2 = str_remove(S_Sequence, "^[-N]+")) %>%
  mutate(S_Sequence2 = str_remove(S_Sequence2, "[-N]+$")) %>%
  mutate(S_Sequence2 = str_remove_all(S_Sequence2, "-+")) %>%
  filter(str_count(S_Sequence2, "N") < (0.1 * str_count(S_Sequence))) #keep sequences with <10% Ns

#check amount of remaining data (still imbalanced as expected)
dfCOVFull %>% group_by(Variant) %>% count()


#CALCULATE DISTANCE MEASURES

#convert to DNA stringset to work with
dfCOVFull$S_Sequence2 <- DNAStringSet(dfCOVFull$S_Sequence2)

#1) ATG PROP

#find letter frequency and put into new columns
dfCOVFull <- cbind(dfCOVFull, as.data.frame(letterFrequency(dfCOVFull$S_Sequence2, letters = c("A", "C","G", "T")))) 
#calculate proportions of each 
dfCOVFull$Aprop <- (dfCOVFull$A) / (dfCOVFull$A + dfCOVFull$T + dfCOVFull$C + dfCOVFull$G)
dfCOVFull$Tprop <- (dfCOVFull$T) / (dfCOVFull$A + dfCOVFull$T + dfCOVFull$C + dfCOVFull$G)
dfCOVFull$Gprop <- (dfCOVFull$G) / (dfCOVFull$A + dfCOVFull$T + dfCOVFull$C + dfCOVFull$G)

dfCOVFull <- dfCOVFull[1:11]

#2) KMER OF 2, 3, and 4

dfCOVFull <- cbind(dfCOVFull, as.data.frame(dinucleotideFrequency(dfCOVFull$S_Sequence2, as.prob = TRUE)))
dfCOVFull <- cbind(dfCOVFull, as.data.frame(trinucleotideFrequency(dfCOVFull$S_Sequence2, as.prob = TRUE)))
dfCOVFull <- cbind(dfCOVFull, as.data.frame(oligonucleotideFrequency(x = dfCOVFull$S_Sequence2, width = 4, as.prob = TRUE)))

#convert back into char
dfCOVFull$S_Sequence2 <- as.character(dfCOVFull$S_Sequence2)


#### TRAINING/VALIDATION DATA SEPERATION----

#ensure repeatable results
set.seed(217)

#20% of data in each group/variant should be validation data
dfValidation <- dfCOVFull %>%
  group_by(Variant) %>%
  sample_frac(0.2)

#the other 80% is used for training
dfTraining <- dfCOVFull %>%
  filter(!id %in% dfValidation$id) %>%
  group_by(Variant) %>%
  sample_frac(0.8)

#check amount of training and validation data
dfValidation %>% group_by(Variant) %>% count()
dfTraining %>% group_by(Variant) %>% count()

#set accuracy function for validation later 
accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100} #overall accuracy
accuracyPerVariant <- function(x){diag(x)/(rowSums(x)) * 100} #accuracy by variant
accuracyWeighted <- function(x){sum(diag(x)/(rowSums(x)))/4 * 100} #accuracy weighted equally for each variant (average of accuracy by variant)


#### RANDOM FOREST 1 (ATG) ----

set.seed(200) #reproducible results

#BUILD CLASSIFIER

gene_classifier <- randomForest::randomForest(x = dfTraining[, 9:11], y = as.factor(dfTraining$Variant), ntree = 50, importance = TRUE)

#EVALUATE

#training
gene_classifier #check error rate estimate
#NOTE: from confusion matrix, we can see there is higher error for the undersampled groups as expected (beta and delta)

#validation
predictValidation <- predict(gene_classifier, dfValidation[, c(1, 9:11)])
tab <- table(observed = dfValidation$Variant, predicted = predictValidation) #validation confusion matrix
tab

#check actual prediction accuracy
accuracy(tab)
accuracyPerVariant(tab)
accuracyWeighted(tab)
#90.2% accuracy, pretty good for first try! but weighted accuracy is quite low due to beta group
#again, as we can see, there is higher error for the undersampled groups (beta and delta)



#### RANDOM FOREST 2 (KMER2) ----

#BUILD CLASSIFIER

#note using default ntree=50 as a simple assessment first
gene_classifier <- randomForest::randomForest(x = dfTraining[, 12:27], y = as.factor(dfTraining$Variant), ntree = 50, importance = TRUE)

#EVALUATE

#training
gene_classifier #check error rate estimate
#NOTE: from confusion matrix, we can see there is higher error for the undersampled groups as expected (beta and delta)

#validation
predictValidation <- predict(gene_classifier, dfValidation[, c(1, 12:27)])
tab <- table(observed = dfValidation$Variant, predicted = predictValidation) #validation confusion matrix
tab

#check actual prediction accuracy
accuracy(tab)
accuracyPerVariant(tab)
accuracyWeighted(tab)
#much higher accuracy!! 97.7% overall and 92.3% weighted. Beta still lagging behing at only 78% accuracy.


#### RANDOM FOREST 3 (KMER3) ----

#BUILD CLASSIFIER

gene_classifier <- randomForest::randomForest(x = dfTraining[, 28:91], y = as.factor(dfTraining$Variant), ntree = 50, importance = TRUE)

#EVALUATE

#training
gene_classifier #check error rate estimate
#NOTE: from confusion matrix, we can see there is higher error for the undersampled groups as expected (beta and delta)

#validation
predictValidation <- predict(gene_classifier, dfValidation[, c(1, 28:91)])
tab <- table(observed = dfValidation$Variant, predicted = predictValidation) #validation confusion matrix
tab

#check actual prediction accuracy
accuracy(tab)
accuracyPerVariant(tab)
accuracyWeighted(tab)
#wow, 100% accuracy for alpha and beta variants! Beta and delta still lagging behind, causing overall accuracy of 98.9 and weighted accuracy to be 95.3. Beta at 85.7%


#### RANDOM FOREST 4 (KMER4) ----

#BUILD CLASSIFIER

gene_classifier <- randomForest::randomForest(x = dfTraining[, 92:347], y = as.factor(dfTraining$Variant), ntree = 50, importance = TRUE)

#EVALUATE

#training
gene_classifier #check error rate estimate
#NOTE: from confusion matrix, we can see there is higher error for the undersampled groups as expected (beta and delta)

#validation
predictValidation <- predict(gene_classifier, dfValidation[, c(1, 92:347)])
tab <- table(observed = dfValidation$Variant, predicted = predictValidation) #validation confusion matrix
tab

#check actual prediction accuracy
accuracy(tab)
accuracyPerVariant(tab)
accuracyWeighted(tab)
#not much difference here, we have reached the limit to how well we can differentiate these variants based on poor class balance for random forest (we will no longe rbe looking at kmer 4)
#interestingly, same confusion matrix as kmer3. Perhaps those few inaccurate sequences have differences that make them difficult to group? Could they represent transition states between alpha-beta (etc.) variants? Very interesting questions for the future and also great way to examine the evolution of COVID spike proteins


#### K-NN () ----

#import library
library(class)

#normalization function
nor <-function(x) { (x -min(x))/(max(x)-min(x))   }

set.seed(200)

#BUILD CLASSIFIER (ATG)

#normalize training and validation data
dfValidationNorm <- as.data.frame(lapply(dfValidation[,c(9:11)], nor))
dfTrainingNorm <- as.data.frame(lapply(dfTraining[,c(9:11)], nor))

#create prediction classifier, use k=9 (various k were tried: recommended sqrt of all data did not work, but sqrt of smallest stat group works well! (9*9=81 == approx beta data size))
pr <- knn(train = dfTrainingNorm, test = dfValidationNorm, cl = dfTraining$Variant, k=9)
  
#confustion matrix
tab <- table(dfValidation$Variant, pr)
tab

#check actual prediction accuracy
accuracy(tab)
accuracyPerVariant(tab)
accuracyWeighted(tab)
#interestingly more poorly predicting for delta than beta...


#BUILD CLASSIFIER (KMER2)

#normalize training and validation data
dfValidationNorm <- as.data.frame(lapply(dfValidation[,c(12:27)], nor))
dfTrainingNorm <- as.data.frame(lapply(dfTraining[,c(12:27)], nor))
#create prediction classifier, use k=9 (various k were tried: recommended sqrt of all data did not work, but sqrt of smallest stat group works well! (9*9=81 == approx beta data size))
pr <- knn(train = dfTrainingNorm, test = dfValidationNorm, cl = dfTraining$Variant, k=9)
#confusion matrix
tab <- table(dfValidation$Variant, pr)
tab
#check actual prediction accuracy
accuracy(tab)
accuracyPerVariant(tab)
accuracyWeighted(tab)
# overall 97.4, 92.2 weighted, beta performed most poorly at 78.6

#BUILD CLASSIFIER (KMER3)

#normalize training and validation data
dfValidationNorm <- as.data.frame(lapply(dfValidation[,c(28:91)], nor))
dfTrainingNorm <- as.data.frame(lapply(dfTraining[,c(28:91)], nor))
#create prediction classifier, use k=9 (various k were tried: recommended sqrt of all data did not work, but sqrt of smallest stat group works well! (9*9=81 == approx beta data size))
pr <- knn(train = dfTrainingNorm, test = dfValidationNorm, cl = dfTraining$Variant, k=9)
#confustion matrix
tab <- table(dfValidation$Variant, pr)
tab
#check actual prediction accuracy
accuracy(tab)
accuracyPerVariant(tab)
accuracyWeighted(tab)
#higher accuracy: 98.8 overall, 96.9 weighted, beta at 92.6 (Better than Rand forest for class imbalance!)




#### VISUALIZATION 1: Representative tree ----

#installation (from https://stats.stackexchange.com/questions/41443/how-to-actually-plot-a-sample-tree-from-randomforestgettree)
options(repos='http://cran.rstudio.org')
have.packages <- installed.packages()
cran.packages <- c('devtools','plotrix','randomForest','tree')
to.install <- setdiff(cran.packages, have.packages[,1])
if(length(to.install)>0) install.packages(to.install)

library(devtools)
if(!('reprtree' %in% installed.packages())){
  install_github('munoztd0/reprtree')
}
for(p in c(cran.packages, 'reprtree')) eval(substitute(library(pkg), list(pkg=p)))

#load libraries
library(randomForest)
library(reprtree)

#draw tree for random forest ATG (for simplicity)
set.seed(500)
dfTrainingModel <- round(dfTraining[, 9:11], 3) #round smaller so we can see in tree visual.
dfTrainingModel$Variant <- as.factor(substr(dfTraining$Variant, 1, 1)) #take only first letter of variants so we can see in visual
model <- randomForest(Variant ~ ., data=dfTrainingModel,  ntree = 100, importance = TRUE, do.trace = 100) #train model
model
reprtree::plot.getTree(model, main="Representative Tree for AGT Proportion") #draw tree

#repeat for KMER2
set.seed(500)
dfTrainingModel <- round(dfTraining[, 12:27], 3)
dfTrainingModel$Variant <- as.factor(substr(dfTraining$Variant, 1, 1))
model <- randomForest(Variant ~ ., data=dfTrainingModel,  ntree = 100, importance = TRUE, do.trace = 100)
model
reprtree::plot.getTree(model, main="Representative Tree for K-mer 2")

#repeat for KMER 3
set.seed(500)
dfTrainingModel <- round(dfTraining[, 28:91], 3)
dfTrainingModel$Variant <- as.factor(substr(dfTraining$Variant, 1, 1))
model <- randomForest(Variant ~ ., data=dfTrainingModel,  ntree = 100, importance = TRUE, do.trace = 100)
model
reprtree::plot.getTree(model, main="Representative Tree for K-mer 3")

