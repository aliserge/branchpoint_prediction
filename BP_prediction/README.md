# Branchpoints prediction

# About project

## Goal: 	
Create ML-architecturre model to locate branchpoints in all introns in human genome. 

## Initial data: 
* [high confidence branchpoints](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4315302/bin/supp_gr.182899.114_Supplemental_TableS1.xlsx) from experimental work [1] for GRCh37/hg19 assembly
* genome assembly [GRCh37/hg19](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/)
* genome annotation [GENCODE v19](https://genome.ucsc.edu/cgi-bin/hgTables) with introns coordinates 

# Dependences
Currently, this work is devoted to the creation of the architecture of the model, so it was performed in Jupiter notebooks. Accordingly, the following dependencies are required to run notebooks:
* [python 3.9.7](https://www.python.org/download/releases/3.0/)
* [jupyter notebook](https://jupyter.org/install)
* standard python libraries: [Biopython 1.79](https://biopython.org/), [numpy 1.20.3](https://numpy.org/doc/stable/index.html), [pandas 1.3.4](https://pandas.pydata.org/), [matplotlib 3.5.2](https://matplotlib.org/1.2.1/examples/pylab_examples/errorbar_demo.html), [seaborn 0.11.2](https://seaborn.pydata.org/)
* python libraries for Machine learning: [scikit-learn 0.24.2](https://scikit-learn.org/stable/), [lightgbm 3.3.2](https://lightgbm.readthedocs.io/en/latest/), [catboost 1.0.5](https://catboost.ai/), [xgboost 1.6.0](https://xgboost.readthedocs.io/en/stable/index.html)

All python libraries can be installed using command `pip install <library name>`

Note: all jupyter notebooks are executable and can be runned using [papermill](https://papermill.readthedocs.io/en/latest/index.html): 

 ```
 papermill input.ipynb output.ipynb --autosave-cell-every 60 --no-progress-bar --no-log-output & 
 disown #& and disown for background
 ```
 
But it should be noted that this will take a lot of time, so the files that take the longest to be considered are uploaded here.

## Note about big files
Unfortunately, some files were too big to upload them to Github, so all necessary data can be found [here](https://drive.google.com/drive/folders/1_i4X-JG8_BPvt-cC5FvP6oNdhgCIrxrF?usp=sharing). Some of the files are archived, to unzip them you can use the command  `gzip -d <filename>`. 

Rem: to correct notebooks work the file structure should be like this: 

![image](https://user-images.githubusercontent.com/83416875/169665459-c559dcb0-6d0f-4e42-b5da-3d60a1116c65.png)

# Main idea of new branchpoint annotator model
The main idea of the model is to use k-mers and distance to 3' ends to predict the position of branchpoints. The reasonableness of such an attempt lies in the fact that a polypyrimidine tract is usually located near the branchpoint (see left picture below), and the branchpoint itself is located at a short distance from the 3' end (see right plot). Moreover, such a model may take into account some other important positions in the intron, which we simply cannot assume.

![image](https://user-images.githubusercontent.com/83416875/169663848-e3b34022-a942-4391-aaa6-c6308f854166.png)
![image](https://user-images.githubusercontent.com/83416875/169662839-95ae92b2-6a7b-45e4-91e8-d229a1356515.png)


The model has ML-architecture with following structure: 

## Determine possible branchpoints using k-mers in (position-n, position+n) area of all position.

To train the model, we need to use both positive and negative examples of branchpoints
* As positive examples (label=1) we used experimental high confidence branchpoints
* As negative examples (label=0) we used random positions from introns that are not high confidence branchpoints

Note: we tried to use several variants of negative examples in the course of work: negative positions in introns near by 5' end or 3' end and random negative positions in introns, but best results the model showed at random negative positions (for example, if we learn model on negative positions near 5' end, it learned only distance to 3' end and no k-mers)

As a result, we obtain the table for training classifier like this: 

![image](https://user-images.githubusercontent.com/83416875/169663868-e919b119-548f-48e6-9ab7-58aa81a9cd81.png)![image](https://user-images.githubusercontent.com/83416875/169662672-22ee89e5-6122-4bd1-a9e5-edb22fecf69b.png)

Than we train classifier and then we can predict branchpoints. 

## Predict branchpoints
To predict positions of branchpoints looked around latest 50 positions of every intron (near to 3’ end) and find position with largest probability of label 1 (branchpoint)

![image](https://user-images.githubusercontent.com/83416875/169663304-dd6c55d9-b069-4708-aa09-1f56ae124151.png)![image](https://user-images.githubusercontent.com/83416875/169663289-acc60a38-c54b-4d85-b2c7-b163404dbf27.png)

On the right plot there is the probability of being a branchpoint for all introns near the 3' end. There is a probability peak, which in our model is called a branchpoint.

# Description of notebooks and results
## 1. Model development and parameters selection

This notebook is dedicated to choosing the best model parameters and the best classifier as well as comparison with existing models. 

### Parameters selection
The main result of the notebook is the determination of the best parameters and classifier: 
* n=70
* k=5
* classifier = Random Forest

The graphs on the basis of which this conclusion was made are presented below.

![image](https://user-images.githubusercontent.com/83416875/169662464-e3b21d01-b9dd-43ec-bdb2-e563132dfde3.png)

Moreover, the notebook shows that other classifiers like SVM and Naive Bayes have performed worse results. 

### Precalculations
After the model parameters are selected, tables with k-mers and their labels can be precalculated to save time. The resulting files are too large for Github, so they are located [here](https://drive.google.com/drive/folders/1_i4X-JG8_BPvt-cC5FvP6oNdhgCIrxrF?usp=sharing) as archived file `Xy for training model (n=70, k=5).tar.gz`, and files in it are compressed too.


### Comparison with other annotators
According to paper [2] there are roc-auc scores (see left plot below) for several BP-annotators like: branchpointer, HSF, Naive Bayes, SVM-BPFinder. We compared the metrics of our dash model with them and found that our results are better (see right boxplot).

![image](https://user-images.githubusercontent.com/83416875/169663761-8de31769-4317-43cc-ade7-66663219bc74.png) ![image](https://user-images.githubusercontent.com/83416875/169663604-51ad4a90-61aa-4085-877f-5c849e4fb4d1.png)


## 2. Branchpoint search evaluation
The goal of this notebook was to learn how to predict branch points in introns and compare predictions with known high confidence branchpoints from [1].

![image](https://user-images.githubusercontent.com/83416875/169664121-45250d78-fc9e-4b34-9759-3f60c04c9e23.png)

The histogram above demonstrates distribution of distance between high confidence branchpoint and pedicted one. As we can see, there are many coincidences , while the spread is not very large. So, the model often predict right positions. 

## 3. Branchpoint search for chromosome 
This notebook is devoted to predict branchpoints in all introns in all chromosomes. The results are presented in folder `predicted BPs for all chromosomes`.

Note: to find branchpoints in all introns in particular chromosome (for example, chromosome 1) you can run notebook with command 
```
#!/bin/bash
chrom="chr1"
papermill 3.\ Branchpoint\ search\ for\ chromosome.ipynb search_$chrom.ipynb --autosave-cell-every 60 --no-progress-bar --no-log-output -p chrom $chrom &
disown
```
or run the same python script 
```
python Banchpoint_search_in_chromosome.py chr1
```

# Future plans
Thus, we got a model that can determine the positions of branch points with good accuracy, however, to further improve it, new features can be added. Next, we want to implement an annotator similar to BBP and a branchpointer that works with vcf files and predicts the pathogenicity of mutations, and compare the results with the second part of our project.

# References
[1] Mercer, T. R., Clark, M. B., Andersen, S. B., Brunck, M. E., Haerty, W., Crawford, J., Taft, R. J., Nielsen, L. K., Dinger, M. E., & Mattick, J. S. (2015). Genome-wide discovery of human splicing branchpoints. Genome research, 25(2), 290–303. https://doi.org/10.1101/gr.182899.114

[2] Signal, B., Gloss, B. S., Dinger, M. E., & Mercer, T. R. (2018). Machine learning annotation of human branchpoints. Bioinformatics (Oxford, England), 34(6), 920–927. https://doi.org/10.1093/bioinformatics/btx688 
