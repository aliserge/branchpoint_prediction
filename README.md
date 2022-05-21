# Beyond standard genetic variant annotation: branchpoints

Authors:
* Sergeeva Alisa
* Veretenenko Irina

# Introduction
A key step in the molecular diagnosis of hereditary diseases is the interpretation of genetic variants found during genome or exome sequencing. At the moment, however, interpretation is quite a challenge, and the cumulative efficiency of molecular diagnostics ranges from 30 to 50%. This is partly due to the lack of good methods for annotating different classes of genetic variants that do not involve a violation of the amino acid sequence of a protein.

This repository is dedicated to solving one of the tasks of a large project on the annotation of genetic variants. It is dedicated to predicting the effects of genetic variants on intra-intron sequences involved in the passage of splicing (primarily, in the region of the branchpoint). 

## About branchpoints
![image](https://user-images.githubusercontent.com/83416875/169659518-e4e0b883-ed98-42cb-b183-96b81bafe59d.png)

The above picture shows the mechanism of splicing, where 5' end of the intron connects to branchpoint (BP) and only after that 3' end cuts of exon. There are several splicing mechanisms that differ in the use of proteins, but in general it works as follows: 
* The 2'-OH group A of the BP attacks the phosphodiester bond at the 5' end of the intron and cleaves it. In this case, an unusual 2'->5'-bond is formed inside the * intron, due to which the intron takes the form of a loop.
* The free 3'-OH group at the end of the 5' end of the exon attacks the A–G bond at the 3' end of the intron.
* This results in the binding of the two exons and the release of the intron, still in the form of a loop.

Thus, a mutation in a branchpoint can disrupt splicing and lead to pathological consequences; therefore, predicting the effects of variants in branchpoints is an urgent task for diagnosing various diseases.

## Task and workflow

The aim of our work was to determine the effects (including the assessment of pathogenicity) of variants located in branch points.

In the course of work, the project was divided into two main parts:
* Variants annotation based on existing annotators (BPP_Branchpointer)
* New predictive model development (BP prediciton)

Description and current results of work on these subtasks are presented in the corresponding folders.

