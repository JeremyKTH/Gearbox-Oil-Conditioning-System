# Gearbox-Oil-Conditioning-System
This is a repo for our Master Thesis. The code mainly done in Matlab and TwinCAT 3 through jupyter Lab.


## Table of contents

<!--ts-->
   * [Introduction](#Introduction)
   * [JupyterLab Folder](#JupyterLab_Folder)
   * [Matlab_Folder](#Matlab_Folder)
   * [Simulink_Folder](#Simulink_Folder)
   * [TwinCAT_Folder](#TwinCAT_Folder)
   * [Practice_Folder](#Practice_Folder)

<!--te-->

<!-- ABOUT THE PROJECT -->
## Introduction
This master thesis contains two research questions:

RQ1: Is about the investigation of system identification of a Gearbox Oil Conditioning System.

RQ2: Is about the further controller implementation using pole placement methods via TwinCAT3.

5 different linear system identification models were used: ARX, ARMAX, BJ, TF, SS.

## JupyterLab_Folder
The purpose for the folder is for data pre-processing. Namely, reorder and extract the useful data information from the original recorded data from the testbed. There are three processed dataset: (SI_47414, SI_48414, SI_49404), which represents the different valve operating range. In each dataset, it is then further splitted into training and testing dataset with size ratio 1:1.
* Data Pre-process
* SI_47414 (± 30% operating range )
* SI_48404 (± 40% operating range )
* SI_49404 (± 50% operating range )

## Matlab_Folder
* Best_30
* Best_40
* Best_50
* Controller_ARMAX
* Controller_SS
* Controller_TF

## Simulink_Folder
* Simulink_ARMAX
* Simulink_SS
* Simulink_TF

## TwinCAT_Folder
* PI
* PI + Pole placement
* PID + LP
* PID + LP + Pole placement

## Practice_Folder
The purpose for this folder is mainly for archiving past practicing files

<!-- CONTACT -->
## Contact

Chieh-Ju Wu (Jeremy) - jeremy.cjwukth@gmail.com
