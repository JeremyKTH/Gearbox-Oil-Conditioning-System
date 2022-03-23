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

![GOCS](https://github.com/JeremyKTH/Gearbox-Oil-Conditioning-System/blob/main/Images/GOCS.jpg)

## JupyterLab_Folder
The original recoreded data from the testbed are cleaned and the useful information are extracted into three separate dataset: (D_47414, D_48414, D_49404), which represents the system data at different valve operating range. In each dataset, it is then further splitted into training and testing dataset with size ratio 1:1.

SI_47414, SI_48414, SI_49404 are the system identification programs developed with JupyterLab in language Matlab for easier reading.

* Data Pre-process
* SI_47414 (± 30% operating range )
* SI_48404 (± 40% operating range )
* SI_49404 (± 50% operating range )

## Matlab_Folder
Best_30, Best_40, Best_50 contains all the best model structure for each system identification models at ± 30%, ± 40%, ± 50% operating range. The best 3 system identification models based on simulation accuracy are selected for research question 2. (ARMAX, SS, TF)

Controller_ARMAX, Controller_SS, Controller_TF contains the controllers developed upon the top 3 models from research question 1.

* Best_30
* Best_40
* Best_50
* Controller_ARMAX
* Controller_SS
* Controller_TF

## Simulink_Folder
The simulink files are used for initial simulation of the system identification models as well as inspecting the controller performances.
* Simulink_ARMAX
* Simulink_SS
* Simulink_TF

## TwinCAT_Folder
The developed controllers from Matlab are transformed into compatible version for TwinCAT3 environment.
* PI
* PI + Pole placement
* PID + LP
* PID + LP + Pole placement

## Practice_Folder
The purpose for this folder is mainly for archiving past practicing files.

<!-- CONTACT -->
## Contact

Chieh-Ju Wu (Jeremy) - jeremy.cjwukth@gmail.com
