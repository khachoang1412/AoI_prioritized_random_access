# Age of Information in Prioritized Random Access

This repository contains the matlab numerical routines in the paper:

K.-H. Ngo, G. Durisi, and A. Graell i Amat, "Age of information in prioritized random access," in 55th Asilomar Conference on Signals, Systems, and Computers, CA, USA, Oct. 2021.

Please, cite the aforementioned paper if you use this code.

## Content of the repository

This repository contains three folders:

1. `/IRSA_implementation/`: an implementation of irregular repetition slotted AHOLA (IRSA) and a computation of the packet loss rate (PLR) via Monte-Carlo simulation. 

2. `/PLR_approximation/`: Various approximations of the PLR:
  - `PLR_DE_singleClass.m`, `PLR_DE_mutiClass.m`: asymptotic PLR as the framelength goes to infinity, computed via density evolution.
  - `PLR_errfloor_singleClass.m`, `PLR_errfloor_mutiClass.m`: PLR approximation in the error floor region, according to Ivanov (2016) and Ivanov (2017).
  - `PLR_waterfall.m`: PLR approximation in the waterfall region, according to Graell i Amat (2018).
  - `PLR_approx_multiClass.m`: PLR approximation proposed in Section IV of Ngo (2021).

3. `/AoI_analysis_optimization/`
  - `eval_PLR_AVP.m`: evaluate the PLR and AoI for given degree distributions and other parameters. This file is used to generate the data plotted in the figures in Ngo (2021).
  - `minimize_power.m`: optimize the degree distributions and the update probability to mimize the average number of transmitted packets per slot (as a proxy to minimize power consumption) while satisfying AoI constraints.
