# Spatial statistics of submicron size clusters of activating and inhibitory Natural Killer cell receptors in the resting state regulate early time signal discrimination

This project includes in-silico signaling modeling and information theory and image data for activating (NKG2D) and inhibtory (KIR2DL1).  

## 📂 Data Folder (`data/`)  

The `data/` directory contains all the datasets used in this case study. It includes raw from references (doi:10.1126/scisignal.aal3606 and doi: 10.1016/j.celrep.2016.04.075), processed after randomization. 

## 📂 src Folder (`src/`)  

The `src/` directory contains spparks code to Monte-Carlo simulation. 

Remianing folder are to run the code for different stages of Kintic proof reading which futher can be analysed after generating data at later time.
---
## 📂 Repository Structure  

```plaintext
├── data/               # Input datasets used for the case study
├── src/                # spparks code for Monte-Carlo simulation
├── Analysis/           # Jupyter notebooks for visualization
├── Kp_Both/            # Kintic proof reading present on both activating and inhibitory receptor ligand 
├── No_Kp/              # Kintic proof reading absent on both activating and inhibitory receptor ligand
├── No_Kp_Act/          # Kintic proof reading absent on activating receptor ligand
├── No_Kp_Inh/          # Kintic proof reading absent on inhibitory receptor ligand 
├── requirements.txt    # Python dependencies
├── README.md           # This file
└── LICENSE             # License file
'''
--


