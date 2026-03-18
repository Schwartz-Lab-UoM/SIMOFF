# **SIMOFF**
SIMulated annealing Objective Function Finder (SIMOFF)
## **Examples**
Example applications are provided for Chinese Hamster Ovary (CHO) cells and Yeast cells. Examples are Jupyter notebooks. Examples can be found in the 'examples' folder.
## **Implementation**
To use SIMOFF with your own data, download the SIMOFF.py file (from 'src' folder) and save this into the same folder as your application notebook. Within the application notebook, type 'import SIMOFF as smf', then smf.simoff() with the correct arguments specified within brackets will allow you to run the SIMOFF function. Ensure all appropriate requirements (requirements.txt) file have been installed and imported into notebook.
SIMOFF has not yet been made into a Python package format. 
## **Table of contents** #
```
├── examples/                                         # Folder containing applications/tutorials using SIMOFF
|  ├── Tutorial_CHO_late_exponential.ipynb            # Tutorial for SIMOFF on a late exponential culture phase CHO model
|  ├── Tutorial_Yeast_faster_grower.ipynb             # Tutorial for SIMOFF on a faster-growing yeast condition
|  ├── Tutorial_Yeast_slower_grower.ipynb             # Tutorial for SIMOFF on a slower-growing yeast condition
├── figures/                                          # Folder with code to generate figures for SIMOFF publication (not yet submitted/accepted)
│   ├── CHO_early_exponential_phase_figs1a_b_2a_b_c_d.ipynb
│   ├── CHO_late_exponential_phase_fig6.ipynb
│   ├── CHO_late_exponential_phase_figs1c_d_2f.ipynb
│   ├── CHO_stationary_phase_figs1e_f.ipynb
│   ├── Yeast_fig_4_equations_1_2.ipynb
│   ├── Yeast_fig_5_a_b_c_d.ipynb
│   ├── Yeast_figs5e_f.ipynb
├── src/                                              # Folder containing source code for SIMOFF project/publication
│   ├── Mismatch_grid.py                              # Create a heatmap-style grid for agreement/disagreement between SIMOFF predictions and experimental criteria
│   ├── Plot_coefficients.py                          # Plot the trajectory of the weighted reaction coefficients over the SIMOFF search
│   ├── SIMOFF.py                                     # SIMOFF function
│   ├── Upset_plot.py                                 # Create Upset-plot style plot to show how well each combination of input reactions explored matched experimental criteria
│   ├── _init_.py                                     # Initialises packages when imported
├── README.md
├── requirements.txt
```
## **Acknowledgements** ##
The authors acknowledge financial support from a Prosperity Partnership grant (EP/V038095/1) funded by EPSRC, BBSRC, and FUJIFILM Diosynth Biotechnologies. We also acknowledge everyone who has helped interpret experimental data for input datasets; those groups that developed the GEMs that are tested in 'Examples' and the more general systems biology field for their support and feedback.
## **FAQ/Troubleshooting** ##
