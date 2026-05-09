# A Framework for Evaluating Replication-Success Metrics

This repository accompanies the manuscript entitled *“A Framework for Evaluating Replication-Success Metrics”*  
It contains the R code, processed data, and supplementary materials required to reproduce the analyses and visualizations presented in the paper.
Readers may find it helpful to explore our companion [Shiny app](https://quantpsych.shinyapps.io/BJMSP_AnatomyofReplicationSuccess/), an interactive diagram that visualizes the main findings of this study, while reading the manuscript.
<p align="center">
  <img src="images/RepSuccess_StackBar_nonnull.png" width="1345" style="margin-left: 20px;" />
</p>

---

## Repository Contents

### R Code
The repository provides original R code used for simulation, data analysis, and visualization.

- The R code in the `original_code` folder reproduces the full simulation workflow.  
  Due to computational demands, R code for **Steps 1.0-3.5** require high-performance computing (HPC) resources and cannot be executed on a standard computer.  
- Fully processed simulation datasets are provided for readers to reproduce the main plots in the manuscript, including ROC curves and stacked bar charts.  
  The demonstration scripts for these analyses are located in the `demo code` folder.

---

### Simulation and Analysis Workflow

**Step 1.0–1.2: Data Generation**  
- Generate original study results.  
- Generate replication study results.  
- Synthesize replication data using meta-analytic Bayes factor (MABF) and random-effects meta-analysis (REMA) methods.  

**Step 2: Data Combination**  
- Combine individual simulation outputs into a single dataset for each of the MABF method as well as the REMA method for subsequent analysis.  

**Step 3.0–3.5: Data Analysis**  
- Analyze original study results.  
- Compute evaluation metrics for ROC curves and stacked bar plots.  

**Step 4: Data Visualization**  
- Visualize original study results using bar plots, histograms, and heatmaps.  
- Visualize replication study results using ROC-like curves and stacked bar plots.      
