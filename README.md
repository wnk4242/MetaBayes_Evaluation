# A Framework for Evaluating Replication-Success Metrics

This repository accompanies the manuscript entitled *“A Framework for Evaluating Replication-Success Metrics”*  
It contains the R code required to reproduce the analyses and visualizations presented in the paper.
Readers may find it helpful to explore our companion [Shiny app](https://quantpsych.shinyapps.io/BJMSP_AnatomyofReplicationSuccess/), an interactive diagram that visualizes the main findings of this study, while reading the manuscript.
<p align="center">
  <img src="images/RepSuccess_StackBar_nonnull.png" width="1345" style="margin-left: 20px;" />
</p>

---

### Simulation and Analysis Workflow

**Step 1.0–1.2: Data Generation**  
- Generate original study results.  
- Generate replication study results.  
- Synthesize replication data using meta-analytic Bayes factor (MABF) and random-effects meta-analysis (REMA) methods.  

**Step 2.0-2.1: Data Combination**  
- Combine individual simulation outputs into a single dataset for each of the MABF method as well as the REMA method for subsequent analysis.  
- Convert two-sided _p_-values into one-sided ones for meta-analysis in Step 2.1

**Step 3.0–3.5: Data Analysis**  
- Analyze original study results.  
- Compute classical classification performance metrics (TPR, FPR, etc) and replicaiton-success classification metrics (TSR, FSR, etc)

**Step 4.0: Data Visualization**  
- Visualize original study and replication study results using various types of plots.

