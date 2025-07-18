#install.packages("DiagrammeR") 
library(DiagrammeR)

#Vertical flow

grViz("
digraph flowchart {
  graph [layout = dot, rankdir = TB]  # Top to Bottom for better left alignment
  
  node [shape=rect, style=filled, fontname=Arial, fontsize=12, width=2.5, height=0.8]

  # Define Nodes
  A [label='Step 1.0\\nOGDG.R', fillcolor=lightblue, align=left]
  B [label='Step 1.1\\nREPDG.R', fillcolor=lightblue, align=left]

  C1 [label='Step 1.2\\nBFbMA_Synth.R', fillcolor=lightblue, align=left]
  C2 [label='Step 1.2\\nEUBF_Synth.R', fillcolor=lightblue, align=left]
  C3 [label='Step 1.2\\nFEMABF_Synth.R', fillcolor=lightblue, align=left]
  C4 [label='Step 1.2\\niBF_Synth.R', fillcolor=lightblue, align=left]
  C5 [label='Step 1.2\\nFEMA_Synth.R', fillcolor=lightblue, align=left]

  D1 [label='Step 2.0\\nBFbMA_Data_Combination.R', fillcolor=green, align=left]
  D2 [label='Step 2.0\\nEUBF_Data_Combination.R', fillcolor=green, align=left]
  D3 [label='Step 2.0\\nFEMABF_Data_Combination.R', fillcolor=green, align=left]
  D4 [label='Step 2.0\\niBF_Data_Combination.R', fillcolor=green, align=left]
  D5 [label='Step 2.0\\nFEMA_Data_Combination.R', fillcolor=green, align=left]

  E1 [label='Step 3.0\\nBFbMA_Data_Analysis.R', fillcolor=lightcoral, align=left]
  E2 [label='Step 3.0\\nEUBF_Data_Analysis.R', fillcolor=lightcoral, align=left]
  E3 [label='Step 3.0\\nFEMABF_Data_Analysis.R', fillcolor=lightcoral, align=left]
  E4 [label='Step 3.0\\niBF_Data_Analysis.R', fillcolor=lightcoral, align=left]
  E5 [label='Step 3.0\\nFEMA_Data_Analysis.R', fillcolor=lightcoral, align=left]

  F1 [label='Step 4.0\\nBFbMA_ROC.R', fillcolor=plum, align=left]
  F2 [label='Step 4.0\\nEUBF_ROC.R', fillcolor=plum, align=left]
  F3 [label='Step 4.0\\nFEMABF_ROC.R', fillcolor=plum, align=left]
  F4 [label='Step 4.0\\niBF_ROC.R', fillcolor=plum, align=left]
  F5 [label='Step 4.0\\nFEMA_ROC.R', fillcolor=plum, align=left]

  G [label='Step 4.1\\nROC_AUC.R', fillcolor=plum, align=left]

  # Define Edges
  A -> B
  B -> C1; B -> C2; B -> C3; B -> C4; B -> C5;
  C1 -> D1; C2 -> D2; C3 -> D3; C4 -> D4; C5 -> D5;
  D1 -> E1; D2 -> E2; D3 -> E3; D4 -> E4; D5 -> E5;
  E1 -> F1; E2 -> F2; E3 -> F3; E4 -> F4; E5 -> F5;
  F1 -> G; F2 -> G; F3 -> G; F4 -> G; F5 -> G;
}
")



#Horizontal flow

grViz("
digraph flowchart {
  graph [layout = dot, rankdir = LR]  # Left to Right flow

  node [shape=rect, style=filled, fontname=Arial, fontsize=12, width=2.5, height=0.8]

  # Define Nodes
  A [label='Step 1.0\\nOGDG.R', fillcolor=lightblue]
  B [label='Step 1.1\\nREPDG.R', fillcolor=lightblue]

  C1 [label='Step 1.2\\nBFbMA_Synth.R', fillcolor=lightblue]
  C2 [label='Step 1.2\\nEUBF_Synth.R', fillcolor=lightblue]
  C3 [label='Step 1.2\\nFEMABF_Synth.R', fillcolor=lightblue]
  C4 [label='Step 1.2\\niBF_Synth.R', fillcolor=lightblue]
  C5 [label='Step 1.2\\nFEMA_Synth.R', fillcolor=lightblue]

  D1 [label='Step 2.0\\nBFbMA_Data_Combination.R', fillcolor=green]
  D2 [label='Step 2.0\\nEUBF_Data_Combination.R', fillcolor=green]
  D3 [label='Step 2.0\\nFEMABF_Data_Combination.R', fillcolor=green]
  D4 [label='Step 2.0\\niBF_Data_Combination.R', fillcolor=green]
  D5 [label='Step 2.0\\nFEMA_Data_Combination.R', fillcolor=green]

  E1 [label='Step 3.0\\nBFbMA_Data_Analysis.R', fillcolor=lightcoral]
  E2 [label='Step 3.0\\nEUBF_Data_Analysis.R', fillcolor=lightcoral]
  E3 [label='Step 3.0\\nFEMABF_Data_Analysis.R', fillcolor=lightcoral]
  E4 [label='Step 3.0\\niBF_Data_Analysis.R', fillcolor=lightcoral]
  E5 [label='Step 3.0\\nFEMA_Data_Analysis.R', fillcolor=lightcoral]

  F1 [label='Step 4.0\\nBFbMA_ROC.R', fillcolor=plum]
  F2 [label='Step 4.0\\nEUBF_ROC.R', fillcolor=plum]
  F3 [label='Step 4.0\\nFEMABF_ROC.R', fillcolor=plum]
  F4 [label='Step 4.0\\niBF_ROC.R', fillcolor=plum]
  F5 [label='Step 4.0\\nFEMA_ROC.R', fillcolor=plum]

  G [label='Step 4.1\\nROC_AUC.R', fillcolor=plum]

  # Define Edges
  A -> B
  B -> C1; B -> C2; B -> C3; B -> C4; B -> C5;
  C1 -> D1; C2 -> D2; C3 -> D3; C4 -> D4; C5 -> D5;
  D1 -> E1; D2 -> E2; D3 -> E3; D4 -> E4; D5 -> E5;
  E1 -> F1; E2 -> F2; E3 -> F3; E4 -> F4; E5 -> F5;
  F1 -> G; F2 -> G; F3 -> G; F4 -> G; F5 -> G;
}
")

