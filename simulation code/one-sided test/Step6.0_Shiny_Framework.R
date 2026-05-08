# Note: I2 is set to 30% in this stacked bar shiny app. In the Step6.0_Shiny_Framework_Curves.R, I2=10%, so expect to see different values for proportions.
library(shiny)
library(plotly)
library(dplyr)
library(shinyjs)

ui <- fluidPage(
  useShinyjs(),
  
  tags$head(
    tags$style(HTML("
            .MathJax_Display {
      text-align: left !important;
      margin-left: 0 !important;
        }
    
/* ===== PAGE BACKGROUND ===== */
body {
  background-image: url('https://www.transparenttextures.com/patterns/cubes.png');
  background-size: auto;
  background-repeat: repeat;
  background-color: #fafafa;
}

/* ===== SIDEBAR DESIGN ===== */
.sidebar {
  background-image: url('https://www.transparenttextures.com/patterns/cubes.png');
  background-color: rgba(255,255,255,0.90);
  padding: 15px;
  border-radius: 5px;
  max-height: 85vh;
  overflow-y: auto;
}

/* Hide the scrollbar */
.sidebar::-webkit-scrollbar {
  width: 8px;
}
.sidebar::-webkit-scrollbar-track {
  background: transparent;
}
.sidebar::-webkit-scrollbar-thumb {
  background-color: transparent;
}

/* ===== MAIN PANEL ===== */
.main-panel {
  padding: 20px;
  background-color: rgba(255,255,255,0.80);
  border-radius: 5px;
}

/* ===== TABSET PANEL STYLING ===== */
.nav-tabs {
  border-bottom: 1px solid transparent;
}

.nav-tabs > li > a {
  background-image: url('https://www.transparenttextures.com/patterns/cubes.png');
  background-color: rgba(255,255,255,0.6);
  border: 1px solid #dcd7d7;
  background-blend-mode: multiply;
  color: black;
  font-weight: bold;
  padding: 10px 20px;
  border-radius: 5px 5px 0 0;
  transition: all 0.15s ease-in-out;
}

/* ACTIVE TAB */
.nav-tabs > li.active > a,
.nav-tabs > li.active > a:focus,
.nav-tabs > li.active > a:hover {
  background-color: rgba(255,255,255,1);
  border: 1px solid #999;
  color: black;
  z-index: 10;
  position: relative;
  top: -2px;
}
        .js-irs-0 .irs-min, .js-irs-0 .irs-max, .js-irs-0 .irs-single {
          display: none !important;
        }
        .bias-label {
          background-color: #d6d6d6;
          padding: 1px 6px;
          border-radius: 4px;
          font-size: 12px;
          color: #222;
        }

        /* ============================
           TRANSPARENT PLOT WINDOWS
           ============================ */

        /* Transparent plotly background */
        .plotly-container,
        .plotly,
        .plot-container,
        svg.main-svg {
          background: transparent !important;
        }

        /* Remove white fill inside axes / panel */
        .plotly .cartesianlayer,
        .plotly .bglayer rect,
        svg .bglayer rect {
          fill: transparent !important;
        }

        /* Make Shiny output boxes transparent */
        #replicationPlot,
        #replicationPlot2,
        #categoryDiagram {
          background-color: transparent !important;
          border: none !important;
          box-shadow: none !important;
        }

  "))
  ),
  
  titlePanel("Evaluation of Multi-Lab Replication Outcomes"),
  
  sidebarLayout(
    sidebarPanel(
      tags$style(HTML("
        .js-irs-0 .irs-min, .js-irs-0 .irs-max, .js-irs-0 .irs-single {
          display: none !important;
        }
        .bias-label {
          background-color: #d6d6d6;
          padding: 1px 6px;
          border-radius: 4px;
          font-size: 12px;
          color: #222;
        }
      ")),
      
      
      # --- Bias level slider ---
      div(style = "margin-top:10px;",
          
          # Title WITH tooltip
          div(
            style = "text-align:center; margin-bottom:5px;",
            HTML('
        <span title="Combined effect of p-hacking and publication bias on original studies">
          <h5>Bias Level <br>
          (p-hacking and publication bias)</h5>
        </span>
      ')
          ),
          
          # Labels WITHOUT tooltips
          div(
            style = "display:flex; justify-content:space-between;
         margin-top:-6px;
         margin-bottom:-10px;
         font-size:11px;",
            
            div("No bias",
                style = "background-color:#e0e0e0; padding:1px 6px;
             border-radius:4px; color:#222;"
            ),
            
            div("High bias",
                style = "background-color:#e0e0e0; padding:1px 6px;
             border-radius:4px; color:#222;"
            )
          ),
          
          # Slider
          
          
          
          sliderInput(
            "bias", label = NULL,
            min = 0, max = 1, value = 0.8,
            step = 0.025, ticks = TRUE
          )
      ),
      # Underlying effect size (non-null only)
      withMathJax(
        div(style = "margin-top:25px;",
            div(
              style = "text-align:center; margin-bottom:5px;",
              HTML('
              <span title="Population effect for both original and replication studies, measured in standardized mean difference">
                   Underlying effect size <span style="font-size:85%;">(\\(\\theta\\))</span>')
            ),
            sliderInput("effect_size", NULL, min=0.1, max=0.5, value=0.1, step=0.05)
        )
      ),
      
      # --- Alpha slider ---
      withMathJax(
        div(style = "margin-top:25px;",
            div(
              style = "text-align:center; margin-bottom:5px;",
              HTML('
              <span title="Nominal significance level for the orignal study">
                   Nominal Significance Level <span style="font-size:85%;">(\\(\\alpha\\))</span>')
            ),
            sliderInput("alpha", NULL, min=0.01, max=0.05, value=0.05, step=0.01)
        )
      ),
      
      
      
      
      div(style = "margin-top:10px;",
          
          # Title with tooltip
          div(
            style = "text-align:center; margin-bottom:5px;",
            HTML('
        <span title="Number of participants per group in the original study">
          Original Study Sample Size<span style="font-size:85%;"> (\\(N_{orig}\\))
        </span>
      ')
          ),
          
          # Slider
          sliderInput(
            "orig_n", label = NULL,
            min = 20, max = 200, value = 20,
            step = 10, ticks = TRUE
          )
      ),
      
      
      # --- Number of Replications ---
      div(style = "margin-top:25px;",
          
          # Title using HTML + MathJax (allows styling)
          div(
            style = "text-align:center; margin-bottom:1px;",
            HTML('
            <span title="Number of direct replications of the original study">
        <span>
          Number of Replications 
          <span style="font-size:85%;">(\\(N_{rep}\\))</span>
        </span>
      ')
          ),
          
          sliderInput(
            "rep_num", label = NULL,
            min = 2, max = 10, value = 2,
            step = 1, ticks = TRUE
          )
      ),
      
      
      # --- Replication Sample Size ---
      div(style = "margin-top:25px;",
          
          # Title WITH tooltip
          div(
            style = "text-align:center; margin-bottom:1px;",
            HTML('
        <span title="Number of participants per group in the replication study">
          Replication Sample Size <span style="font-size:85%;">(\\(n_{rep}\\))
        </span>
      ')
          ),
          
          # Slider
          sliderInput(
            "rep_n", label = NULL,
            min = 40, max = 400, value = 40,
            step = 10, ticks = TRUE
          )
      ),
      
      div(
        style = "margin-top:30px; text-align:left; color:#4a8dc4; font-size:14px;",
        HTML('
    &copy; 2025 <a href="https://wnk4242.github.io" target="_blank" style="color:#4a8dc4; text-decoration:none;">Naike Wang, PhD</a><br>
    Last updated November, 2025
  ')
      ),
      width = 2
    ),
    
    mainPanel(
      width = 10,
      fluidRow(
        tabsetPanel(
          id = "tabs", 
          tabPanel(
            "Null Effect", value = "null",
            
            fluidRow(
              # LEFT column: Plot only
              column(
                width = 7,
                plotlyOutput("replicationPlot", height = "700px")
              ),
              
              # RIGHT column: All explanatory text
              column(
                width = 5,
                div(style = "margin-top: 18px; font-size: 14px; text-align: justify; max-width: 600px; padding-left: 100px;",
                    HTML('
              <p>
                This conceptual diagram summarizes the main findings of my study 
                <i><a href="https://osf.io/preprints/psyarxiv/a42c6_v1" target="_blank">
                    “Seeing Beyond Replication Success: A Framework for Evaluating Multi-Lab Replication Outcomes Using Meta-Analytic Bayes Factors.”
                  </a></i> 
                The simulation code, data files, and supplementary files are available in this 
                <a href="https://github.com/wnk4242/MetaBayes_Evaluation" target="_blank">
                  Github repository</a>.
              </p>


           <p>The study uses large-scale simulations to examine how several methodological factors (shown in the left sidebar) influence the replication outcomes of original experimental studies. Moving beyond the conventional binary classification of replication as simply a success or a failure, we further distinguish outcomes into true and false successes and true and false failures (detailed definitions are provided in the Study Overview tab).</p>
          
          <p>This diagram shows how changes in these study-level methodological factors influence replication outcomes when the underlying true effect is null or spurious. Some key observations include:</p>
          <ul>
            <li><b>Higher nominal significance levels</b> (e.g., \\( \\alpha = 0.05 \\)) combined with strong bias (e.g., p-hacking and publication bias) substantially increase the false positive rate in original studies. This leads to greater inconsistency with replication results and a larger proportion of replication failures—most of which are false failures.</li>
            <li><b>Lower nominal significance levels</b> (e.g., \\( \\alpha = 0.01 \\)) keep the false positive rate low even under high bias. As a result, original and replication results are more consistent, leading to more true replication successes and fewer failures.</li>
            <li><b>Replication outcomes are sensitive to the nominal significance level</b> used in the original studies. With \\( \\alpha = 0.05 \\), most failures are false failures. With \\( \\alpha = 0.01 \\), most successes are true successes.</li>
            <li><b>Increasing the original study\'s sample size</b> has only a small effect on replication outcomes under the assumption that the true effect is null. It slightly increases the proportion of false failures but does not change the outcome patterns substantially.</li>
          </ul>

        ')
                )
              )
            )
          ),
          
          tabPanel("Non-Null Effect", value = "nonnull",
                   
                   fluidRow(
                     # LEFT column: Plot only
                     column(
                       width = 7,
                       plotlyOutput("replicationPlot2", height = "700px")
                     ),
                     
                     # RIGHT column: All explanatory text
                     column(
                       width = 5,
                       div(style = "margin-top: 18px; font-size: 14px; text-align: justify; max-width: 600px; padding-left: 100px;",
                           HTML('
                           </p>Replication outcomes follow a more complex pattern when the underlying effect deviates from null, especially when it is a small positive one (e.g., <i>d</i> = 0.2). 
                           Some key observations include:</p>
          <ul>
            <li> <b>In the absence of bias (p-hacking and publication bias)</b>, a small original study using a nominal significance level of 
                  \\( \\alpha = 0.05 \\) will generate many false negatives. If replications also use small designs 
                  (few replications and/or small sample sizes), they will likewise produce many false negatives. 
                  Although these results appear consistent and are labeled as "replication successes," most of them are actually false successes. </li>
            <li> <b>Under strong bias</b>, the true positive rate in original studies becomes inflated, even when the original sample size is very small 
                (e.g., \\(N_{orig}=20\\)). When replications rely on small designs, they tend to 
                produce many false negatives. These inconsistencies lead to many failed replications, which in 
                this case are true replication failures. However, these true failures do not indicate that the underlying 
                effect is null; rather, they happen because the replications simply lack sufficient 
                statistical power.</li>
            <li> <b>When the bias is strong</b>, a large original study that uses a more lenient significance cutoff, together with replications based on large 
            sample sizes, will mostly yield true successes. However, the magnitude of the estimated effect should be interpreted cautiously: 
            under strong bias, original studies are more likely to produce inflated effect size estimates. Thus, the effect estimated by the replications
            are more trustworthy.</li>

          </ul>
        ')
                       )
                     )
                   )
          ),
          
          tabPanel("Instructions",
                   div(style = "margin: 15px; font-size: 14px;",
                       HTML("
                          
                          <h4><b>How to use this app</b></h3>
                          
                          <p>
                          This app illustrates how different study-level methodological factors influence replication outcomes
                          under two general situations:
                          </p>
                          
                          <ul>
                            <li><b>Null effect</b>: The underlying effect is null or spurious (θ = 0).</li>
                            <li><b>Non-null Effect</b>: The underlying effect is genuine and positive (θ > 0).</li>
                          </ul>
                          
                          <p>
                          The left sidebar contains several sliders. Each slider represents a study-level factor that affects 
                          the behavior of the original study, the replication studies, or both.  
                          The diagrams update automatically as you adjust these sliders.
                          </p>
                          
                          <br>
                          
                          <h4><b>Slider descriptions</b></h3>
                          
                          <h5><b>Bias level</b></h4>
                          <p>
                          The level of bias mechanism represents the degree of p-hacking and publication bias applied to the original study results. Higher values indicate stronger p-hacking and more severe publication bias.
                          </p>
                          
                          <ul>
                            <li><b>Null effect scenario:</b> Increasing bias inflates the original study\'s false positive rate.  
                                A higher level of bias leads to more <i>false successes</i> and <i>false failures</i> in the replication outcomes, keeping other variables constant.</li>
                          
                            <li><b>Non-null scenario:</b> Bias inflates the original study\'s true positive rate (but also distorts effect sizes).  
                                Stronger bias increases the chance of <i>true successes</i> or <i>true failures</i>, keeping other variables constant.</li>
                          </ul>
                      
                          
                          <h5><b>Underlying effect size (θ)</b></h4>
                          <p>
                          This slider is active only in the <b>Non-Null Effect</b> tab.  
                          It represents the population effect size in standardized mean difference.
                          </p>
                          
                          <ul>
                            <li><b>Null effect scenario:</b> This slider is disabled because the effect is fixed at θ = 0.</li>
                          
                            <li><b>Non-null scenario:</b> Increasing θ raises both the original and replication studies\' true positive rates, keeping other variables the same.  
                                When θ is small, both original and replication studies may fail to detect the effect, producing many <i>false successes</i>.  
                                As θ increases, true successes become more common, especially when the replication design is strong (i.e., more replications with larger sample sizes).</li>
                          </ul>
                          
                          
                          <h5><b>Nominal significance level (α)</b></h4>
                          <p>
                          This slider represents the nominal significance cutoff used in the <b>original study</b> (e.g., 0.01 or 0.05).  
                          </p>
                          
                          <ul>
                            <li><b>Null effect scenario:</b> A larger α increases the false positive rate in original studies.  
                                As the nominal significance level reaching α = 0.05, replication failures become more common (keeping the replication design the same), and most failures are <i>false failures</i>.</li>
                          
                            <li><b>Non-null scenario:</b> A larger α increases the true positive rate in original studies.
                                A smaller α value increases the false negative rate and shifts the replication outcomes toward <i>false successes</i> and <i>false failure</i> keeping other variables constant.
                          </ul>
                          

                          
                          
                          <h5><b>Original study sample size (N<sub>orig</sub>)</b></h4>
                          <p>
                          This slider controls how many participants per group are in the original study.  
                          </p>
                          
                          <ul>
                            <li><b>Null effect scenario:</b> Increasing sample size slightly raises the chance of obtaining borderline significant 
                                false positives, leading to more <i>false failures</i>.</li>
                          
                            <li><b>Non-null scenario:</b> Larger samples increase the original study\'s true positive rate, keeping other variables constant.  
                                Increasing original study's sample size leads to more <i>true successes</i> and <i>true failures</i>, keeping the replication design the same.</li>
                          </ul>
                          

                          
                          <h5><b>Number of replications (N<sub>rep</sub>)</b></h4>
                          <p>
                          This slider controls how many direct replication studies are conducted.  
                          </p>
                          
                          <ul>
                            <li><b>Null effect scenario:</b> More replications slightly reduce the false positive rate in replications because the false positive rate is already very low.  
                                </li>
                          
                            <li><b>Non-null scenario:</b> More replications make the meta-analytic Bayes factor more reliable, increasing the true positive rate in replicaitons.  
                                This increases <i>true successes</i> and <i>false failures</i> and reduces <i>false successes</i> and <i>true failures</i>.</li>
                          </ul>
                          

                          
                          <h5><b>Replication sample size (n<sub>rep</sub>)</b></h4>
                          <p>
                          This slider controls the number of participants per group used in each replication study.  

                          </p>
                          
                          <ul>
                            <li><b>Null effect scenario:</b> A larger replication sample size slightly reduces false positives as the false positive rate is already low.</li>
                          
                            <li><b>Non-null scenario:</b> Larger replication sample sizes make the meta-analytic Bayes factor more reliable, increasing the true positive rate in replicaitons.   
                                As the sample size increases, the replication study becomes much more accurate, increasing <i>true successes</i> and <i>false failures</i> and reducing <i>false successes</i> and <i>true failures</i>.</li>
                          </ul>
                          
                          <br>
                          
                          

                        <h4><b>Summary</b></h3>
                          <p>
                          Each slider controls a different methodological factor.  
                          Together, they determine the pattern of replication outcomes in both the null and non-null scenarios.
                          The diagrams help illustrate how different treatment effects, severity of bias, and changes in the design of original and replication studies can vary what “replication success” or “replication failure” actually means.
                          </p>"))
          ),
          
          
          tabPanel("Study Overview",
                   div(style = "margin: 15px; font-size: 14px;",
                       HTML("
        <b>Study overview</b><br>               
        This study introduces a new replication success classification framework for multi-lab replication studies that goes beyond 
        the traditional binary distinction of replication success and failure. The framework distinguishes between true and false 
        replication successes and failures, revealing the underlying composition of observed replication results and how they are 
        shaped by study-level factors in both the original and replication studies. Using a large-scale simulation, we evaluated four
        meta-analytic Bayes factor (MABF) methods as tools for determining replication success within this framework. Results showed
        that the MABF methods were generally more effective at ruling out spurious effects than at confirming small true effects. 
        The framework illustrates how bias and methodological choices influence the underlying composition of 
        replication outcomes. These findings provide a more nuanced understanding of replication results.

        <br><br>
        
        <b>Research questions</b><br>
        The study answers the following question, among others: <i>When original and replication findings are consistent (replication success), how likely are they both to be correct or incorrect?</i> and <i>When the original and replication findings are inconsistent (replication failure), how likely are the original results to be correct and how likely are the replication results to be correct?</i>
        <br><br>
        
        <b>Simulation design</b><br>
         In the simulation, original studies using a two-group between-subjects design are generated with varying nominal significance levels
         and sample sizes. A bias mechanism is introduced to impose different degrees of p-hacking and publication bias on these original 
         studies. Replication studies are then simulated using the same experimental design but with varying numbers of replications and 
         sample sizes, and they are not influenced by the bias mechanism. The presence of an effect in the original study is evaluated 
         using null-hypothesis significance testing, while the existence of the effect in the replications is evaluated using meta-analytic 
         Bayesian hypothesis testing through meta-analytic Bayes factors.
        
        <br><br>
         
        <b>Replication outcome classification</b><br>
        Replication success is defined by the consistency between the original study and replication results: 
        if they align, the replication is considered successful; otherwise, it is classified as unsuccessful. 
        Replication outcomes are categorized into four types: true success (TS), false failure (FF), true failure (TF), and false success (FS).
        This diagram illustrates how replication outcomes are classified as true or false successes and failures.
      ")
                   ),
                   
                   # Your diagram
                   plotlyOutput("categoryDiagram"),
                   
                   # More paragraphs below the diagram
                   div(style = "margin: 15px; font-size: 14px;",
                       HTML("
        <br>
        <br>
        <b>When the underlying effect is null:</b>
        <ul style='margin-top:4px;'>
          <li><b>True success:</b> Both the original and replication results consistently indicate the effect is null.</li>
          <li><b>False success:</b> Although the results are consistent, both incorrectly indicate the effect deviates from null.</li>
          <li><b>False failure:</b> The results are inconsistent. Only the replication correctly indicates the effect is null.</li>
          <li><b>True failure:</b> The results are inconsistent. Only the original study indicates the effect is null.</li>
        </ul>     
        <b>When underlying effect deviates from null (positive):</b>
        <ul style='margin-top:4px;'>
          <li><b>True success:</b> Both the original and replication results indicate the effect deviates from null.</li>
          <li><b>False success:</b> Although the results are consistent, neither indicates the effect deviates from null.</li>
          <li><b>False failure:</b> The results are inconsistent. Only the replication correctly indicates the effect deviates from null.</li>
          <li><b>True failure:</b> The results are inconsistent. Only the original study indicates the effect deviates from null.</li>
        </ul>   
        <br>
 <b>Mathematical relationships of the outcome classifications</b><br>
        The proportions in the diagrams reflect general (but not exact) patterns observed in simulation results. 
        The changes in the proportions are governed by several fundamental statistical rules and mathematical relationships:
        <br><br>
        <ul>
          <li>
            <b>Success-Failure decomposition</b><br>
            Successes = True successes + False successes<br>
            Failures = True failures + False failures
          </li>
          <br>
          <li>
            <b>When the true effect is null</b><br>
            False positives = False successes + True failures<br>
            True negatives = True successes + False failures
          </li>
          <br>
          <li>
            <b>When the true effect is non-null</b><br>
            True positives = True successes + False failures<br>
            False negatives = True failures + False successes
          </li>
        </ul>
      ")
                   ),
                   div(style = "margin: 15px; font-size: 14px;",
                       HTML("
              <b>Statistical derivation of outcome proportions</b><br>
              The outcome proportions shown in the diagrams are derived directly from the statistical models
              implemented in the simulation. This section documents the exact formulas used to compute the
              original and replication study outcomes.
              
              <br><br>
              
              <b>1. Original study outcomes</b><br>
              
              Let <i>&alpha;<sub>0</sub></i> denote the nominal significance level of the original study and
              <i>&theta;</i> denote the true standardized effect size.
              
              <br><br>
              
              <b>Nominal false positive rate</b><br>
              Under the null hypothesis and in the absence of bias, the false positive rate of the original
              study equals the nominal Type I error rate:
              <br><br>
              
              \\[
              FP_{orig, nominal} = \\alpha_0
              \\]
              
              <br><br>
              
              <b>Original study power (true effect present)</b><br>
              Let <i>n<sub>0</sub></i> denote the per-group sample size of the original study.
              The critical value for a one-sided z test is:
              <br><br>
              
              \\[
              z_{crit} = \\Phi^{-1}(1 - \\alpha_0)
              \\]
              
              <br><br>
              
              The nominal power of the original study is:
              <br><br>
              
              \\[
              TP_{orig, nominal}
              = 1 - \\Phi\\left(z_{crit} - \\theta \\sqrt{\\frac{n_0}{2}}\\right)
              \\]
              
              <br><br>
              
              The corresponding false negative rate is:
              <br><br>
              
              \\[
              FN_{orig} = 1 - TP_{orig}
              \\]
              
              <br><br>
              
              <b>Bias-adjusted original study outcomes (optional stopping)</b><br>
              Bias is introduced through optional stopping with probability <i>b</i>.
              When <i>b = 0</i>, original study outcomes are unbiased:
              <br><br>
              
              \\[
              TP_{orig} = TP_{orig, nominal}, \\quad FP_{orig} = \\alpha_0
              \\]
              
              <br><br>
              
              When <i>b &gt; 0</i>, original study outcomes are a mixture of unbiased and biased processes:
              <br><br>
              
              \\[
              TP_{orig} = (1 - b) \\cdot TP_{orig, nominal} + b \\cdot \\pi_{1,OS}
              \\]
              
              \\[
              FP_{orig} = (1 - b) \\cdot \\alpha_0 + b \\cdot \\pi_{0,OS}
              \\]
              
              <br><br>
              
              where <i>&pi;<sub>1, OS</sub></i> and <i>&pi;<sub>0, OS</sub></i> denote the probabilities of rejecting the null hypothesis under optional stopping when the
          true effect is non-null and null, respectively.
              
              <br><br>
              
              The true negative rate of the original study is:
              <br><br>
              
              \\[
              TN_{orig} = 1 - FP_{orig}
              \\]
              
              <br><br>
              
              <b>2. Replication study outcomes (meta-analysis)</b><br>
              
              Replication studies are evaluated using a two-sided meta-analytic hypothesis test with
              significance level <i>&alpha;<sub>meta</sub> = 0.05</i>. Bias affects only the original study;
              replication studies are unbiased.
              
              <br><br>
              
              <b>Meta-analytic variance structure</b><br>
              
              Let <i>n<sub>r</sub></i> denote the per-study replication sample size and <i>k</i> the number of
              replications. The within-study sampling variance is:
              <br><br>
              
              \\[
              v = \\frac{2}{n_r} + \\frac{\\theta^2}{4 n_r}
              \\]
              
              <br><br>
              
              Between-study heterogeneity is fixed by setting <i>I<sup>2</sup> = 0.30</i>:
              <br><br>
              
              \\[
              \\tau^2 = \\frac{0.30 \\cdot v}{1 - 0.30}
              \\]
              
              <br><br>
              
              The meta-analytic variance of the pooled estimate is:
              <br><br>
              
              \\[
              \\bar{v} = \\frac{v + \\tau^2}{k}
              \\]
              
              <br><br>
              
              <b>Meta-analytic power (true effect present)</b><br>
              
              The noncentrality parameter of the meta-analytic test statistic is:
              <br><br>
              
              \\[
              \\lambda = \\frac{\\theta}{\\sqrt{\\bar{v}}}
              \\]
              
              <br><br>
              
              The critical value for a two-sided test is:
              <br><br>
              
              \\[
              z_{crit, meta} = \\Phi^{-1}\\left(1 - \\frac{\\alpha_{meta}}{2}\\right)
              \\]
              
              <br><br>
              
              The probability that the meta-analysis detects a true effect is:
              <br><br>
              
              \\[
              TP_{rep}
              = 1 - \\Phi(z_{crit, meta} - \\lambda)
              + \\Phi(-z_{crit, meta} - \\lambda)
              \\]
              
              <br><br>
              
              The corresponding false negative rate is:
              <br><br>
              
              \\[
              FN_{rep} = 1 - TP_{rep}
              \\]
              
              <br><br>
              
              <b>Meta-analytic outcomes under the null</b><br>
              
              When the true effect is null (<i>&theta; = 0</i>), the meta-analytic false positive rate equals
              the nominal level:
              <br><br>
              
              \\[
              FP_{rep} = \\alpha_{meta}, \\quad
              TN_{rep} = 1 - \\alpha_{meta}
              \\]
              
              <br><br>
              
              <b>3. Joint replication outcome classification</b><br>
              
              Replication outcomes are classified by combining original and replication results through joint
              probabilities.
              
              <br><br>
              
              <b>When the true effect is non-null (&theta; &gt; 0):</b>
              <br><br>
              
              \\[
              TS = TP_{orig} \\times TP_{rep}
              \\]
              
              \\[
              FS = FN_{orig} \\times FN_{rep}
              \\]
              
              \\[
              FF = FN_{orig} \\times TP_{rep}
              \\]
              
              \\[
              TF = TP_{orig} \\times FN_{rep}
              \\]
              
              <br><br>
              
              <b>When the true effect is null (&theta; = 0):</b>
              <br><br>
              
              \\[
              TS = TN_{orig} \\times TN_{rep}
              \\]
              
              \\[
              FS = FP_{orig} \\times FP_{rep}
              \\]
              
              \\[
              FF = FP_{orig} \\times TN_{rep}
              \\]
              
              \\[
              TF = TN_{orig} \\times FP_{rep}
              \\]
              
              <br><br>
              
              All outcome proportions are normalized to sum to one before visualization.
              "))
          )
          
          
        )
      )
    )
    
    
    
  )
  
  
)

server <- function(input, output, session) {
  # ================================
  # Enable/Disable effect size slider
  # ================================
  observe({
    if (input$tabs == "null") {
      shinyjs::disable("effect_size")
    } else {
      shinyjs::enable("effect_size")
    }
  })
  
  # =========================================================
  # OPTIONAL STOPPING FUNCTION (COPY FROM SECOND APP)
  # Bias affects ORIGINAL study only
  # =========================================================
  qrp_optional_stopping <- function(theta, alpha, q,
                                    n0 = 30,
                                    max_mult = 2,
                                    step = 5,
                                    n_sim = 2000) {
    
    if (q == 0) {
      return(list(pi1_os = NA, pi0_os = NA))
    }
    
    n_max <- n0 + q * (n0 * max_mult - n0)
    n_seq <- seq(n0, floor(n_max), by = step)
    
    zcrit <- qnorm(1 - alpha)
    
    simulate_one <- function(mu) {
      for (ni in n_seq) {
        z <- rnorm(1, mean = mu * sqrt(ni / 2), sd = 1)
        if (z > zcrit) return(TRUE)
      }
      FALSE
    }
    
    pi1_os <- mean(replicate(n_sim, simulate_one(theta)))
    pi0_os <- mean(replicate(n_sim, simulate_one(0)))
    
    list(pi1_os = pi1_os, pi0_os = pi0_os)
  }
  
  
  make_shapes <- function(b, alpha, orig_n, rep_n, rep_num) {
    
    offset <- -0.08
    right_offset <- -0.15
    
    # ============================================================
    # 1. ORIGINAL STUDY — FP_orig / TN_orig  (CORRECT)
    # ============================================================
    
    # Nominal type-I error
    alpha0 <- alpha
    
    # Nominal false-positive rate
    pi0_nominal <- alpha0
    
    # Optional stopping (bias affects ORIGINAL only)
    if (b == 0) {
      FP_orig <- pi0_nominal
    } else {
      os <- qrp_optional_stopping(theta = 0, alpha = alpha0, q = b, n0 = orig_n)
      FP_orig <- (1 - b) * pi0_nominal + b * os$pi0_os
    }
    
    TN_orig <- 1 - FP_orig
    
    
    # ============================================================
    # 2. REPLICATION STUDY — FP_rep / TN_rep  (CLEAN)
    # ============================================================
    
    alpha_meta <- 0.05
    
    FP_rep <- alpha_meta
    TN_rep <- 1 - alpha_meta
    
    
    
    # ============================================================
    # 3. 2×2 JOINT LOGIC FOR OUTCOME DECOMPOSITION
    # ============================================================
    
    # Your conceptual logic:
    TS <- TN_orig * TN_rep          # true success
    FS <- FP_orig * FP_rep          # false success
    FF <- FP_orig * TN_rep          # false failure
    TF <- TN_orig * FP_rep          # true failure
    
    # normalize for safety
    total <- TS + FS + FF + TF
    TS <- TS / total
    FS <- FS / total
    FF <- FF / total
    TF <- TF / total
    
    # success/failure bar
    success <- TS + FS
    failure <- FF + TF
    
    true_success  <- TS
    false_success <- FS
    false_failure <- FF
    true_failure  <- TF
    
    
    # ============================================================
    # 4. BUILD PLOTLY SHAPES / TRACES (UNCHANGED)
    # ============================================================
    
    shapes <- list()
    annotations <- list()
    hover_traces <- list()
    
    COLORS <- list(
      "False positives" = "#f5a167",
      "True negatives"  = "#70C170",
      "Successes"       = "#52bf68",
      "Failures"        = "#E5734E",
      "True successes"  = "#5BBBE4",
      "False successes" = "#33A852",
      "False failures"  = "#F7DD6C",
      "True failures"   = "#45629E"
    )
    
    
    # ==========================
    # ORIGINAL STUDIES BAR
    # ==========================
    y0 <- 0.58 + offset
    
    for (pair in list(
      list("False positives", FP_orig, "Proportion of studies rejecting the null hypothesis"),
      list("True negatives",  TN_orig, "Proportion of studies not rejecting the null hypothesis")
    )) {
      
      label <- pair[[1]]; val <- pair[[2]]; tooltip <- pair[[3]]
      y1 <- y0 + val * 0.30
      
      # shape
      shapes <- append(shapes, list(list(
        type="rect", x0=0.05, x1=0.15, y0=y0, y1=y1,
        fillcolor=COLORS[[label]], line=list(color = COLORS[[label]])
      )))
      
      # hover trace
      hover_traces <- append(hover_traces, list(list(
        x=c(0.05, 0.15, 0.15, 0.05),
        y=c(y0, y0, y1, y1),
        type="scatter", mode="none",
        hoverinfo="text",
        text=paste0("<b>", label, "</b><br>", tooltip,
                    "<br>Proportion: ", sprintf("%.1f%%", val*100)),
        fill="toself", fillcolor="rgba(0,0,0,0)",
        hoveron="fills", showlegend=FALSE
      )))
      
      # annotation
      annotations <- append(annotations, list(list(
        x=0.10, y=(y0+y1)/2, text=label,
        showarrow=FALSE, font=list(color="white", size=12)
      )))
      
      y0 <- y1
    }
    
    
    # ==========================
    # REPLICATION STUDIES BAR
    # ==========================
    y0 <- 0.58 + offset
    
    for (pair in list(
      list("False positives", FP_rep, "Proportion of studies showing support for the presence of an effect"),
      list("True negatives",  TN_rep, "Proportion of studies showing support for the absence of an effect")
    )) {
      
      label <- pair[[1]]; val <- pair[[2]]; tooltip <- pair[[3]]
      y1 <- y0 + val * 0.30
      
      shapes <- append(shapes, list(list(
        type="rect", x0=0.25, x1=0.35, y0=y0, y1=y1,
        fillcolor=COLORS[[label]], line=list(color = COLORS[[label]])
      )))
      
      hover_traces <- append(hover_traces, list(list(
        x=c(0.25, 0.35, 0.35, 0.25), y=c(y0, y0, y1, y1),
        type="scatter", mode="none", hoverinfo="text",
        text=paste0("<b>", label, "</b><br>", tooltip,
                    "<br>Proportion: ", sprintf("%.1f%%", val*100)),
        fill="toself", fillcolor="rgba(0,0,0,0)",
        hoveron="fills", showlegend=FALSE
      )))
      
      annotations <- append(annotations, list(list(
        x=0.30, y=(y0+y1)/2, text=label,
        showarrow=FALSE, font=list(color="white", size=12)
      )))
      
      y0 <- y1
    }
    
    
    # ==========================
    # SUCCESS / FAILURE BAR
    # ==========================
    y0 <- 0.07 + offset
    bar_x <- 0.13; bar_width <- 0.12
    
    for (pair in list(
      list("Failures", failure, "Proportion of original and replication studies showing inconsistent results"),
      list("Successes", success, "Proportion of studies showing consistent results")
    )) {
      label <- pair[[1]]; val <- pair[[2]]; tooltip <- pair[[3]]
      y1 <- y0 + val * 0.35
      
      shapes <- append(shapes, list(list(
        type="rect",
        x0=bar_x, x1=bar_x+bar_width, y0=y0, y1=y1,
        fillcolor=COLORS[[label]], line=list(color = COLORS[[label]])
      )))
      
      hover_traces <- append(hover_traces, list(list(
        x=c(bar_x, bar_x+bar_width, bar_x+bar_width, bar_x),
        y=c(y0,y0,y1,y1),
        type="scatter", mode="none", hoverinfo="text",
        text=paste0("<b>", label, "</b><br>", tooltip,
                    "<br>Proportion: ", sprintf("%.1f%%", val*100)),
        fill="toself", fillcolor="rgba(0,0,0,0)",
        hoveron="fills", showlegend=FALSE
      )))
      
      annotations <- append(annotations, list(list(
        x=bar_x+bar_width/2, y=(y0+y1)/2, text=label,
        showarrow=FALSE, font=list(color="white", size=12)
      )))
      
      y0 <- y1
    }
    
    
    # ==========================
    # ANATOMY OF OUTCOMES BAR
    # ==========================
    y0 <- 0.04 + right_offset
    
    for (pair in list(
      list("True failures",  true_failure,
           "The results are inconsistent. Only the original study indicates the effect is null."),
      list("False failures", false_failure,
           "The results are inconsistent. Only the replication correctly indicates the effect is null."),
      list("False successes", false_success,
           "Although the results are consistent, both incorrectly indicate the effect deviates from null."),
      list("True successes", true_success,
           "Both the original and replication results consistently indicate the effect is null.")
    )) {
      
      label <- pair[[1]]; val <- pair[[2]]; tooltip <- pair[[3]]
      y1 <- y0 + val * 0.75
      
      shapes <- append(shapes, list(list(
        type="rect", x0=0.6, x1=0.8, y0=y0, y1=y1,
        fillcolor=COLORS[[label]], line=list(color = COLORS[[label]])
      )))
      
      hover_traces <- append(hover_traces, list(list(
        x=c(0.6,0.9,0.9,0.6),
        y=c(y0,y0,y1,y1),
        type="scatter", mode="none", hoverinfo="text",
        text=paste0("<b>", label, "</b><br>", tooltip,
                    "<br>Proportion: ", sprintf("%.1f%%", val*100)),
        fill="toself", fillcolor="rgba(0,0,0,0)",
        hoveron="fills", showlegend=FALSE
      )))
      
      annotations <- append(annotations, list(list(
        x=0.7, y=(y0+y1)/2, text=label,
        showarrow=FALSE, font=list(color="white", size=11)
      )))
      
      y0 <- y1
    }
    
    
    # return EVERYTHING
    list(shapes = shapes,
         annotations = annotations,
         hover_traces = hover_traces)
  }
  
  
  make_shapes2 <- function(b, alpha, orig_n, rep_n, rep_num) {
    
    offset <- -0.08
    right_offset <- -0.15
    
    # ============================================================
    # 1. ORIGINAL STUDY — TP_orig / FN_orig (CORRECT)
    # ============================================================
    
    theta <- input$effect_size
    alpha0 <- alpha
    
    zcrit <- qnorm(1 - alpha0)
    power_nominal <- 1 - pnorm(zcrit - theta * sqrt(orig_n / 2))
    
    if (b == 0) {
      TP_orig <- power_nominal
    } else {
      os <- qrp_optional_stopping(theta = theta, alpha = alpha0, q = b, n0 = orig_n)
      TP_orig <- (1 - b) * power_nominal + b * os$pi1_os
    }
    
    FN_orig <- 1 - TP_orig
    
    # ============================================================
    # 2. REPLICATION STUDY — TP_rep / FN_rep (CORRECT)
    # ============================================================
    
    alpha_meta <- 0.05
    
    v <- (2 / rep_n) + (theta^2 / (4 * rep_n))
    tau2 <- 0.30 * v / (1 - 0.30)   # I² fixed at 0.30
    v_dot <- (v + tau2) / rep_num
    
    lambda <- theta / sqrt(v_dot)
    zcrit_meta <- qnorm(1 - alpha_meta / 2)
    
    TP_rep <-
      1 - pnorm(zcrit_meta - lambda) +
      pnorm(-zcrit_meta - lambda)
    
    FN_rep <- 1 - TP_rep
    
    
    # ============================================================
    # 3. TRUE 2×2 JOINT OUTCOME CONSTRUCTION
    # ============================================================
    
    TS_raw <- TP_orig * TP_rep    # true success
    FS_raw <- FN_orig * FN_rep    # false success
    FF_raw <- FN_orig * TP_rep    # false failure
    TF_raw <- TP_orig * FN_rep    # true failure
    
    eps <- 1e-6
    TS_raw <- max(TS_raw, eps)
    FS_raw <- max(FS_raw, eps)
    FF_raw <- max(FF_raw, eps)
    TF_raw <- max(TF_raw, eps)
    
    total <- TS_raw + FS_raw + FF_raw + TF_raw
    TS <- TS_raw / total
    FS <- FS_raw / total
    FF <- FF_raw / total
    TF <- TF_raw / total
    
    
    
    
    
    true_success  <- TS
    false_success <- FS
    false_failure <- FF
    true_failure  <- TF
    
    success <- TS + FS
    failure <- TF + FF
    
    
    # ============================================================
    # 5. PLOTLY SHAPES / HOVERS / ANNOTATIONS  (UNCHANGED)
    # ============================================================
    
    shapes <- list()
    annotations <- list()
    hover_traces <- list()
    
    COLORS <- list(
      "False negatives" = "#f5a167",
      "True positives"  = "#70C170",
      "Successes"       = "#52bf68",
      "Failures"        = "#E5734E",
      "True successes"  = "#5BBBE4",
      "False successes" = "#33A852",
      "False failures"  = "#F7DD6C",
      "True failures"   = "#45629E"
    )
    
    
    # ORIGINAL STUDY BAR
    y0 <- 0.58 + offset
    for (pair in list(
      list("False negatives", FN_orig, "Proportion of studies not detecting the effect"),
      list("True positives",  TP_orig, "Proportion of studies detecting the effect")
    )) {
      
      label <- pair[[1]]; val <- pair[[2]]; tooltip <- pair[[3]]
      y1 <- y0 + val * 0.30
      
      shapes <- append(shapes, list(list(
        type="rect", x0=0.05, x1=0.15, y0=y0, y1=y1,
        fillcolor=COLORS[[label]],line=list(color = COLORS[[label]])
      )))
      
      hover_traces <- append(hover_traces, list(list(
        x=c(0.05,0.15,0.15,0.05),
        y=c(y0,y0,y1,y1),
        type="scatter", mode="none", hoverinfo="text",
        text=paste0("<b>",label,"</b><br>",tooltip,
                    "<br>Proportion: ",sprintf("%.1f%%", val*100)),
        fill="toself", fillcolor="rgba(0,0,0,0)",
        hoveron="fills", showlegend=FALSE
      )))
      
      annotations <- append(annotations, list(list(
        x=0.10, y=(y0+y1)/2, text=label,
        showarrow=FALSE, font=list(color="white",size=12)
      )))
      
      y0 <- y1
    }
    
    
    # REPLICATION STUDY BAR
    y0 <- 0.58 + offset
    for (pair in list(
      list("False negatives", FN_rep, "Proportion of replication studies missing the effect"),
      list("True positives",  TP_rep, "Proportion of replication studies detecting the effect")
    )) {
      
      label <- pair[[1]]; val <- pair[[2]]; tooltip <- pair[[3]]
      y1 <- y0 + val * 0.30
      
      shapes <- append(shapes, list(list(
        type="rect", x0=0.25, x1=0.35, y0=y0, y1=y1,
        fillcolor=COLORS[[label]], line=list(color = COLORS[[label]])
      )))
      
      hover_traces <- append(hover_traces, list(list(
        x=c(0.25,0.35,0.35,0.25),
        y=c(y0,y0,y1,y1),
        type="scatter", mode="none", hoverinfo="text",
        text=paste0("<b>",label,"</b><br>",tooltip,
                    "<br>Proportion: ",sprintf("%.1f%%", val*100)),
        fill="toself", fillcolor="rgba(0,0,0,0)",
        hoveron="fills", showlegend=FALSE
      )))
      
      annotations <- append(annotations, list(list(
        x=0.30, y=(y0+y1)/2, text=label,
        showarrow=FALSE, font=list(color="white",size=12)
      )))
      
      y0 <- y1
    }
    
    
    # SUCCESS / FAILURE BAR
    y0 <- 0.07 + offset
    bar_x <- 0.13; bar_width <- 0.12
    
    for (pair in list(
      list("Failures", failure, "Proportion of original and replication studies showing inconsistent results"),
      list("Successes", success, "Proportion of studies showing consistent results")
    )) {
      
      label <- pair[[1]]; val <- pair[[2]]; tooltip <- pair[[3]]
      y1 <- y0 + val * 0.35
      
      shapes <- append(shapes, list(list(
        type="rect", x0=bar_x, x1=bar_x+bar_width, y0=y0, y1=y1,
        fillcolor=COLORS[[label]], line=list(color = COLORS[[label]])
      )))
      
      hover_traces <- append(hover_traces, list(list(
        x=c(bar_x,bar_x+bar_width,bar_x+bar_width,bar_x),
        y=c(y0,y0,y1,y1),
        type="scatter", mode="none", hoverinfo="text",
        text=paste0("<b>",label,"</b><br>",tooltip,
                    "<br>Proportion: ",sprintf("%.1f%%", val*100)),
        fill="toself", fillcolor="rgba(0,0,0,0)",
        hoveron="fills", showlegend=FALSE
      )))
      
      annotations <- append(annotations, list(list(
        x=bar_x+bar_width/2, y=(y0+y1)/2, text=label,
        showarrow=FALSE, font=list(color="white",size=12)
      )))
      
      y0 <- y1
    }
    
    
    # TS / FS / FF / TF BAR
    y0 <- 0.04 + right_offset
    
    for (pair in list(
      list("True failures",  true_failure,  "Original detects effect, replication misses it"),
      list("False failures", false_failure, "Original misses effect, replication detects it"),
      list("False successes", false_success, "Both fail to detect the true effect"),
      list("True successes", true_success,  "Both correctly detect the true effect")
    )) {
      
      label <- pair[[1]]; val <- pair[[2]]; tooltip <- pair[[3]]
      y1 <- y0 + val * 0.75
      
      shapes <- append(shapes, list(list(
        type="rect", x0=0.6, x1=0.8, y0=y0, y1=y1,
        fillcolor=COLORS[[label]], line=list(color = COLORS[[label]])
      )))
      
      hover_traces <- append(hover_traces, list(list(
        x=c(0.6,0.9,0.9,0.6),
        y=c(y0,y0,y1,y1),
        type="scatter", mode="none", hoverinfo="text",
        text=paste0("<b>",label,"</b><br>",tooltip,
                    "<br>Proportion: ",sprintf("%.1f%%", val*100)),
        fill="toself", fillcolor="rgba(0,0,0,0)",
        hoveron="fills", showlegend=FALSE
      )))
      
      annotations <- append(annotations, list(list(
        x=0.70, y=(y0+y1)/2, text=label,
        showarrow=FALSE, font=list(color="white",size=11)
      )))
      
      y0 <- y1
    }
    
    list(shapes=shapes, annotations=annotations, hover_traces=hover_traces)
  }
  
  
  
  
  
  
  
  output$replicationPlot <- renderPlotly({
    elements <- make_shapes(input$bias, input$alpha, input$orig_n, input$rep_n, input$rep_num)
    
    outlines <- list(
      list(type="rect", x0=0.04, x1=0.16, y0=0.49, y1=0.81,
           line=list(color="black", width=1.2, dash="dot"), fillcolor="rgba(0,0,0,0)"),
      list(type="rect", x0=0.24, x1=0.36, y0=0.49, y1=0.81,
           line=list(color="black", width=1.2, dash="dot"), fillcolor="rgba(0,0,0,0)"),
      list(type="rect", x0=0.115, x1=0.265, y0=-0.02, y1=0.35,
           line=list(color="black", width=1.2, dash="dot"), fillcolor="rgba(0,0,0,0)"),
      list(type="rect", x0=0.59, x1=0.81, y0=-0.15, y1=0.67,
           line=list(color="black", width=1.2, dash="dot"), fillcolor="rgba(0,0,0,0)")
    )
    
    connectors <- list(
      list(type="line", x0=0.10, y0=0.49, x1=0.10, y1=0.43,
           line=list(color="black", width=1.5, dash="dot")),
      list(type="line", x0=0.30, y0=0.49, x1=0.30, y1=0.43,
           line=list(color="black", width=1.5, dash="dot")),
      list(type="line", x0=0.10, y0=0.43, x1=0.30, y1=0.43,
           line=list(color="black", width=1.5, dash="dot")),
      list(type="line", x0=0.19, y0=0.43, x1=0.19, y1=0.38,
           line=list(color="black", width=1.5, dash="dot")),
      list(type="path", path="M 0.1875 0.379 L 0.19 0.363 L 0.1925 0.379 Z",
           fillcolor="black", line=list(width=0)),
      list(type="line", x0=0.265, y0=0.35, x1=0.59, y1=0.67,
           line=list(color="black", width=1.5, dash="dot")),
      list(type="line", x0=0.265, y0=-0.02, x1=0.59, y1=-0.15,
           line=list(color="black", width=1.5, dash="dot"))
    )
    
    arrows <- list(
      list(type="line", x0=0.10, y0=0.945, x1=0.10, y1=0.84,
           line=list(color="black", width=1.2)),
      list(type="line", x0=0.30, y0=0.945, x1=0.30, y1=0.84,
           line=list(color="black", width=1.2)),
      list(type="path", path="M 0.097 0.84 L 0.10 0.830 L 0.103 0.84 Z",
           fillcolor="black", line=list(width=0)),
      list(type="path", path="M 0.297 0.84 L 0.30 0.830 L 0.303 0.84 Z",
           fillcolor="black", line=list(width=0))
    )
    
    annotations <- append(elements$annotations, list(
      
      # --- Original Study ---
      list(
        x = 0.1, y = 0.965,
        text = "<b>Original Study</b>",
        showarrow = FALSE,
        bordercolor = "black", borderwidth = 1, borderpad = 4,
        hovertext = "Two-group between-subjects design",
        hoverlabel = list(bgcolor = "white", font = list(color = "black")),
        hoverarrow = FALSE
      ),
      
      # --- Replications ---
      list(
        x = 0.3, y = 0.965,
        text = "<b>Replications</b>",
        showarrow = FALSE,
        bordercolor = "black", borderwidth = 1, borderpad = 4,
        hovertext = "Multiple identical replication studies",
        hoverlabel = list(bgcolor = "white", font = list(color = "black")),
        hoverarrow = FALSE
      ),
      
      # --- Other labels ---
      list(x=0.04, y=0.91, text="<i>Null hypothesis testing</i>",
           showarrow=FALSE, font=list(size=11), align="right"),
      
      list(x=0.23, y=0.91, text="<i>Meta-analytic Bayesian<br>hypothesis testing</i>",
           showarrow=FALSE, font=list(size=11), align="right"),
      
      list(x=0.07, y=0.18, text="<b>Replication<br>Outcome<br>Classification</b>", showarrow=FALSE),
      
      list(
        x=0.7, y=0.72,
        text="<b>Anatomy of Outcomes</b>",
        showarrow=FALSE,
        font = list(size = 12)
      ),
      
      list(
        x=0.68, y=0.97,
        text="<b>When true effect is null</b>",
        showarrow=FALSE,
        font = list(size = 15)
      )
    ))
    
    
    p <- plot_ly()
    
    # Correct: add traces with do.call
    for (trace in elements$hover_traces) {
      p <- do.call(add_trace, c(list(p), trace))
    }
    
    p %>% layout(
      xaxis = list(showticklabels = FALSE, showgrid = FALSE,
                   zeroline = FALSE, showline = FALSE, range = c(0, 1)),
      yaxis = list(showticklabels = FALSE, showgrid = FALSE,
                   zeroline = FALSE, showline = FALSE, range = c(-0.2, 1)),
      shapes = c(elements$shapes, outlines, connectors, arrows),
      annotations = annotations,
      width = 1200, height = 800,
      plot_bgcolor = "white", paper_bgcolor = "white"
    )
  })
  
  output$replicationPlot2 <- renderPlotly({
    elements <- make_shapes2(input$bias, input$alpha, input$orig_n, input$rep_n, input$rep_num)
    
    # Frame outlines and connection lines
    outlines <- list(
      list(type="rect", x0=0.04, x1=0.16, y0=0.49, y1=0.81,
           line=list(color="black", width=1.2, dash="dot"), fillcolor="rgba(0,0,0,0)"),
      list(type="rect", x0=0.24, x1=0.36, y0=0.49, y1=0.81,
           line=list(color="black", width=1.2, dash="dot"), fillcolor="rgba(0,0,0,0)"),
      list(type="rect", x0=0.115, x1=0.265, y0=-0.02, y1=0.35,
           line=list(color="black", width=1.2, dash="dot"), fillcolor="rgba(0,0,0,0)"),
      list(type="rect", x0=0.59, x1=0.81, y0=-0.15, y1=0.67,
           line=list(color="black", width=1.2, dash="dot"), fillcolor="rgba(0,0,0,0)")
    )
    
    connectors <- list(
      list(type="line", x0=0.10, y0=0.49, x1=0.10, y1=0.43,
           line=list(color="black", width=1.5, dash="dot")),
      list(type="line", x0=0.30, y0=0.49, x1=0.30, y1=0.43,
           line=list(color="black", width=1.5, dash="dot")),
      list(type="line", x0=0.10, y0=0.43, x1=0.30, y1=0.43,
           line=list(color="black", width=1.5, dash="dot")),
      list(type="line", x0=0.19, y0=0.43, x1=0.19, y1=0.38,
           line=list(color="black", width=1.5, dash="dot")),
      list(type="path", path="M 0.1875 0.379 L 0.19 0.363 L 0.1925 0.379 Z",
           fillcolor="black", line=list(width=0)),
      list(type="line", x0=0.265, y0=0.35, x1=0.59, y1=0.67,
           line=list(color="black", width=1.5, dash="dot")),
      list(type="line", x0=0.265, y0=-0.02, x1=0.59, y1=-0.15,
           line=list(color="black", width=1.5, dash="dot"))
    )
    
    arrows <- list(
      list(type="line", x0=0.10, y0=0.945, x1=0.10, y1=0.84,
           line=list(color="black", width=1.2)),
      list(type="line", x0=0.30, y0=0.945, x1=0.30, y1=0.84,
           line=list(color="black", width=1.2)),
      list(type="path", path="M 0.097 0.84 L 0.10 0.830 L 0.103 0.84 Z",
           fillcolor="black", line=list(width=0)),
      list(type="path", path="M 0.297 0.84 L 0.30 0.830 L 0.303 0.84 Z",
           fillcolor="black", line=list(width=0))
    )
    
    annotations <- append(elements$annotations, list(
      list(x=0.1, y=0.965, text="<b>Original Study</b>", showarrow=FALSE,
           bordercolor="black", borderwidth=1, borderpad=4,hovertext = "Two-group between-subjects design",
           hoverlabel = list(bgcolor = "white", font = list(color = "black")),
           hoverarrow = FALSE),
      list(x=0.3, y=0.965, text="<b>Replications</b>", showarrow=FALSE,
           bordercolor="black", borderwidth=1, borderpad=4,hovertext = "Multiple idential replication studies",
           hoverlabel = list(bgcolor = "white", font = list(color = "black")),
           hoverarrow = FALSE),
      
      list(x=0.04, y=0.91, text="<i>Null hypothesis testing</i>",
           showarrow=FALSE, font=list(size=11), align="right"),
      
      list(x=0.23, y=0.91, text="<i>Meta-analytic Bayesian<br>hypothesis testing</i>",
           showarrow=FALSE, font=list(size=11), align="right"),
      list(
        x=0.7, y=0.72,
        text="<b>Anatomy of Outcomes</b>",
        showarrow=FALSE,
        font = list(size = 12)
      ),
      list(x=0.07, y=0.18, text="<b>Replication<br>Outcome<br>Classification</b>", showarrow=FALSE),
      list(x=0.68, y=0.97,
           text="<b>When true effect deviates from null</b>",
           showarrow=FALSE,
           font = list(size = 15))
    ))
    
    p <- plot_ly()
    for (trace in elements$hover_traces) {
      p <- do.call(add_trace, c(list(p), trace))
    }
    
    p %>% layout(
      xaxis = list(showticklabels = FALSE, showgrid = FALSE,
                   zeroline = FALSE, showline = FALSE, range = c(0, 1)),
      yaxis = list(showticklabels = FALSE, showgrid = FALSE,
                   zeroline = FALSE, showline = FALSE, range = c(-0.2, 1)),
      shapes = c(elements$shapes, outlines, connectors, arrows),
      annotations = annotations,
      width = 1200, height = 800,
      plot_bgcolor = "white", paper_bgcolor = "white"
    )
  })
  
  
  
  
  
  
  
  output$categoryDiagram <- renderPlotly({
    
    # ================
    # Define function
    # ================
    make_diagram <- function(x_offset = 0, mode = "theta0") {
      # -------- Label setup --------
      if (mode == "theta0") {
        # Left diagram: θ = 0  (Effect not exist)
        label_data <- data.frame(
          row = rep(1:4, each = 3),
          col = rep(1:3, times = 4),
          text = c(
            "True negative", "True negative", "True success",
            "False positive", "False positive", "False success",
            "False positive", "True negative", "False failure",
            "True negative", "False positive", "True failure"
          ),
          fill = c(
            "#70C170", "#70C170", "#5BBBE4",
            "#f5a167", "#f5a167", "#33A852",
            "#f5a167", "#70C170", "#F7DD6C",
            "#70C170", "#f5a167", "#45629E"
          )
        )
      } else if (mode == "thetaExists") {
        # Right diagram: θ ≠ 0  (Effect exists)
        label_data <- data.frame(
          row = rep(1:4, each = 3),
          col = rep(1:3, times = 4),
          text = c(
            "True positive", "True positive", "True success",
            "False negative", "False negative", "False success",
            "False negative", "True positive", "False failure",
            "True positive", "False negative", "True failure"
          ),
          fill = c(
            "#70C170", "#70C170", "#5BBBE4",
            "#f5a167", "#f5a167", "#33A852",
            "#f5a167", "#70C170", "#F7DD6C",
            "#70C170", "#f5a167", "#45629E"
          )
        )
      }
      
      # -------- Layout coordinates --------
      x_pos <- c(0.15, 0.5, 0.85) + x_offset
      y_pos <- c(0.85, 0.62, 0.39, 0.16)
      
      label_data$x <- x_pos[label_data$col]
      label_data$y <- y_pos[label_data$row]
      
      shapes <- list()
      annotations <- list()
      
      # -------- Draw boxes + labels --------
      for (i in 1:nrow(label_data)) {
        shapes[[i]] <- list(
          type = "rect",
          x0 = label_data$x[i] - 0.12, x1 = label_data$x[i] + 0.12,
          y0 = label_data$y[i] - 0.08, y1 = label_data$y[i] + 0.08,
          fillcolor = label_data$fill[i],
          line = list(color = "black")
        )
        annotations[[i]] <- list(
          x = label_data$x[i],
          y = label_data$y[i],
          text = label_data$text[i],
          showarrow = FALSE,
          font = list(color = "white", size = 12)
        )
      }
      
      # -------- Column headers --------
      headers <- data.frame(
        x = c(0.15, 0.5, 0.85) + x_offset,
        y = 1.05,
        text = c("Original", "Replication", "Outcome")
      )
      
      for (i in 1:nrow(headers)) {
        annotations[[length(annotations) + 1]] <- list(
          x = headers$x[i],
          y = headers$y[i],
          text = paste0("<b>", headers$text[i], "</b>"),
          showarrow = FALSE,
          font = list(size = 14),
          bordercolor = "black",
          borderwidth = 1,
          borderpad = 4
        )
      }
      
      # -------- Add × and : symbols --------
      for (i in 1:4) {
        annotations[[length(annotations) + 1]] <- list(
          x = 0.325 + x_offset,
          y = y_pos[i],
          text = "<b>×</b>",
          showarrow = FALSE,
          font = list(size = 16, color = "black")
        )
        annotations[[length(annotations) + 1]] <- list(
          x = 0.675 + x_offset,
          y = y_pos[i],
          text = "<b>:</b>",
          showarrow = FALSE,
          font = list(size = 16, color = "black")
        )
      }
      
      list(shapes = shapes, annotations = annotations)
    }
    
    # =============================
    # Draw both diagrams (side-by-side)
    # =============================
    left  <- make_diagram(x_offset = 0, mode = "theta0")         # θ = 0 (no effect)
    right <- make_diagram(x_offset = 1.3, mode = "thetaExists")  # θ ≠ 0 (effect exists)
    
    # =============================
    # Add θ labels above each diagram
    # =============================
    theta_headers <- data.frame(
      x = c(0.5, 1.8),
      y = 1.18,
      text = c("When effect is null", "When effect is non-null")
    )
    
    for (i in 1:nrow(theta_headers)) {
      left$annotations[[length(left$annotations) + 1]] <- list(
        x = theta_headers$x[i],
        y = theta_headers$y[i],
        text = paste0("<b>", theta_headers$text[i], "</b>"),
        showarrow = FALSE,
        font = list(size = 14),
        align = "center"
      )
    }
    
    # =============================
    # Combine and render both
    # =============================
    plot_ly() %>%
      layout(
        shapes = c(left$shapes, right$shapes),
        annotations = c(left$annotations, right$annotations),
        xaxis = list(showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE, range = c(0, 3)),
        yaxis = list(showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE, range = c(0, 1.3)),
        width = 1400,
        height = 450,
        plot_bgcolor = "white",
        paper_bgcolor = "white"
      )
  })
  
  output$png_display <- renderImage({
    
    rep_num <- as.numeric(input$png_rep_num)
    rep_n   <- as.numeric(input$png_rep_n)
    es      <- input$png_es   # "0.2" or "0.5" as a string
    
    file_name <- sprintf(
      "ROC_plot_%d_repN%d_ES%s_onesided.png",
      rep_num,
      rep_n,
      es
    )
    
    list(
      src = file.path("roc4shiny", file_name),
      contentType = "image/png",
      width = 1000,
      alt = file_name
    )
    
  }, deleteFile = FALSE)
  
}

shinyApp(ui, server)


#############
#The following code generate the category diagram
# library(shiny)
# library(plotly)
# 
# ui <- fluidPage(
#   plotlyOutput("categoryDiagram", height = "500px")
# )
# 
# server <- function(input, output, session) {
#   
#   output$categoryDiagram <- renderPlotly({
#     
#     make_diagram <- function(x_offset = 0, mode = "theta0") {
#       
#       if (mode == "theta0") {
#         label_data <- data.frame(
#           row = rep(1:4, each = 3),
#           col = rep(1:3, times = 4),
#           text = c(
#             "True negative", "True negative", "True success",
#             "False positive", "False positive", "False success",
#             "False positive", "True negative", "False failure",
#             "True negative", "False positive", "True failure"
#           ),
#           fill = c(
#             "#70C170", "#70C170", "#5BBBE4",
#             "#f5a167", "#f5a167", "#33A852",
#             "#f5a167", "#70C170", "#F7DD6C",
#             "#70C170", "#f5a167", "#45629E"
#           ),
#           stringsAsFactors = FALSE
#         )
#       } else if (mode == "thetaExists") {
#         label_data <- data.frame(
#           row = rep(1:4, each = 3),
#           col = rep(1:3, times = 4),
#           text = c(
#             "True positive", "True positive", "True success",
#             "False negative", "False negative", "False success",
#             "False negative", "True positive", "False failure",
#             "True positive", "False negative", "True failure"
#           ),
#           fill = c(
#             "#70C170", "#70C170", "#5BBBE4",
#             "#f5a167", "#f5a167", "#33A852",
#             "#f5a167", "#70C170", "#F7DD6C",
#             "#70C170", "#f5a167", "#45629E"
#           ),
#           stringsAsFactors = FALSE
#         )
#       }
#       
#       x_pos <- c(0.15, 0.5, 0.85) + x_offset
#       y_pos <- c(0.85, 0.62, 0.39, 0.16)
#       
#       label_data$x <- x_pos[label_data$col]
#       label_data$y <- y_pos[label_data$row]
#       
#       shapes <- list()
#       annotations <- list()
#       
#       for (i in seq_len(nrow(label_data))) {
#         shapes[[i]] <- list(
#           type = "rect",
#           x0 = label_data$x[i] - 0.12,
#           x1 = label_data$x[i] + 0.12,
#           y0 = label_data$y[i] - 0.08,
#           y1 = label_data$y[i] + 0.08,
#           fillcolor = label_data$fill[i],
#           line = list(color = "black")
#         )
#         
#         annotations[[i]] <- list(
#           x = label_data$x[i],
#           y = label_data$y[i],
#           text = label_data$text[i],
#           showarrow = FALSE,
#           font = list(color = "white", size = 12)
#         )
#       }
#       
#       headers <- data.frame(
#         x = c(0.15, 0.5, 0.85) + x_offset,
#         y = 1.05,
#         text = c("Original study", "Replication study", "Determination"),
#         stringsAsFactors = FALSE
#       )
#       
#       for (i in seq_len(nrow(headers))) {
#         annotations[[length(annotations) + 1]] <- list(
#           x = headers$x[i],
#           y = headers$y[i],
#           text = paste0("<b>", headers$text[i], "</b>"),
#           showarrow = FALSE,
#           font = list(size = 14),
#           bordercolor = "black",
#           borderwidth = 1,
#           borderpad = 4
#         )
#       }
#       
#       for (i in 1:4) {
#         annotations[[length(annotations) + 1]] <- list(
#           x = 0.325 + x_offset,
#           y = y_pos[i],
#           text = "<b>×</b>",
#           showarrow = FALSE,
#           font = list(size = 16, color = "black")
#         )
#         
#         annotations[[length(annotations) + 1]] <- list(
#           x = 0.675 + x_offset,
#           y = y_pos[i],
#           text = "<b>:</b>",
#           showarrow = FALSE,
#           font = list(size = 16, color = "black")
#         )
#       }
#       
#       list(shapes = shapes, annotations = annotations)
#     }
#     
#     left  <- make_diagram(x_offset = 0,   mode = "theta0")
#     right <- make_diagram(x_offset = 1.3, mode = "thetaExists")
#     
#     theta_annotations <- list(
#       list(
#         x = 0.5, y = 1.18,
#         text = "<b><i>θ</i> = 0</b>",
#         showarrow = FALSE,
#         font = list(size = 14)
#       ),
#       list(
#         x = 1.8, y = 1.18,
#         text = "<b><i>θ</i> = 0.2 or 0.5</b>",
#         showarrow = FALSE,
#         font = list(size = 14)
#       )
#     )
#     
#     plot_ly() %>%
#       layout(
#         shapes = c(left$shapes, right$shapes),
#         annotations = c(left$annotations, right$annotations, theta_annotations),
#         xaxis = list(
#           showticklabels = FALSE,
#           showgrid = FALSE,
#           zeroline = FALSE,
#           range = c(0, 2.5)
#         ),
#         yaxis = list(
#           showticklabels = FALSE,
#           showgrid = FALSE,
#           zeroline = FALSE,
#           range = c(0, 1.3)
#         ),
#         width = 1400,
#         height = 500,
#         plot_bgcolor = "white",
#         paper_bgcolor = "white"
#       ) %>%
#       config(
#         toImageButtonOptions = list(
#           format = "png",
#           filename = "categoryDiagram",
#           width = 1400,
#           height = 500,
#           scale = 5
#         )
#       )
#   })
# }
# 
# shinyApp(ui, server)
# 
# 
