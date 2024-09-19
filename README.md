# Prediction of selective and synergistic drug combinations for relapsed AML

We develop a systematic combinatorial design strategy that uses machine learning to prioritise the most promising targeted drug combinations for relapsed/refractory AML (**RR-AML**) patients using single-cell transcriptomics and single-agent response profiles measured in primary patient samples. By utilizng the established target-based Normalized Single-cell Enrichment (**t-NSE**) score, we can quantitatively compare the co-inhibition effects of drug combinations among various cell types and prioritize combinations that exhibit high synergy and potency in co-inhibiting AML cells, while showing non-synergistic effects in non-malignant cells. The following figure illustrates the workflow of the drug combination prediction and testing pipeline.

<p align = "center">
    <img src="./Figures/Workflow1.png" alt="Workflow..." style="width:90%; height: 90%" height="90%" width="90%"/>
</p>


# Requirements

To install all dependencies, the version of R should be >= 4.0.5. The required packages and their recommended versions are as follows:
```r
[1] ModelMetrics_1.2.2.2 caret_6.0-92         parallel_4.0.5     
[4] xgboost_1.5.1.1      Seurat_4.1.1          ggplot2_3.3.6       
[7] readr_2.1.2          copykat_1.0.8        HGNChelper_0.8.1    
[10] dplyr_1.0.9          lhs_1.2.0          ParamHelpers_1.14.1                

```



# Contact information
For any questions please contact **Yingjia Chen** (yingjia.chen@helsinki.fi)

# Copyright and license
Code copyright *Prediction of selective and synergistic drug combinations for relapsed AML*

License <https://github.com/yingjchen/RR-AML/blob/main/LICENSE>

# Acknowledgements
The original version of the code for developing XGBoost models was forked from https://github.com/IanevskiAleksandr/scComb. Many thanks to Aleksandr Ianevski for making his code available.