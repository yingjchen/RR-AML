# Prediction of selective and synergistic drug combinations for relapsed AML

We develop a systematic combinatorial design strategy that uses machine learning to prioritise the most promising targeted drug combinations for relapsed/refractory AML (**RR-AML**) patients using single-cell transcriptomics and single-agent response profiles measured in primary patient samples. By utilizng the established target-based Normalized Single-cell Enrichment (**t-NSE**) score, we can quantitatively compare the co-inhibition effects of drug combinations among various cell types and prioritize combinations that exhibit high synergy and potency in co-inhibiting AML cells, while showing non-synergistic effects in non-malignant cells. The following figure illustrates the workflow of the drug combination prediction and testing pipeline.
![Workflow](./Figures/Workflow1.png)



# Contact information
For any questions please contact **Yingjia Chen** (yingjia.chen@helsinki.fi)

# Copyright and license
Code copyright *Prediction of selective and synergistic drug combinations for relapsed AML*

License <https://github.com/yingjchen/RR-AML/blob/main/LICENSE>

# Acknowledgements
The original version of the code for developing XGBoost models was forked from https://github.com/IanevskiAleksandr/scComb. Many thanks to Aleksandr Ianevski for making his code available.
