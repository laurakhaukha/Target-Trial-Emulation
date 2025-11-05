# Target-Trial-Emulation
Causal Effect of Statin Therapy on Recurrent Cardiovascular Disease Using Target Trial Emulation

# Overview 👩🏽‍💻 
This repository reproduces an observational causal inference study evaluating the effect of statin initiation on recurrent cardiovascular disease (CVD) events using a Target Trial Emulation (TTE) framework and messy, synthetic data. The analysis uses pooled logistic regression to estimate ITT and IPW-based per-protocol effects, and includes bootstrap-based cumulative risk curve estimation.

# Key Methods 🔨🔨
- Target Trial Emulation (Hernán & Robins framework)
- Descriptive baseline covariate analysis
- Intention-to-Treat (ITT) effect estimation via pooled logistic regression
- Per-Protocol (PP) effect estimation using:
  Inverse Probability Weighting (IPW) for adherence
- Weighted pooled logistic regression
- Bootstrap-based cumulative risk estimation

Results (Brief Summary) 📊
- ITT estimate showed no statistically significant effect of baseline statin initiation.
- PP risk curves indicated a short-term protective effect (months 5–18) under sustained adherence.
- Effects weakened toward month 30, suggesting limited long-term impact or confounding.


