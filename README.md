# Branch Description: experiment-wage-heterogeneity

This research extension leverages an SVAR framework to investigate the heterogeneous effects of **Price Markup Shocks** on the real wage growth of different occupational groups (**White-collar vs. Blue-collar**) in the US labor market.

## 1. Occupational Classification
I classify the labor force into **White-collar** and **Blue-collar** workers based on the CPS Occupational Codes (OCC). The logic is as follows:
- **White-collar (Label: 1)**: Includes management, professional, technical, administrative support, and service occupations.
- **Blue-collar (Label: 0)**: Includes production, maintenance, operators, and manual labor occupations.

## 2. Preliminary Findings
Based on the quarterly data (1987 - 2012), our analysis yields several key insights:

### Case 1: after hp-filter
#### 1. Divergent Immediate Impact at Lag 0 ("Predatory" vs. "Redistributive")
At the moment the shock occurs (**Lag 0**), I observe a striking divergence between the two groups:
* **Blue-collar workers** exhibit a significant **negative correlation** (approx. -0.24). This suggests a unexpected price markups immediately erode the real wage growth of manual labor.
* **White-collar workers** show a **positive correlation** (approx. +0.08). This indicates a "bonus" effect.

#### 2. Heterogeneous Response Dynamics (Reaction Speed vs. Long-term Compensation)
The dynamic path over subsequent quarters reveals a clear trade-off between reaction speed and total recovery:
* **Blue-collar: High Sensitivity, Rapid Adjustment.** While hit hardest initially, blue-collar wages show a sharp rebound by Lag 1. This reflects shorter contract cycles or more frequent wage re-negotiations in response to market volatility.
* **White-collar: High Rigidity, Full Recovery.** White-collar wages experience a delayed "slump" (Lag 1) but achieve the highest positive correlation by **Lag 4** (approx. +0.21).

### Case 2: without hp-filter
#### 1. Synchronized Immediate Impact at Lag 0
* **Blue-collar workers** exhibit a **negative correlation** (approx. -0.11). This suggests a unexpected price markups immediately erode the real wage growth of manual labor.
* **White-collar workers** also show a **negative correlation** (approx. -0.12).

#### 2. Heterogeneous Response Dynamics (Volatility vs. Structural Rigidity)
The dynamic path reveals a sharp contrast in how different occupations response to the shock:
* **Blue-collar: High volatility.** Blue-collar wages show a "zig-zag" pattern with high sensitivity. They rebound rapidly and achieve a significant compensatory spike by Lag 3 (approx. +0.17), reflecting shorter contract cycles.
* **White-collar: High Rigidity** White-collar wages follow a much smoother, lagged path. They experience a prolonged but shallow slump, only returning to a modest positive correlation by Lag 4 (approx. +0.05), confirming the structural wage rigidity of professional roles.