# Research Log: Stage 1

**Last Updated:** 2026-01-10
**Status:** Completed baseline SVAR identification and preliminary wage heterogeneity experiments.

## Repository Structure

```text
.
├── main.jl                 # Main Entry: Automates data prep and baseline estimation
├── README.md               # This file: Research log and directory guidance
├── .gitignore              # Git Ignore: Excludes large .dat files and temp outputs
│
├── src/                    # Source Code Directory
│   ├── Toolbox/            # [Engine] Low-level algorithms for VAR, QR rotation & statistics
|   ├── models/             # [Core] SVAR logic and sign restriction algorithms
│   ├── data/               # [Core] Macroeconomic data cleaning
│   └── wage_shock_analysis/# [Experimental] CPS micro-data processing and alignment
│
├── data/                   # Data Directory (Tracks structure, excludes .dat files)               
│
└── results/                # Output Directory
    ├── method1/            # Preliminary results from Identification Method 1
    └── method2/            # Analysis results from Identification Method 2
        └── wage_heterogeneity/ # [Exp Branch] White-collar vs. Blue-collar analysis
```

## Branching & Project Milestones

### **`main` Branch**
* **Status**: Stable baseline SVAR framework.
* **Scope**: Contains the primary codebase for macroeconomic variable estimation and standard structural identification.

### **`experiment-wage-heterogeneity` Branch (Stage 1 Finished)**
* **Focus**: Investigation of wage disparities under price markup shocks.
* **Key Achievements**:
    * **Data Integration**: Successfully integrated CPS quarterly micro-data spanning 1987–2012.
    * **Classification**: Implemented a occupational classification system (White-collar vs. Blue-collar) based on Census OCC codes.
* **Key Insight**: 
    * Identified a **"predatory" immediate impact** on blue-collar wages (significant negative correlation at Lag 0).
    * Observed **delayed compensation** patterns in white-collar wages (positive correlation/recovery at later lags).

> **Note**: For Stage 1 specific occupational definitions and the latest **"Wage Correlation Response"** plots, please switch to the research branch via:
> `git checkout experiment-wage-heterogeneity`

## Quick Start
1. **Environment**: Julia 1.9+ required.
2. **Execution**: Run julia main.jl from the project root.
3. **Outputs**: Baseline macro IRF plots are saved under results/.