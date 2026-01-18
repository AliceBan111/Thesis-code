# Branch Description: v1-stable-methods
**Date**: 2026-01-09

### Contents of this Branch

#### 1. Implementation of Method 1 & Method 2
* **Method 1**: Nekarda & Ramey (2020)
* **Method 2**: Bils et al. (2018)

#### 2. Data Alignment & Preprocessing Logic
* **Standardized look_back mechanism**: Implemented logic to retrieve data one quarter prior to the official start date.
* **Sample Consistency**: This ensures that the quarterly growth rate for the initial observation (**1964 Q1**) is calculated correctly rather than being marked as `Missing`, maintaining a consistent sample size across different estimation methods for the VAR model.

#### 3. VAR Diagnostics and Visualization
* **Automated Export**: All VAR and IRF plots are automatically organized and saved into their respective directories: `results/method1/` and `results/method2/`.