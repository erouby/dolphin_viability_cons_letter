# 🐬 Dolphin Viability Conservation Letter

## 🌊 Overview

This repository contains the complete analytical framework and reproducible code for the research presented in **"Longevity collapse in dolphins: critical seven-year drop in life expectancy reveals accelerating conservation emergency in the Bay of Biscay"** published in Conservation Letters.

Our study demonstrates the first evidence of declining viability in common dolphins (*Delphinus delphis*) in the Bay of Biscay, documenting a dramatic 7-year decline in female life expectancy (from 24 to 17 years) between 1997-2019, corresponding to a 2.4% reduction in population growth rate.

## 🧬 Background

Population Viability Analysis (PVA) using stranding data represents a novel approach for monitoring wide-ranging cetacean populations where traditional Capture-Mark-Recapture methods are unfeasible. This research showcases how demographic data from stranded specimens can provide early warning indicators of population decline before abundance trends become detectable.

## 📊 Key Findings

- **📉 Female longevity declined from 24 to 17 years** between 1997-2019
- **⚠️ 2.4% reduction in population growth rate** over the study period  
- **🚨 Bay of Biscay identified as demographic sink** despite stable abundance estimates
- **⏰ Early warning system** for proactive conservation management

## 🗂️ Repository Structure

```
dolphin_viability_cons_letter/
├── 📁 data/                     # Demographic and stranding data
├── 📁 model/                    # Stan model files (.stan)
│   ├── M1_null.stan            # Null model (no effects)
│   ├── M2_cov.stan             # Covariate model
│   ├── M3_random.stan          # Random effects model
│   ├── M4_trend.stan           # Temporal trend model
│   ├── M5_random_trend.stan    # Random + trend effects
│   ├── M6_cov_random.stan      # Covariate + random effects
│   ├── M7_cov_trend.stan       # Covariate + trend effects
│   └── M8_cov_random_trend.stan # Full model
├── 📁 scripts/                  # Analysis workflow scripts
│   ├── 000_prep.R              # Data preparation
│   ├── 001_analysis.R          # Main demographic analysis
│   ├── 002_model_selection.R   # Model comparison
│   └── 003_loop.R              # Batch processing
├── 📁 source/                   # Custom functions
│   ├── functions_M7.R          # Model 7 helper functions
│   └── get_growth_rates.R      # Population growth calculations
├── 📄 00_run.R                 # Master script
├── 📄 01_convergence.R         # Model convergence checks
├── 📄 01_get_asm_logistic.R    # Age at sexual maturity analysis
├── 📄 02_results.R             # Results compilation
├── 📄 03_show_raw_age_at_death_data.R # Data visualization
└── 📄 README.md               # This file
```

## 🔧 Prerequisites

### Required Software
- **R** (≥ 4.3.1)
- **Stan** (≥ 2.21.8) via RStan
- **Git** (for cloning repository)

### Required R Packages
```r
# Core analysis packages
install.packages(c("rstan", "loo", "ggplot2", "dplyr"))

# Additional packages for specific analyses
install.packages(c("tidyr", "gridExtra", "cowplot"))
```

## 🚀 Getting Started - Step by Step

### Step 1: Clone the Repository 📥
```bash
git clone https://github.com/erouby/dolphin_viability_cons_letter.git
cd dolphin_viability_cons_letter
```

### Step 2: Setup R Environment 🔧
```r
# Open R in the project directory
# Install required packages if not already installed
source("00_run.R")  # This will check and install dependencies
```

### Step 3: Run Complete Analysis 🏃‍♂️
```r
# Execute the master script to reproduce all analyses
source("00_run.R")
```

This will sequentially run:
1. **Data preparation** (`scripts/000_prep.R`)
2. **Model compilation** (Stan models in `model/` directory)
3. **Model selection** (`scripts/002_model_selection.R`)
4. **Results generation** (`02_results.R`)
5. **Figure creation** (`03_show_raw_age_at_death_data.R`)

### Step 4: Check Model Convergence 📈
```r
source("01_convergence.R")
```

### Step 5: Generate Sexual Maturity Analysis 📊
```r
source("01_get_asm_logistic.R")
```

## 📋 Detailed Workflow

### 🔄 Model Selection Process
The analysis compares 8 different demographic models:
1. **M1**: Null model (no temporal effects)
2. **M2**: Covariate effects only  
3. **M3**: Random year effects
4. **M4**: Linear temporal trend
5. **M5**: Random effects + trend
6. **M6**: Covariates + random effects
7. **M7**: Covariates + trend ⭐ **(Best model)**
8. **M8**: Full model (all effects)

Model selection uses Leave-One-Out Information Criterion (LOOIC) and WAIC.

### 📈 Key Outputs
- **Survivorship curves** by year and sex
- **Population growth rate trends** over time
- **Age at sexual maturity** estimates
- **Model convergence diagnostics**
- **Demographic parameter estimates**

## 🔬 Using This Framework with Your Own Data

To adapt this analysis for your own cetacean demographic data:

### Step 1: Data Structure 📋
Ensure your data follows the same structure as in `data/` folder:
- Age-at-death records
- Sex information
- Year of stranding/death
- Reproductive status (for females)
- Body length measurements

### Step 2: Modify Parameters ⚙️
Update the following in your analysis scripts:
- Species-specific biological parameters
- Study period dates
- Geographic region information
- Sample size considerations

### Step 3: Model Adaptation 🔧
The Stan models are flexible and can accommodate:
- Different age ranges
- Various demographic covariates
- Alternative study periods
- Species-specific mortality patterns

## 📚 Associated Publication

**Citation:**
> Rouby, E., Plard, F., Ridoux, V., Mauchamp, A., Dabin, W., Spitz, J., & Authier, M. (2025). 
> Longevity collapse in dolphins: critical seven-year drop in life expectancy reveals accelerating 
> conservation emergency in the Bay of Biscay. *Conservation Letters*, [DOI pending].

## 📊 Data Availability

- **Stranding data**: Available through PelaMSFD R package
- **Processed demographic data**: Included in this repository
- **Raw biological samples**: Contact Observatoire Pelagis for access

## 👥 Authors & Contact

**Lead Author:** Etienne Rouby  
📧 Email: etienne.rouby@colorado.edu  
🏛️ Affiliations: 
- Institute of Arctic and Alpine Research, University of Colorado Boulder
- Woods Hole Oceanographic Institution
- Observatoire Pelagis, La Rochelle Université

**Corresponding Author:** Etienne Rouby

## 🙏 Acknowledgments

- **French National Stranding Network** for data collection
- **Observatoire Pelagis** (UAR 3462 CNRS-La Rochelle Université) for data curation
- **CNRS SEE-Life program** for long-term support
- **French Ministry of Ecological Transition** for funding
- **Field teams and interns** for sample collection and processing

## ⚖️ License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🚨 Conservation Implications

This research provides critical evidence for:
- **Immediate action needed** in the Bay of Biscay
- **Early warning system** for cetacean populations
- **Policy-relevant indicators** for MSFD implementation
- **Proactive conservation** rather than reactive management

---

🐬 *"In the race between dolphin decline and conservation action, early warning systems like this can mean the difference between recovery and extinction."*

**Last updated:** January 2025  
**Repository maintained by:** Etienne Rouby
