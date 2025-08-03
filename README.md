## Differential gene expression of CCN2 Expression in Ovarian Lesions

### Overview
This R script performs a differential gene expression analysis comparing CCN2 expression between normal (Incidental Fimbriae) and transformed (Incidental STIC) ovarian lesions in epithelial cells. The workflow includes data loading, preprocessing, filtering, statistical testing, and visualization.

---

### Requirements

**R Packages**
- `tidyverse` (including `dplyr`, `ggplot2`, `tidyr`, `readr`)
- `readxl` (for Excel file handling)
- `stringr` (for string manipulation)
- `broom` (for tidy statistical output)

Install missing packages with:
```r
install.packages(c("tidyverse", "readxl", "stringr", "broom"))
```

---

### Input Files

- `Kader et al_SuppTable S8_Q3 norm_countmatrix.xlsx`  
  Sheet: `TargetCountMatrix` (contains gene expression data)
- `data_clinical_sample.xlsx`  
  Contains clinical metadata (e.g., lesion types, ROI types)

---

### Workflow Steps

#### 1. Data Loading & Preprocessing

**Target Count Matrix:**  
Loads and filters columns to retain only "Full ROI" samples.
```r
FullROI_data <- target_counts %>% select(TargetName, contains("| Full ROI"))
```

**Clinical Data:**  
- Standardizes sample IDs (replaces `-` with `|`, cleans "Full ROI" labels)
- Filters samples for:
  - Lesion types: Incidental Fimbriae (normal) and Incidental STIC (transformed)
  - ROI type: Epithelial

#### 2. Sample Matching & Filtering

- Matches clinical samples with expression data using `intersect()`
- **Safety check:** Stops if no matches are found

```r
matched_columns <- intersect(filtered_data$`#Sample Identifier`, colnames(FullROI_data))
```

#### 3. CCN2-Specific Analysis

- Extracts CCN2 expression values and merges with clinical status
- Output: A tidy dataframe (`ccn2_data`) with columns:  
  **SampleID**, **Expression**, **status** (normal/transformed)

#### 4. Statistical Testing

- Performs a t-test to compare CCN2 expression between groups
- Calculates p-value and formats it for plotting

#### 5. Visualization

- Generates a publication-ready boxplot:
  - **Boxes:** Median/IQR expression by status
  - **Points:** Individual sample expression (jittered)
  - **Annotation:** p-value from t-test
- **Customization:**
  - Color scheme: Green (normal) vs. orange (transformed)
  - Minimalist theme with bold title and visible axis lines

```r
ggsave("CCN2_expression_boxplot (STIC).png", width = 8, height = 6, dpi = 300)
```

---

### Key Outputs

- **Data Files** (commented out by default):
  - `FullROI_only.csv`: Filtered expression matrix
  - `ccn2_filtered.csv`: CCN2 expression with clinical status
- **Plot:**
  - `CCN2_expression_boxplot (STIC).png`: Boxplot comparing expression

---

### Interpretation

- **Statistical Result:**  
  High p-value (e.g., p = 0.85) suggests no significant difference in CCN2 expression between normal and transformed lesions.
- **Biological Implication:**  
  CCN2 may not be a key driver of early transformation in STIC lesions. Further validation (e.g., protein-level assays) is recommended.

---

### Customization

- **Lesion Types:** Modify the `filter()` step to include/exclude groups (e.g., Incidental p53sig)
- **Genes:** Replace `CCN2` with other targets in `filter(TargetName == ...)`
- **Plot:** Adjust colors, labels, or themes in the `ggplot2` code

---

### Troubleshooting

- **Error: "No matching columns":** Verify sample ID formats match between clinical and expression data
- **Missing Packages:** Install required packages (see Requirements)
- **File Paths:** Update `setwd()` to your project directory

---

### References

Kader et al **Multimodal spatial profiling reveals immune suppression and microenvironment remodeling in Fallopian tube precursors to high-grade serous ovarian carcinoma.**
**PMCID: PMC12130810, PMID: 39704522** 

---

### Contact

For questions, contact: **Aishwarya Mahale** at **aishamahale30@gmail.com**  
