# Contributing to rnaseqDegy

Thank you for your interest in contributing to **rnaseqDegy**!  
We welcome contributions from researchers, students, developers, and anyone interested in improving bulk RNA‑seq workflows.

This document describes the contribution process, coding standards, and expectations for all contributors.

---

## 🧭 How to Contribute

### 1. **Fork the Repository**
Navigate to:  
**https://github.com/ebareke/rnaseqDegy**  
Click **Fork** to create your own copy under your GitHub account.

### 2. **Clone Your Fork**
```bash
git clone https://github.com/<your-username>/rnaseqDegy.git
cd rnaseqDegy
```

### 3. **Create a New Feature Branch**
```bash
git checkout -b feature/my-new-feature
```
Use meaningful branch names (e.g., `bugfix/pca-labels`, `feature/add-kegg-support`).

### 4. **Install Development Environment**
You'll need:
- R ≥ 4.2
- RStudio (recommended)
- **devtools**, **roxygen2**, **BiocManager**

Install dependencies:
```r
install.packages(c("devtools", "roxygen2"))
BiocManager::install(c("DESeq2","edgeR","limma","clusterProfiler","enrichplot"))
```

### 5. **Build & Load Package Locally**
```r
devtools::load_all()
```
This loads the package without installing it.

### 6. **Run Tests** (if tests exist)
```r
devtools::test()
```

### 7. **Document Your Changes**
Whenever you update functions:
```r
devtools::document()
```
This will update **NAMESPACE** and `.Rd` files.

### 8. **Commit Your Changes**
```bash
git add .
git commit -m "Added new PCA auto-scaling feature"
```
Write **clear and concise commits**.

### 9. **Push and Open a Pull Request**
```bash
git push origin feature/my-new-feature
```
Then go to GitHub and open a **Pull Request** (PR).

Provide:
- A clear description of the change
- Motivation (why this matters)
- Any visual examples (plots, error messages, etc.)
- Impact on downstream functionality

---

## 🧪 Coding Guidelines

### ✔ R Style
- Follow **tidyverse** style conventions
- Use meaningful variable names (`vst_mat`, `norm_counts`, etc.)
- Avoid deeply nested logic where possible
- Add inline comments for clarity

### ✔ Roxygen Documentation
Every exported function **must** have:
- Title
- Description
- Arguments
- Value
- Examples where appropriate

### ✔ Avoid Hard-Coded Paths
All paths must be function arguments or constructed from `output_dir`.

### ✔ Reproducibility
Ensure:
- Deterministic output where possible
- Documented random seeds
- Clear logging (`message()` calls)

---

## 🧩 Bug Reports & Issues
Please use the GitHub Issue Tracker:  
👉 https://github.com/ebareke/rnaseqDegy/issues

Include:
- Operating system
- R version
- Package version
- Error logs
- Minimal reproducible example if possible

---

## 🤝 Code of Conduct
All contributors must follow the community guidelines described in `CODE_OF_CONDUCT.md`.

We strive to provide a welcoming and inclusive research environment.

---

## 📝 Licensing
By contributing, you agree that your contributions will be licensed under the **MIT License**, the same license governing rnaseqDegy.

---

## 🙏 Acknowledgements
Thank you for contributing to **rnaseqDegy** and supporting open, reproducible, high‑quality RNA‑seq analysis.

If you have suggestions or want to propose a feature, do not hesitate to open an issue!
