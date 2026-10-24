# Data Structure

Each directory contains data used to run DM[H]N for a different study:

- **`ZIM`**: Zimbabwe cohort
  - Burkhart Colorado, A. S., Lazzaro, A., Neff, C. P., et al. (2024). Differential effects of antiretroviral treatment on immunity and gut microbiome composition in people living with HIV in rural versus urban Zimbabwe. *Microbiome, 12*, 18. https://doi.org/10.1186/s40168-023-01718-4
  - Burkhart Colorado, A. S., Nusbacher, N. M., O'Connor, J., Marden, T., Higgins, J., Neff, C. P., Fiorillo, S., Campbell, T. B., Borok, M., Boyd, K., Sterrett, J., Palmer, B. E., & Lozupone, C. (2026). The impact of western versus agrarian diet consumption on gut microbiome composition and immune dysfunction in people living with HIV in rural and urban Zimbabwe. *Microbiome, 14*(1), 179. https://doi.org/10.1186/s40168-026-02410-z

- **`DM`**: Diet Modification cohort
  - O'Connor, J. B., Fouquier, J., Neff, C. P., Sterrett, J. D., Marden, T., Fiorillo, S., Siebert, J. C., Schneider, J., Nusbacher, N., Noe, A. T., Fennimore, B., Higgins, J., Campbell, T. B., Palmer, B. E., & Lozupone, C. (2025). Agrarian diet improves metabolic health in HIV-positive men with *Prevotella*-rich microbiomes: Results from a randomized trial. *mSystems, 10*(12), e0118525. https://doi.org/10.1128/msystems.01185-25

- **`IBD`**: Inflammatory Bowel Disease cohort
  - Lloyd-Price, J., Arze, C., Ananthakrishnan, A. N., et al. (2019). Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases. *Nature, 569*, 655–662. https://doi.org/10.1038/s41586-019-1237-9

Each cohort directory contains a `participants` file with information about individual participants. Because DM[H]N operates on a per-sample basis, participant information is associated with the corresponding samples. Each cohort directory may also contain additional files specific to that study.
