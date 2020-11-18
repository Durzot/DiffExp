test_that("refine_beta_names works", {

  data <- data.frame(genotype=sample(c("I", "II", "III"), size=100, replace=T),
                     treatment=sample(c("A", "B", "C"), size=100, replace=T),
                     status=sample(c("S1", "S2", "S5"), size=100, replace=T),
                     stringsAsFactors=T)

  design <- ~ genotype + treatment
  beta_names <- refine_beta_names(design, data)
  expect_equal(beta_names, c("Intercept", "genotype_II_vs_I", "genotype_III_vs_I", "treatment_B_vs_A", 
                             "treatment_C_vs_A"))

  design <- ~ genotype +  genotype:treatment
  beta_names <- refine_beta_names(design, data)
  expect_equal(beta_names, c("Intercept","genotype_II_vs_I",
                             "genotype_III_vs_I","genotypeI.treatment_B_vs_A",
                             "genotype_II_vs_I.treatment_B_vs_A","genotype_III_vs_I.treatment_B_vs_A",
                             "genotypeI.treatment_C_vs_A","genotype_II_vs_I.treatment_C_vs_A",
                             "genotype_III_vs_I.treatment_C_vs_A"))


  design <- ~ genotype:treatment:status
  beta_names <- refine_beta_names(design, data)
  expect_equal(beta_names, c("Intercept", "genotypeI.treatmentA.statusS1",
                             "genotype_II_vs_I.treatmentA.statusS1",
                             "genotype_III_vs_I.treatmentA.statusS1",
                             "genotypeI.treatment_B_vs_A.statusS1",
                             "genotype_II_vs_I.treatment_B_vs_A.statusS1",
                             "genotype_III_vs_I.treatment_B_vs_A.statusS1",
                             "genotypeI.treatment_C_vs_A.statusS1",
                             "genotype_II_vs_I.treatment_C_vs_A.statusS1",
                             "genotype_III_vs_I.treatment_C_vs_A.statusS1",
                             "genotypeI.treatmentA.status_S2_vs_S1",
                             "genotype_II_vs_I.treatmentA.status_S2_vs_S1",
                             "genotype_III_vs_I.treatmentA.status_S2_vs_S1",
                             "genotypeI.treatment_B_vs_A.status_S2_vs_S1",
                             "genotype_II_vs_I.treatment_B_vs_A.status_S2_vs_S1",
                             "genotype_III_vs_I.treatment_B_vs_A.status_S2_vs_S1",
                             "genotypeI.treatment_C_vs_A.status_S2_vs_S1",
                             "genotype_II_vs_I.treatment_C_vs_A.status_S2_vs_S1",
                             "genotype_III_vs_I.treatment_C_vs_A.status_S2_vs_S1",
                             "genotypeI.treatmentA.status_S5_vs_S1",
                             "genotype_II_vs_I.treatmentA.status_S5_vs_S1",
                             "genotype_III_vs_I.treatmentA.status_S5_vs_S1",
                             "genotypeI.treatment_B_vs_A.status_S5_vs_S1",
                             "genotype_II_vs_I.treatment_B_vs_A.status_S5_vs_S1",
                             "genotype_III_vs_I.treatment_B_vs_A.status_S5_vs_S1",
                             "genotypeI.treatment_C_vs_A.status_S5_vs_S1",
                             "genotype_II_vs_I.treatment_C_vs_A.status_S5_vs_S1",
                             "genotype_III_vs_I.treatment_C_vs_A.status_S5_vs_S1"))
})
