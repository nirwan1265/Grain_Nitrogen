library(vroom)
library(dplyr)

### Read the table:
table <- vroom("results/tables/Annotation_table/Final_Summary_GWAS_Results_MLM_MLMM_0.05_BLINK_FarmCPU_7_GLM_25.csv")
annotation <- vroom("results/tables/Annotation_table/annotation_Final_summary_GWAS_Results_MLM_MLMM_0.05_BLINK_FarmCPU_7_GLM_25.txt")

# Join the annotation to the table:
table <- table %>%
  left_join(annotation)

# Website
table$website <- paste0("http://maizegdb.org/gene_center/gene/",table$GeneID)

write.csv(table,"full_annotation_MLM_MLMM_BLINK_FarmCPU_GLM.csv", row.names = FALSE, quote = TRUE)
