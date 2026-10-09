devtools::load_all()

devtools::document()

library(protqc)

## example data
example_prot_data_path <- system.file("extdata", "proteomics_pipeline_data_template.csv", package = "protqc")
example_prot_metadata_path <- system.file("extdata", "proteomics_pipeline_meta_template.csv", package = "protqc")

## Calculate QC metrics
prot_result <- protqc::qc_conclusion(example_prot_data_path, example_prot_metadata_path)

prot_result

## Generate report
report_template <- system.file("extdata", "quartet_template.docx", package = "protqc")
report_template

prot_result$conclusion

generate_protein_report(qc_result=prot_result, report_template = report_template,report_dir = './test/')

data("reference_dataset", package = "protqc")
