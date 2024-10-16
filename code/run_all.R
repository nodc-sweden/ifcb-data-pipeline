# Load required package
library(rmarkdown)

# Define parameters for each classifier and year combination
params_list <- list(
  list(classifier = "Baltic", 
       years = c(2022, 2023, 2024), 
       remove_flagged_data = c("bubbles", "incomplete", "near land"), 
       multiyear_delivery = TRUE),
  list(classifier = "Skagerrak-Kattegat", 
       years = c(2022, 2023, 2024), 
       remove_flagged_data = c("bubbles", "incomplete", "near land"), 
       multiyear_delivery = TRUE),
  list(classifier = "Tångesund", 
       years = c(2016), 
       remove_flagged_data = c("bubbles", "incomplete"), 
       multiyear_delivery = FALSE)
)

# Loop through each classifier and render the RMarkdown file with the list of years and other parameters
for (params in params_list) {
  
  # Render the RMarkdown file with different parameters
  render(
    input = "run_all.Rmd",
    output_file = paste0("output/html_reports/ifcb_data_export_report_",
                         params$classifier, "_", paste(params$years, collapse = "_"), "_", format(Sys.time(), "D%Y%m%dT%H%M"), ".html"),
    params = list(
      classifier = params$classifier,
      year = params$years,  # Pass the list of years directly
      threshold = "opt",  # Assuming the threshold is "opt" for all
      remove_flagged_data = params$remove_flagged_data,  # Custom remove_flagged_data
      regional = TRUE,
      multiyear_delivery = params$multiyear_delivery,  # Custom multiyear_delivery
      f1_threshold = 0.9  # Assuming this remains constant
    ),
    encoding = "UTF-8"
  )
}
