### Superseded by iRfcb::ifcb_get_dashboard_metadata()


#' Download IFCB metadata from the dashboard API
#'
#' @param base_url Character. Base URL to the IFCB Dashboard (e.g. "https://ifcb-data.whoi.edu/").
#' @param dataset_name Optional character. Dataset slug (e.g. "RV_Svea") to retrieve metadata for a specific dataset.
#' If NULL, all available metadata are downloaded.
#' @param quiet Logical. If TRUE, suppresses progress messages. Default is FALSE.
#'
#' @return A data frame containing the exported metadata.
#' @examples
#' \dontrun{
#'   # Download all metadata
#'   metadata_all <- ifcb_get_dashboard_metadata("https://ifcb-data.whoi.edu/")
#'
#'   # Download metadata for a specific dataset
#'   metadata_svea <- ifcb_get_dashboard_metadata("https://ifcb-data.whoi.edu/", "mvco")
#' }
#' @export
ifcb_get_dashboard_metadata <- function(base_url, dataset_name = NULL, quiet = FALSE) {
  # Ensure base_url has no trailing slash
  base_url <- sub("/+$", "", base_url)
  
  # Build API URL
  api_url <- paste0(base_url, "/api/export_metadata/")
  if (!is.null(dataset_name) && nzchar(dataset_name)) {
    # Safely encode dataset name (handles special chars)
    dataset_name <- utils::URLencode(dataset_name, reserved = TRUE)
    api_url <- paste0(api_url, dataset_name)
  }
  
  if (!quiet) message("Fetching metadata from: ", api_url)
  
  # Perform GET request
  response <- tryCatch(
    httr::GET(api_url, httr::accept("text/csv")),
    error = function(e) stop("Failed to connect to IFCB Dashboard API: ", e$message)
  )
  
  # Check response status
  if (httr::status_code(response) != 200) {
    stop("API request failed [", httr::status_code(response), "]: ", api_url)
  }
  
  # Extract and parse CSV
  csv_content <- httr::content(response, as = "text", encoding = "UTF-8")
  df <- tryCatch(
    read.csv(text = csv_content, stringsAsFactors = FALSE),
    error = function(e) stop("Failed to parse CSV content: ", e$message)
  )
  
  if (!quiet) {
    message("Successfully retrieved ", nrow(df), " records",
            if (!is.null(dataset_name)) paste0(" for dataset '", dataset_name, "'"), ".")
  }
  
  return(df)
}

#' Download IFCB bin list from the dashboard API
#'
#' @param base_url Character. Base URL to the IFCB Dashboard
#'   (e.g. "https://ifcb-data.whoi.edu/").
#' @param quiet Logical. If TRUE, suppresses progress messages. Default is FALSE.
#'
#' @return A data frame containing the bin list returned by the API.
#' @examples
#' \dontrun{
#'   bins <- ifcb_list_dashboard_bins("https://ifcb-data.whoi.edu/")
#'   head(bins)
#' }
#' @export
ifcb_list_dashboard_bins <- function(base_url, quiet = FALSE) {
  # Ensure base_url has no trailing slash
  base_url <- sub("/+$", "", base_url)
  
  # Construct full API URL
  api_url <- paste0(base_url, "/api/list_bins")
  
  if (!quiet) message("Fetching bin list from: ", api_url)
  
  # Perform GET request
  response <- tryCatch(
    httr::GET(api_url, httr::accept("application/json")),
    error = function(e) stop("Failed to connect to IFCB Dashboard API: ", e$message)
  )
  
  # Check status code
  if (httr::status_code(response) != 200) {
    stop("API request failed [", httr::status_code(response), "]: ", api_url)
  }
  
  # Parse JSON content
  json_content <- httr::content(response, as = "text", encoding = "UTF-8")
  parsed_data <- tryCatch(
    jsonlite::fromJSON(json_content, flatten = TRUE),
    error = function(e) stop("Failed to parse JSON content: ", e$message)
  )
  
  # Convert to data frame if necessary
  df <- if (is.data.frame(parsed_data)) parsed_data else as.data.frame(parsed_data)
  
  if (!quiet) message("Successfully retrieved ", nrow(df), " bins.")
  
  return(df)
}
