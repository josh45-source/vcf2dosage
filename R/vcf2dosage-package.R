#' @keywords internal
"_PACKAGE"

#' @importFrom rlang .data %||% abort inform warn
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr bind_cols bind_rows mutate select filter summarise arrange
#'   n distinct pull across left_join rename group_by ungroup
#' @importFrom purrr map map_chr map_dbl map_int map_lgl compact
#' @importFrom glue glue
#' @importFrom cli cli_h1 cli_h2 cli_h3 cli_ul cli_alert_success
#'   cli_alert_info cli_alert_warning cli_alert_danger cli_abort
#'   cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom data.table fread
#' @importFrom stats complete.cases setNames
#' @importFrom utils head tail read.table write.csv write.table
NULL
