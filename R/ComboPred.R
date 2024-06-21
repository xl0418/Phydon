#' @title ComboPred
#' @description This function estimates the maximum growth rate of a species making use of both gRodon and phylogenetic methods by regression.
#' @param input_df The dataframe that contains the gRodon prediction, phylogenetic distance, and phylogenetic prediction.
#'
#' @noRd



combopred <- function(input_df) {
  test_data <- data.frame(log_gRodon = log(input_df$gRodonpred),
                          phy_distance = input_df$phy_distance)
  # if input_df has tmp column
  if("temp" %in% colnames(input_df)) {
    test_data$pred_fitted <- stats::predict(reg_model_tmp, test_data, type = "response")
  } else {
    test_data$pred_fitted <- stats::predict(reg_model, test_data, type = "response")
  }

  input_df$combopred <- input_df$gRodonpred * (test_data$pred_fitted > 0.5) + input_df$phylopred * (test_data$pred_fitted <= 0.5)
  return(input_df)

}
