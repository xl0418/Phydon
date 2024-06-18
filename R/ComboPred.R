#' @title ComboPred
#' @description This function estimates the maximum growth rate of a species making use of both gRodon and phylogenetic methods by regression.
#' @param input_df The dataframe that contains the gRodon prediction, phylogenetic distance, and phylogenetic prediction.
#' @param reg_model The regression model that is trained on the gRodon and phylogenetic predictions.
#'
#' @noRd



combopred <- function(input_df, reg_model) {
  test_data <- data.frame(gRodon = input_df$gRodonpred,
                          phy_distance = input_df$phy_distance)
  test_data$pred_fitted <- stats::predict(reg_model, test_data, type = "response")

  input_df$combopred <- input_df$gRodonpred * (test_data$pred_fitted > 0.5) + input_df$phylopred * (test_data$pred_fitted <= 0.5)
  return(input_df)

}
