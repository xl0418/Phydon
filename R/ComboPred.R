#' @title ComboPred
#' @description This function estimates the maximum growth rate of a species making use of both gRodon and phylogenetic methods by regression.
#' @param input_df The dataframe that contains the gRodon prediction, phylogenetic distance, and phylogenetic prediction.
#' @noRd



combopred <- function(input_df, regression_mode="arithmetic_mean") {
  # find the rows with any of gRodonpred, phy_distance, or phylopred is missing
  missing_rows <- which(is.na(input_df$gRodonpred) | is.na(input_df$phy_distance) | is.na(input_df$phylopred))
  if (length(missing_rows) > 0) {
    input_df_missing <- input_df[missing_rows,]
    input_df_missing$combopred <- NA
    # find the rows with all of gRodonpred, phy_distance, and phylopred are not missing
    input_df <- input_df[-missing_rows,]
  }



  if(nrow(input_df) > 0) {
    test_data <- data.frame(log_gRodon = log(input_df$gRodonpred),
                            phy_distance = input_df$phy_distance)
    # if input_df has tmp column
    if("temp" %in% colnames(input_df)) {
      test_data$pred_fitted <- stats::predict(reg_model_tmp, test_data, type = "response")
    } else {
      test_data$pred_fitted <- stats::predict(reg_model, test_data, type = "response")
    }
    if(regression_mode == "arithmetic_mean") {
      input_df$combopred <- input_df$gRodonpred * (test_data$pred_fitted) + input_df$phylopred * (1 - test_data$pred_fitted)
    } else if(regression_mode == "geometric_mean") {

      input_df$combopred <- input_df$gRodonpred ^ (test_data$pred_fitted) * input_df$phylopred ^ (1 - test_data$pred_fitted)
    }
  }

  if (length(missing_rows) > 0) {
    input_df <- rbind(input_df, input_df_missing)
  }


  return(input_df)

}
