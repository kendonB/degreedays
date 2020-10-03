# JT REVIEW: add author, description/purpose of code, creation date and last modified date

plot_spl1 <- function(reg_lm, y_lab = "Test",
                      title_main = NULL,
                      var_select = NULL,
                      text_scale_fact = 1,
                      y_lim = NULL,
                      x_lim = c(-5, 35),
                      x_trans = 0,
                      from_tp = 0,
                      interaction_str = NULL,
                      interaction_evaluation_value = NULL,
                      ...){
  x1 <- seq(x_lim[1], x_lim[2], length.out = 1000)
  
  varnames <- rownames(reg_lm$coefficients) %>%
    grep(pattern = "(Intercept)", x = ., invert = TRUE, value = TRUE)
  # Filter to get the variables we want.
  if(!is.null(var_select)){
    for(elem in var_select){
      varnames <- grep(elem, varnames, value = TRUE)
    }
  }
  
  # Detect interaction terms
  if(any(str_detect(varnames, "_x_"))){
    interaction_vars <- varnames %>% 
      stringr::str_subset("_x_") %>% 
      stringr::str_split("_x_") %>%
      unlist() %>% 
      unique() %>%   
      .[!str_detect(., var_select)]
    if(is.null(interaction_evaluation_value)){
      # Add to the coefficients using the interaction term's mean value
      interaction_evaluation_value <- tibble::as_tibble(reg_lm$cX) %>% 
        dplyr::select(period_length, interaction_vars) %>% 
        summarize_at(.vars = vars(-period_length),
                     list(mean = ~sum(.)/sum(period_length))) %>% 
        as.matrix() %>% 
        .[1, ] %>% 
        as.double()
    } else if (length(interaction_evaluation_value) == 1){
      interaction_evaluation_value <- rep(interaction_evaluation_value, 
                                          length(interaction_vars))
    }
    varnames <- varnames %>% 
      .[!str_detect(., "_x_")]
  } else {
    interaction_vars <- character(0)
  }
  
  # Create the newdata variables.
  newdata <- data.frame(matrix(ncol = length(reg_lm$coefficients),
                               nrow = length(x1), 0))
  names(newdata) <- rownames(reg_lm$coefficients)
  
  dd_knots <- get_temps_bottom(temp_names = varnames)
  dd_knots <- dd_knots[dd_knots != Inf]
  if(length(dd_knots) != length(varnames)){
    stop("length(dd_knots) != length(varnames). Check plot_dd.")
  }
  for(i in 1:length(varnames)){
    # Assumes they're in the correct order.
    newdata[,varnames[i]] <- x1
    if(dd_knots[i] > -Inf){
      newdata[,varnames[i]][x1 < dd_knots[i]] <- 0
      newdata[,varnames[i]][x1 >= dd_knots[i]] <-
        newdata[,varnames[i]][x1 >= dd_knots[i]] - dd_knots[i]
    }
  }
  
  if(from_tp > 0){
    x_trans <- dd_knots[min(from_tp + 1, length(dd_knots))]
  }
  
  for(i in 1:length(varnames)){
    # Assumes they're in the correct order.
    if(dd_knots[i] < x_trans & dd_knots[i] > -Inf){
      newdata[,varnames[i]] <- newdata[,varnames[i]] - x_trans + dd_knots[i]
    } else if (dd_knots[i] < x_trans) {
      newdata[,varnames[i]] <- newdata[,varnames[i]] - x_trans
    }
  }
  tmp <- list(...)
  if(length(tmp) > 0){
    for(i in 1:length(tmp)){
      newdata[[names(tmp)[i]]] <- tmp[[i]]
    }
  }
  if(length(interaction_vars) > 0){
    var_select <- paste0(var_select, paste0("|_x_", gsub("\\^", "", var_select)))
    for(j in interaction_vars %>% seq_along()){
      var <- interaction_vars[j]
      var <- paste0("^", var, paste0("_x_|_x_", gsub("\\^", "", var)))
      
      main_variables_interacted <- names(newdata) %>% 
        str_subset(var_select) %>% 
        str_subset(var) %>% 
        str_subset("_x_") %>% 
        str_replace(var, "") %>% 
        str_replace("_x_", "")
      
      interaction_columns_logical <- stringr::str_detect(names(newdata), var) &
        stringr::str_detect(names(newdata), var_select) &
        stringr::str_detect(names(newdata), "_x_")
      
      bare_columns_logical <- stringr::str_detect(names(newdata), var_select) & 
        !stringr::str_detect(names(newdata), "_x_") & 
        stringr::str_detect(names(newdata), paste0(main_variables_interacted, collapse = "|"))
      names(newdata)[bare_columns_logical]
      #newdata[,var] <- interaction_value[i]
      newdata[, interaction_columns_logical] <- newdata[, bare_columns_logical]*interaction_evaluation_value[j]
    }
  }
  pred_out <- predict(reg_lm, newdata,
                      interval = "confidence")
  pred_out$lwr <- Re(pred_out$lwr)
  pred_out$upr <- Re(pred_out$upr)
  pred_out$x <- x1
  pred_out$lwr[is.nan(pred_out$lwr)] <- pred_out$fit[is.nan(pred_out$lwr)]
  pred_out$upr[is.nan(pred_out$upr)] <- pred_out$fit[is.nan(pred_out$upr)]
  if(length(interaction_vars) > 0){
    pred_out <- rbind(pred_out, 0)
    pred_out$x[NROW(pred_out)] <- x_trans
  }
  env <- new.env(parent = globalenv())
  env$pred_out <- pred_out
  env$text_scale_fact <- text_scale_fact
  env$y_lab <- y_lab
  env$title_main <- title_main
  env$y_lim <- y_lim
  env$x_lim <- x_lim
  with(env, {
    graph_out <- pred_out %>% ggplot2::ggplot(ggplot2::aes(x = x)) +
      ggplot2::geom_line(ggplot2::aes(y = fit), size = 1, colour = "black") +
      ggplot2::geom_ribbon(ggplot2::aes(ymax = upr, ymin = lwr), alpha = 0.3, fill = "#E6E6E6") +
      ggplot2::geom_hline(ggplot2::aes(yintercept=0), linetype = "dashed", size = .5, color = "grey50") +
      ggplot2::theme(panel.border = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(colour = "black", size = .5),
                     axis.text.y=ggplot2::element_text(size = 12*text_scale_fact),
                     axis.title.y=ggplot2::element_text(size = 12*text_scale_fact,
                                                        vjust = 0.2),
                     axis.title.x=ggplot2::element_blank()
      ) +
      ggplot2::ylab(y_lab) +
      ggplot2::ggtitle(title_main) +
      ggplot2::coord_cartesian(xlim = x_lim) +
      ggplot2::scale_x_continuous(expand = c(0, 0))
    if(!is.null(y_lim)){
      graph_out <- graph_out + ggplot2::coord_cartesian(xlim = x_lim,
                                                        ylim = y_lim)
    }
    graph_out$df <- pred_out
    graph_out
  })
  env$graph_out
}

