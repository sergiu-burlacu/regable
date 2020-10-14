form_reg <- function(y, ..., add = NULL) {
  lhs <- enquo(y)
  rhs <- enquos(...)
  rhs <- paste(map(rhs, as_label), collapse = " + ")
  form <- paste(as_label(lhs), "~", rhs)
  if (!is.null(add)) {
    form <- paste(form, add)
  }
  return(as.formula(form))
}

base_reg <- function(data = data ,y, ..., add = NULL, fun = lm) {
  lhs <- enquo(y)
  rhs <- enquos(...)
  rhs <- paste(map(rhs, as_label), collapse = " + ")
  form <- paste(as_label(lhs), "~", rhs)
  if (!is.null(add)) {
    form <- paste(form, add)
  }
  fun(data = data, formula = as.formula(form))
}

base_reg2 <- function(data = data, y, ..., add = NULL, fun = lm) {
  lhs <- y
  rhs <- c(...) %>% paste(collapse = "+")
  form <- paste(lhs, "~", rhs)
  if (!is.null(add)) {
    form <- paste(form, add)
  }
  fun(data = data, formula = as.formula(form))
}
