has_attribute <- function(key) {
  function(x)
    expectation(key %in% names(attributes(x)),
                paste("does not have attribute", key))
}

is_null <- function() {
  function(x)
    expectation(is.null(x), "is not NULL")
}

