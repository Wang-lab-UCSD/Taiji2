# make the first letter uppercase and other letters as lower case
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(tolower(x), 1, 1))
  x
}

