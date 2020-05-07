

#' Correct chemical formulas
#'
#' Orders elements in chemical formulas alphabetically and adds 1 for elements
#' that occur only once.
#'
#' @param chem.forms Vector of character strings with chemical formulas.
#' @param isotopes dataframe with stable isotopes (see
#' \code{\link[enviPat]{isotopes}})
#'
#' @return Vector of character strings with corrected chemical formulas
#'
#' @import enviPat
#' @importFrom utils data
#' @export
#'
#' @examples
#' library(enviPat)
#' data("isotopes")
#'
#' correct_chem_formula(chem.forms = c("C44H92NO6P", "C6H12O6"),
#'                      isotopes = isotopes)

correct_chem_formula <- function(chem.forms,
                                 isotopes) {

    elements.l = check_chemform(isotopes = isotopes,
                                chemforms = chem.forms,
                                get_list = TRUE)
    chem.forms.sorted = vapply(
        X = elements.l,
        FUN = function(x) {
            if (length(x) == 0)
                return(NA)
            x = x[order(names(x))]
            paste(vapply(
                X = seq_len(length(x)),
                FUN = function(y) {
                    paste0(names(x)[y], x[y])
                },
                FUN.VALUE = character(1)
            ),
            collapse = "")
        },
        FUN.VALUE = character(1)
    )
    return(chem.forms.sorted)
}
