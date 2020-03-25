
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

#    data("isotopes",
#         package = "enviPat",
#         envir = environment())

    elements.l = check_chemform(isotopes = isotopes,
                                chemforms = chem.forms,
                                get_list = TRUE)
    chem.forms.sorted = sapply(elements.l, function(x) {
        if (length(x) == 0) return(NA)
        x = x[order(names(x))]
        paste(sapply(1:length(x), function(y) {
            paste0(names(x)[y], x[y])}), collapse = "")
    })
    return(chem.forms.sorted)
}
