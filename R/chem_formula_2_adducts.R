

#' Get adducts and isotopes for chemical formulas
#'
#' Calculates chemical formulas, mz values and theoretical abundances of
#' isotope peaks of adducts of chemical formulas. Based on functions in R
#' package enviPat.
#'
#' @param chem.forms Vector of character strings with chemical formulas.
#' @param isotopes dataframe with stable isotopes (see
#' \code{\link[enviPat]{isotopes}})
#' @param adducts dataframe with common adducts (see
#' \code{\link[enviPat]{adducts}}).
#' @param adduct.names Character vector with adducts to be calculated
#' (see \code{\link[enviPat]{adducts}} for common adducts). If set to NULL,
#' ion.mode needs to be specified.
#' @param ion.mode Character of ionisation mode ("positive", "negative"). If
#' positive adduct.names is set to c("M+H", "M+Na", "M+K"), if negative
#' adduct.names is set to c("M-H", "M+Cl").
#' @param rel_to Integer (0, 1, 2, 3 or 4). Probability definition
#' (see \code{\link[enviPat]{isopattern}}).
#' @param threshold Numeric (1-100). Probability below which isotope peaks can
#' be omitted (see \code{\link[enviPat]{isopattern}}).
#' @param verbose Logical. Print out additional information.
#'
#' @return data frame with the following information about each feature:
#' \itemize{
#'   \item chemical.form.adduct: chemical formula of adduct
#'   \item chemical.form: original chemical formula
#'   \item adduct.name: adduct used
#'   \item charge: charge of adduct
#'   \item chemical.form.isotope: chemical formula of isotope of adduct
#'   \item m.z: theoretical mz value of isotope
#'   \item abundance: theoretical abundance of isotope peak. Most abundant
#'   isotope has value of 100
#'   \item id: unique identifier (<chemical.form>.<chemical.form.isotope>)
#' }
#'
#' @import enviPat
#' @importFrom stats na.omit
#' @export
#'
#' @examples
#' library(enviPat)
#' data("isotopes")
#' data("adducts")
#'
#' chem_formula_2_adducts(chem.forms = "C6H12O6",
#'                        adduct.names = c("M+H", "M+Na", "M+2H"),
#'                        adducts = adducts,
#'                        isotopes = isotopes,
#'                        threshold = 5)

chem_formula_2_adducts <- function(chem.forms,
                                   isotopes,
                                   adducts,
                                   adduct.names = NULL,
                                   ion.mode = "positive",
                                   rel_to = 0,
                                   threshold = 20,
                                   verbose = TRUE) {
    if (is.null(adduct.names)) {
        if (is.null(ion.mode)) {
            stop("no adduct names given!")
        } else {
            if (ion.mode == "positive") {
                adduct.names = c("M+H", "M+Na", "M+K")
            } else if (ion.mode == "negative") {
                adduct.names = c("M-H", "M+Cl")
            } else {
                stop("ion.mode needs to be 'positive' or 'negative'!")
            }
        }
    }
    if (verbose) {
        print("using the following adducts:")
        print(paste(adduct.names, collapse = ", "))
    }

    ## restrict to unique chemical formula
    chem.forms = unique(chem.forms)

    ## remove NAs
    chem.forms = na.omit(chem.forms)

    ## check chemical formula
    res.check = check_chemform(isotopes, chem.forms)

    ## remove formula without mass
    ind = which(res.check$monoisotopic_mass != -9999)
    chem.forms = unique(res.check[ind, "new_formula"])
    if (verbose)
        print(paste(length(chem.forms), "unique chemical formula"))

    info.patterns = NULL
    for (a in adduct.names) {
        if (verbose)
            print(paste("get adducts for", a))
        info.patterns = rbind(
            info.patterns,
            get_info_adduct(
                chem.forms = chem.forms,
                adduct.name = a,
                adducts = adducts,
                isotopes = isotopes,
                rel_to = rel_to,
                threshold = threshold
            )
        )
    }

    return(info.patterns)
}


#' Adduct calculations
#'
#' Based on functions in R package enviPat.
#'
#' @param chem.forms Vector of character strings with chemical formulas.
#' @param adduct.name Character. Name of adduct to be calculated (see
#' \code{\link[enviPat]{adducts}} for possible adduct names).
#' @param adducts dataframe with common adducts (see
#' \code{\link[enviPat]{adducts}}).
#' @param isotopes dataframe with stable isotopes (see
#' \code{\link[enviPat]{isotopes}}).
#' @param threshold Numeric (1-100). Probability below which isotope peaks can
#' be omitted (see \code{\link[enviPat]{isopattern}}).
#' @param rel_to Integer (0, 1, 2, 3 or 4). Probability definition
#' (see \code{\link[enviPat]{isopattern}}).
#'
#' @return data frame with the following information about each feature:
#' \itemize{
#'   \item chemical.form.adduct: chemical formula of adduct
#'   \item chemical.form: original chemical formula
#'   \item adduct.name: adduct used
#'   \item charge: charge of adduct
#'   \item chemical.form.isotope: chemical formula of isotope of adduct
#'   \item m.z: theoretical mz value of isotope
#'   \item abundance: theoretical abundance of isotope peak. Most abundant
#'   isotope has value of 100
#'   \item id: unique identifier (<chemical.form>.<chemical.form.isotope>)
#' }
#'
#' @import enviPat
#' @keywords internal

get_info_adduct <- function(chem.forms,
                            adduct.name = "M+H",
                            adducts,
                            isotopes,
                            threshold = 20,
                            rel_to = 0) {
    adducts$Formula_add[which(adducts$Formula_add == "FALSE")] = NA
    adducts$Formula_ded[which(adducts$Formula_ded == "FALSE")] = NA

    ## extract info for given adduct
    info.adduct = adducts[which(adducts$Name == adduct.name),]
    if (nrow(info.adduct) == 0) {
        stop(paste("adduct.name", adduct.name, "not found!"))
    }
    if (nrow(info.adduct) > 1) {
        stop(paste(adduct.name, "not unique!"))
    }

    ## multiplication and addition of adduct
    chem.forms.adducts = multiform(formula1 = chem.forms,
                                   fact = info.adduct$Mult)
    if (!is.na(info.adduct$Formula_add)) {
        chem.forms.adducts = mergeform(formula1 = chem.forms.adducts,
                                       formula2 = info.adduct$Formula_add)
    }
    if (!is.na(info.adduct$Formula_ded)) {
        chem.forms.adducts = subform(formula1 = chem.forms.adducts,
                                     formula2 = info.adduct$Formula_ded)
    }

    ## sort by element
    chem.forms.adducts.sorted = correct_chem_formula(chem.forms =
                                                         chem.forms.adducts,
                                                     isotopes = isotopes)

    info = data.frame(
        chemical.form = chem.forms,
        chemical.form.adduct = chem.forms.adducts.sorted,
        adduct.name = rep(adduct.name, length(chem.forms)),
        charge = info.adduct$Charge,
        stringsAsFactors = FALSE
    )

    ## isotopes
    info.isotopes = get_info_isotopes(
        chem.forms = chem.forms.adducts.sorted,
        isotopes = isotopes,
        rel_to = rel_to,
        threshold = threshold
    )
    info.final = merge(info, info.isotopes,
                       all.y = TRUE)

    ## modify m.z depending on charge
    info.final$m.z = info.final$m.z / abs(info.final$charge)

    ## define ids
    info.final$id = paste(info.final$chemical.form,
                          info.final$chemical.form.isotope,
                          sep = ".")
    return(info.final)
}

# library(enviPat)
# data("isotopes")
# data("adducts")
#
# get_info_adduct(chem.forms = "C6H12O6",
#                 adduct.name = "M+H",
#                 adducts = adducts,
#                 isotopes = isotopes)



#' Isotope information
#'
#' Based on functions in R package enviPat.
#'
#' @param chem.forms Vector of character strings with chemical formulas.
#' @param isotopes dataframe with stable isotopes (see
#' \code{\link[enviPat]{isotopes}}).
#' @param threshold Numeric (1-100). Probability below which isotope peaks can
#' be omitted (see \code{\link[enviPat]{isopattern}}).
#' @param rel_to Integer (0, 1, 2, 3 or 4). Probability definition
#' (see \code{\link[enviPat]{isopattern}}).
#'
#' @return data frame with the following information about each feature:
#' \itemize{
#'   \item chemical.form.adduct: chemical formula of adduct
#'   \item chemical.form: original chemical formula
#'   \item adduct.name: adduct used
#'   \item charge: charge of adduct
#'   \item chemical.form.isotope: chemical formula of isotope of adduct
#'   \item m.z: theoretical mz value of isotope
#'   \item abundance: theoretical abundance of isotope peak. Most abundant
#'   isotope has value of 100
#'   \item id: unique identifier (<chemical.form>.<chemical.form.isotope>)
#' }
#'
#' @import enviPat
#' @keywords internal

get_info_isotopes <- function(chem.forms,
                              isotopes,
                              rel_to = 0,
                              threshold = 20) {
    patterns = isopattern(
        isotopes = isotopes,
        chemforms = chem.forms,
        rel_to = rel_to,
        threshold = threshold,
        verbose = FALSE
    )
    ind.error = grep("error", patterns)
    if (length(ind.error) > 0)
        patterns = patterns[-ind.error]

    if (length(patterns) == 0)
        return(NULL)

    ## change isotope names
    isotopes.names = isotopes[-grep("\\[|^D$", isotopes$element),]
    isotopes.names$element = vapply(
        X = isotopes.names$isotope,
        FUN = function(x) {
            pattern = "^[0-9]*"
            m = regexpr(pattern, x)
            pattern.found = regmatches(x, m)
            gsub(pattern.found, paste0("[", pattern.found, "]"), x)
        },
        FUN.VALUE = character(1)
    )
    isotopes.names$abundance = rep(1, nrow(isotopes.names))
    isotopes.names = unique(isotopes.names)
    rownames(isotopes.names) = isotopes.names$isotope

    info.patterns = NULL
    for (i in seq_len(length(patterns))) {
        temp = patterns[[i]]
        chem.forms.iso = NULL
        col.atoms = setdiff(colnames(temp), c("m/z", "abundance"))
        for (r in seq_len(nrow(temp))) {
            x = temp[r, col.atoms]
            ind = which(x > 0)
            cf = paste(vapply(
                X = ind,
                FUN = function(y) {
                    paste0(unique(isotopes.names[which(
                        isotopes.names$isotope == col.atoms[y]), "element"]),
                           x[y])
                },
                FUN.VALUE = character(1)
            ),
            collapse = "")
            chem.forms.iso[r] = correct_chem_formula(chem.forms = cf,
                                                     isotopes = isotopes.names)
        }

        info.iso = data.frame(
            chemical.form.adduct = rep(names(patterns)[i],
                                       nrow(temp)),
            chemical.form.isotope = chem.forms.iso,
            temp[, seq_len(2), drop = FALSE],
            stringsAsFactors = FALSE
        )
        info.patterns = rbind(info.patterns, info.iso)
    }

    return(info.patterns)

}

# library(enviPat)
# data("isotopes")
# data("adducts")
#
#get_info_isotopes(chem.forms = "C6H12O6",
#                  isotopes = isotopes,
#                  threshold = 1)
