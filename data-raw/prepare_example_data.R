## example data
##
## source:
## Metabolights MTBLS162
## three quality control samples
## reduced mz range (containing some of the standards)

library(preprocessHighResMS)
library(SummarizedExperiment)
library(enviPat)
data(isotopes)
data(adducts)

## features of standards
chem.forms = c("C40H81N2O6P1", "C44H92N1O6P1", "C45H94N1O6P1")
info.features = chem_formula_2_adducts(chem.forms = chem.forms,
                                       isotopes = isotopes,
                                       adducts = adducts,
                                       verbose = TRUE,
                                       adduct.names = c("M+H", "M+NH4"),
                                       rel_to = 0,
                                       threshold = 20)

usethis::use_data(info.features,
                  overwrite = TRUE)


## SE
data("info.features")
res.dir = tempdir()
mzml.files = dir(system.file("extdata",
                        package = "preprocessHighResMS"),
                 full.names = TRUE)

sapply(mzml.files[1:2],
       extract_feature_intensity,
       scanrange = c(1, 2),
       info.features = info.features,
       ppm = 20,
       res.dir = res.dir)

feature.files = dir(path = res.dir,
                    pattern = ".rds",
                    full.names = TRUE)
se.example = combine_feature_intensities(feature.files = feature.files,
                                         verbose = TRUE)

#peak = assays(se.example)$peak.detection

usethis::use_data(se.example,
                  overwrite = TRUE)

# #data("se.example")
# se.example = remove_features(se = se.example,
#                              assay = "mz",
#                              method = "identical.peaks")
#
#
# #data("se.example")
# dat = prepare_data_for_annotation(se = se.example)
#
#
#
# data("se.example")
# data("info.features")
#
# dat = prepare_data_for_annotation(se = se.example)
# hits.m = find_hits(info.features = info.features,
#                    dat = dat,
#                    ppm = 20)
#
# prior.prob = compute_prior_prob(hits.m = hits.m,
#                                 info.features = info.features,
#                                 dat = dat,
#                                 ppm = 20)
#
#
# data("info.features")
# add.m = generate_connectivity_matrix(info.features = info.features,
#                                      type = "adducts")
#
# iso.m = generate_connectivity_matrix(info.features = info.features,
#                                      type = "isotopes",
#                                      ratio = FALSE)
#
#
# data("se.example")
# data("info.features")
#
# dat = prepare_data_for_annotation(se = se.example)
# hits.m = find_hits(info.features = info.features,
#                    dat = dat,
#                    ppm = 20)
#
# prior.prob = compute_prior_prob(hits.m = hits.m,
#                                 info.features = info.features,
#                                 dat = dat,
#                                 ppm = 20)
# add.m = generate_connectivity_matrix(info.features = info.features,
#                                      type = "adducts")
#
# set.seed(20200402)
# post.prob = compute_posterior_prob(prior.prob = prior.prob,
#                                            dat = dat,
#                                            add.m = add.m,
#                                            delta.add = 0.1)
#
# info.assigned.use = assign_features(post.prob = post.prob,
#                                     dat = dat)
#
# se.example.met = summarize_metabolites(info.assigned = info.assigned.use,
#                                        info.features = info.features,
#                                        se = se.example)
#
