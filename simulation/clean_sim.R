# Removes cases with a missing or corrupted model file
exp = 'finalc'
sim_init()

caseData = getExperimentResults(exp)
expFiles = vapply(caseData, '[[', 'output', FUN.VALUE = '')
validMsk = file.exists(expFiles)
print(sum(!validMsk))

deleteCases(exp = exp, sim = FALSE, pattern = names(expFiles)[!validMsk], fixed = TRUE)
