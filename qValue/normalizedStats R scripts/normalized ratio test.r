# Reads in test data sets, construct linear regression models, and do
# ANOVA to compare two models and determine if KOnp/WTnp == KOpp/WTpp

# Part 1: unpaired (label-free) data test

setwd('C:/Users/knoh1/Documents/qValue_ken/normalizedStats R scripts/artificial test data (unpaired)')

s = read.table('sampleData_S.txt',sep='\t',header=TRUE) # p < 0.05
ns = read.table('sampleData_NS.txt',sep='\t',header=TRUE) # n.s.
s2 = read.table('sampleData_S2.txt',sep='\t',header=TRUE) # p < 0.05
ns2 = read.table('sampleData_NS2.txt',sep='\t',header=TRUE) # n.s.

dat = ns2

# construct models with and without interaction
model1 = lm(log(peakArea) ~ ko + np, data = dat, na.action = na.omit)
model2 = lm(log(peakArea) ~ ko*np, data = dat, na.action = na.omit) # ko*np == ko + np + ko:np

# compare models, view results, and get p value for comparison
a = anova(model1, model2)
a
a[['Pr(>F)']][2]

# Part 2: paired (SILAC) data test
# data from FM still comes in as peak areas from two separate channels,
# so you will need to construct the paired ratios (e.g. [Ch1.rep2.timepoint4]/[Ch2.rep2.timepoint4])
# before testing. The test data examples below are provided as paired ratios.

setwd('C:/Users/knoh1/Documents/qValue_ken/normalizedStats R scripts/artificial test data (paired)')

s = read.table('sampleData_S.txt',sep='\t',header=TRUE) # p < 0.05
ns = read.table('sampleData_NS.txt',sep='\t',header=TRUE) # n.s.
s2 = read.table('sampleData_S2.txt',sep='\t',header=TRUE) # p < 0.05
ns2 = read.table('sampleData_NS2.txt',sep='\t',header=TRUE) # n.s.

dat = ns2

# construct models with and without interaction
model = lm(log(ratio) ~ np, data = dat, na.action = na.omit)

# compare models, view results, and get p value for comparison
a = anova(model)
a
a[['Pr(>F)']][1]
