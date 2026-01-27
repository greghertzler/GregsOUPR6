# Open file in RStudio and choose Source with Echo, or
# Type 'source("demo/Adoption_EntryOption.R",echo=TRUE)', or
# Type 'demo(Adoption_EntryOption)'
# R6 object
A <- Analytical$new()
# option envelope
A$OptionEnvelope(phi=1)
# entry decision
A$DecisionThreshold()
# break even
A$DecisionThreshold(y=5)
# entry subsidy
A$DecisionThreshold(b=10)
# first passage time from exit to entry
A$DecisionThreshold(y=0,b=0)
A$sync_zyxt_stoch()
A$PassageTimePercentiles()
# less uncertainty
A$DecisionThreshold(sigma=5)
A$sync_zyxt_stoch()
A$PassageTimePercentiles()
