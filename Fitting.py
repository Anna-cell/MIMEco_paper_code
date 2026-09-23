# Script to adapt resource to mimeco/CarveMe formats
import pandas as pd
import cobra
import mimeco

###################
### Condition 1 ###
###################

#Adapt YPD medium to CarveMe gapfilling format
YPD = pd.read_csv('media/YPD_media_BiGG.csv', sep=',', header=0, index_col=0)
medium = ["YPD"]*len(YPD)
compound = list(YPD.index)

YPD_CarveMe = pd.DataFrame({"medium":medium, "compound" :compound})
YPD_CarveMe.to_csv("media/YPD_gapfill.tsv", sep="\t", index=False)

###################
### Condition 4 ###
###################

# Medium LB+Cam broth.
# LB broth from CarveMe and add chloramphenicol

media_db = pd.read_csv('media/media_db.tsv', sep='\t', header=0)

met = []
for i in media_db.index:
    if media_db.loc[i, "medium"] == 'M9':
        met.append(media_db.loc[i, "compound"])
met.append("cm")
met.append("thm")
met.append("btn")

flux = [10]*len(met) #In machado's gapfill, max uptake for these coumpounds is 10.

min_medium = pd.DataFrame(flux, index=met, columns=["flux"])
print(min_medium)

min_medium.to_csv('media/min_medium.csv', index=True)

#e_coli deleted models

e_coli = cobra.io.read_sbml_model('models/iJO1366.xml')

with e_coli:
    e_coli.reactions.get_by_id('THRD_L').bounds = (0, 0)
    cobra.io.write_sbml_model(e_coli, "e_coli_noIle.xml")

with e_coli:
    e_coli.reactions.get_by_id('DAPDC').bounds = (0, 0)
    cobra.io.write_sbml_model(e_coli, "e_coli_noLys.xml")

#Smetana

# Using partially constrained medium worked with mimeco but it will not be possible to build the same condition in smetana. 
# I need a medium where metabolites are only described in presence/absence.
# So I will look for minimal medium with YPD + other metabolites at lox flux.
min_medium_untreated = pd.read_csv('media/min_medium.csv', sep=',', header=0, index_col=0)
print(min_medium_untreated)

e_coli_min_med, constrained_medium_dict = mimeco.utils.restrain_medium(e_coli, min_medium_untreated, 
                                                                       undescribed_metabolites_constraint = 'blocked', 
                                                                       undescribed_met_lb=-0.0001)

solo_growth_coli = e_coli_min_med.optimize().objective_value
print(solo_growth_coli)

NoIle = cobra.io.read_sbml_model('models/e_coli_noIle.xml')
NoIle.solver = "gurobi"
NoIle.id = "e.coli_NoIle"

NoIle_min_med, constrained_medium_dict = mimeco.utils.restrain_medium(NoIle, min_medium_untreated, 
                                                                       undescribed_metabolites_constraint = 'partially_constrained', 
                                                                       undescribed_met_lb=-0.01)
solo_growth_coli = NoIle_min_med.optimize().objective_value
print(solo_growth_coli)
print(constrained_medium_dict)

minimal_medium_cobra = cobra.medium.minimal_medium(NoIle_min_med, min_objective_value = 0.034, minimize_components=True)
print(minimal_medium_cobra)