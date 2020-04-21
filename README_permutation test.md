# Permutation test-Python
"<font color=Red>Permutation test</font>" is a method to investigate the relationship between average frequency of disease-specific genotype combinations which made up by only 3 SNPs (fi) and the prevalence of disease in different ethnicities (P).<br> <br> Three different models were simulated by the process developed by ourselvese with a script written in python language called the complete random model, constraint model, constraint lost and mixed model respectively.
## Complete random model
#python completely_random.py.txt 100 20 # 100 times permutation, 20pair data
## Constraint model
#python constraint.py.txt 100 20 20# 100 times permutation, 20pair data, 20 polarization
## Constraint lost and mixed model
#python constraint_lost_mixed.py.txt 100 20 20 10 5# 100 times permutation, 20pair data, 20 polarization,10 lossing proportion, 5 mixing  proportion
