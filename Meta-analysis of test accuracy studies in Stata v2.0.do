
/********************************************************************************************
*																							*
* TITLE:			Meta-analysis of diagnostic test accuracy studies						*
* AUTHOR:			Yemisi Takwoingi														*
* DATE CREATED:		01/10/2013								   								*
* DATE MODIFIED:	25/08/2020								   								*
* PURPOSE: 			Bivariate model for meta-analysis of CT and MRI for diagnosis of CAD	*
*									 														*
* ACKNOWLEDGEMENT:	Some of the code in here has been informed by Roger Harbord's metandi 	*
* 					which performs meta-analysis of a single test without options for 		*
*					meta-regression 														*
*																							*
********************************************************************************************/


/* Set your working directory to the appropriate drive where you saved the file schuetz.csv. 
Replace "U:\Handbook 2020" with the appropriate path for you. */
cd "U:\Handbook 2020"

**************************** 1. DATA IMPORT ****************************************************************

*** Read in the data from the .csv file ***
insheet using "schuetz.csv", comma clear

*** Produce a summary of the dataset to check data import was ok ***
describe 

*** Convert the string variable 'test' to a numeric variable named 'testtype' ***
encode test, gen(testtype) 

*** List the numeric value assigned to each test ***
label list testtype

**********************2. META-ANALYSIS OF A SINGLE TEST - METANDI ******************************************

*** Install gllamm and metandi if needed (to run each command statement remove the *) ***
* ssc install gllamm, replace	
* ssc install metandi, replace

*** Use metandi to meta-analyse studies that evaluated CT
metandi tp fp fn tn if testtype==1

*** To obtain a SROC plot as well as parameter estimates, add the plot option to the metandi statement as follows ***
metandi tp fp fn tn if testtype==1, plot

*** There is no option with metandi to modify the plot but this can be done using metandiplot. ***
metandiplot if testtype==1

/* If the optional variables tp fp fn tn are included in the command line, the plot also includes estimates 
of sensitivity and specificity from each study. */
metandiplot tp fp fn tn if testtype==1

/* The default is to scale the plot symbol by the sample size of each study. To make the symbols all the same size, 
specify constant weights, e.g. [aw=1]. Try some other options too. 
The command below will produce a plot using constant weights for the plot symbol, remove the confidence region and SROC curve,
as well as draw a 50% prediction region around the summary point. */
metandiplot tp fp fn tn if testtype==1 [aw=1], conf(off) curve(off) predlevel(50)

*** Here’s another example including some twoway graph options. ***
metandiplot tp fp fn tn  [aw=1], curve(off) legend(off) title(SROC plot of CT for CAD) scheme(s1mono)


**********************3. META-ANALYSIS OF A SINGLE TEST - meqrlogit ******************************************

/*Set up the data
Generate 5 new variables of type long. We need these before we can reshape the data.
•	n1 is number diseased
•	n0 is number without disease
•	true1 is number of true positives
•	true0 is the number of true negatives
•	recordid is the unique identifier for each study (and test if a study evaluated more than one test). _n will generate a sequence of numbers.*/

gen long n1=tp+fn
gen long n0=fp+tn
gen long true1=tp
gen long true0=tn
gen long recordid= _n

*** Reshape the data from wide to long format ***
reshape long n true, i(recordid) j(sens)

*** Sort data by study and sens variables to cluster the 2 records per study together ***
sort study sens
gen byte spec=1-sens

*** Perform meta-analysis for CT ***
meqrlogit true sens spec if testtype==1, nocons|| study_id: sens spec, /// 	
nocons cov(un) binomial(n) refineopts(iterate(3)) intpoints(5) variance 

/* To find the covariance between the expected (mean) logit sensitivity and expected logit specificity, 
display contents of the variance-covariance matrix: */
matrix list e(V)  

*** Display the coefficient vector ***
matrix list e(b)        
mmm
**********************4. DISPLAY SUMMARY ESTIMATES *******************************************************
*** Drop the program in case it is already in Stata’s memory. *** 
capture program drop renamematrix

program define renamematrix, eclass
matrix mb = e(b)
matrix mv = e(V)
matrix colnames mb = logitse:_cons logitsp:_cons 
matrix colnames mv = logitse:_cons logitsp:_cons 
matrix rownames mv = logitse:_cons logitsp:_cons 
ereturn post mb mv
end

renamematrix

_diparm logitse, label(Sensitivity) invlogit 
_diparm logitsp, label(Specificity) invlogit 
_diparm logitse logitsp, label(DOR) ci(log) function(exp(@1+@2)) derivative(exp(@1+@2) exp(@1+@2))
_diparm logitse logitsp, label(LR+) ci(log) function(invlogit(@1)/(1-invlogit(@2))) derivative(exp(@2-@1)*invlogit(@1)^2/invlogit(@2) exp(@2)*invlogit(@1))
_diparm logitse logitsp, label(LR-) ci(log) function((1-invlogit(@1))/invlogit(@2)) derivative(exp(-@1)*invlogit(@1)^2/invlogit(@2) exp(-@1-@2)*invlogit(@1))  


***IMPORTANT NOTE: If you need to run metandi again remember you have modified the dataset so use the original data.


********************** 5. META-REGRESSION  - USING MEQRLOGIT WITH A COVARIATE ***********************************************

*** Create dummy variables for the covariate testtype ***
gen seCT=0
gen spCT=0
gen seMRI=0
gen spMRI=0
replace seCT=1 if testtype==1 & sens==1
replace spCT=1 if testtype==1 & spec==1
replace seMRI=1 if testtype==2 & sens==1
replace spMRI=1 if testtype==2 & spec==1
														
*** Meta-analysis of CT ***
meqrlogit true sens spec if testtype==1, nocons || study_id: sens spec, nocons cov(un) binomial(n)  refineopts(iterate(3)) intpoints(5) variance 

*** Meta-analysis of MRI ***
meqrlogit true sens spec if testtype==2, nocons || study_id: sens spec, nocons cov(un) binomial(n)  refineopts(iterate(3)) intpoints(5) variance 

*** Fit the model without the covariate ***
meqrlogit true sens spec, nocons || study_id: sens spec, nocons cov(un) ///
binomial(n)  refineopts(iterate(3)) intpoints(5) variance nolr

*** Store the estimates of the log likelihood from the model above for doing the likeilhood ratio test later ***
estimates store A

*** Add covariate terms to the model for both logit sensitivity and logit specificity. This model assumes equal variances for both tests. ***
meqrlogit true seCT seMRI spCT spMRI, nocons || study_id: sens spec, nocons cov(un) binomial(n)  ///
refineopts(iterate(3)) intpoints(5) variance nolr

estimates store B

/* Perform a likelihood ratio test comparing the model (A) without covariate with the model (B) that includes the covariate testtype and 
assumes equal variances for each test. Use the stored values in A and B. */
lrtest A B

*** Assume sensitivity is the same for CT and MRI but allow specificity to vary with test type ***
meqrlogit true sens spCT spMRI, nocons || study_id: sens spec, nocons cov(un) binomial(n)  ///
refineopts(iterate(3)) intpoints(5) variance nolr

estimates store C

/* Perform a likelihood ratio test comparing the model (C) that assumes sensitivity is the same for CT and MRI but allows specificity to vary with test type 
with the model (B) that allows both sensitivity and specificity to vary with testtype. Use the stored values in B and C. */
lrtest B C 

*** Assume specificity is the same for CT and MRI but allow sensitivity to vary with test type ***
meqrlogit true seCT seMRI spec, nocons || study_id: sens spec, nocons cov(un) binomial(n)  ///
refineopts(iterate(3)) intpoints(5) variance nolr

estimates store D

/* Perform a likelihood ratio test comparing the model (D) that assumes specificity is the same for CT and MRI but allows sensitivity to vary with test type 
with the model (B) that allows both sensitivity and specificity to vary with testtype. Use the stored values in B and D. */
lrtest B D 

/* Having looked at the variances assumption of equal variances may not be appropriate.
There was a marked difference in the variances for the logit specificities in the 2 meta-analyses so
let's fit a model with separate variances for the logits for each test */

*** Fit model (E) with covariate testtype and separate variances for the tests ***
meqrlogit true seCT seMRI spCT spMRI, nocons || study_id: seCT spCT, nocons cov(un) ///
|| study_id seMRI spMRI, nocons cov(un) binomial(n)  refineopts(iterate(3)) intpoints(5) variance nolr

estimates store E

/* Perform a likelihood ratio test comparing the model (B) that assumes equal variances with the model (E) that allows for separate variances for 
each test. Use the stored values in B and E. */
lrtest B E

/* Perform a likelihood ratio test comparing the model (A) without covariate  with the model (E) that includes the covariate testtype 
and allows for separate variances for each test. Use the stored values in A and E. */
lrtest A E

/* To find the covariance between the estimated mean logit sensitivity and mean logit specificity for each test, 
display contents of the variance-covariance matrix. */
matrix list e(V)  

*** Delete the program from Stata's memory if it exists already ***
capture program drop renamematrix

*** Rename the elements of the coefficient and variance matrices ***
program define renamematrix, eclass
	matrix mb = e(b)
	matrix mv = e(V)
	matrix colnames mb = logitseCT:_cons logitseMRI:_cons logitspCT:_cons logitspMRI:_cons 
	matrix colnames mv = logitseCT:_cons logitseMRI:_cons logitspCT:_cons logitspMRI:_cons 
	matrix rownames mv = logitseCT:_cons logitseMRI:_cons logitspCT:_cons logitspMRI:_cons 
	ereturn post mb mv
end

*** Run the program ***
renamematrix
 
*** Display summary estimates by taking the inverse logits of the  mean logit sensitivity and mean logit specificity for each test ***
_diparm logitseCT, label(Sensitivity CT) invlogit 
_diparm logitseMRI, label(Sensitivity MRI) invlogit 
_diparm logitspCT, label(Specificity CT) invlogit 
_diparm logitspMRI, label(Specificity MRI) invlogit 

*** Display other summary estimates derived using functions of the mean logit sensitivities and mean logit specificities ***	
_diparm logitseCT logitspCT, label(LR+ CT) ci(log) function(invlogit(@1)/(1-invlogit(@2))) derivative(exp(@2-@1)*invlogit(@1)^2/invlogit(@2) exp(@2)*invlogit(@1))
_diparm logitseMRI logitspMRI, label(LR+ MRI) ci(log) function(invlogit(@1)/(1-invlogit(@2))) derivative(exp(@2-@1)*invlogit(@1)^2/invlogit(@2) exp(@2)*invlogit(@1))
_diparm logitseCT logitspCT, label(LR- CT) ci(log) function((1-invlogit(@1))/invlogit(@2)) derivative(exp(-@1)*invlogit(@1)^2/invlogit(@2)  exp(-@1-@2)*invlogit(@1))  
_diparm logitseMRI logitspMRI, label(LR- MRI) ci(log) function((1-invlogit(@1))/invlogit(@2)) derivative(exp(-@1)*invlogit(@1)^2/invlogit(@2)  exp(-@1-@2)*invlogit(@1)) 



***********************************************Additional analyses ***************************************************************************************

*** Fit model (F) assuming same sensitivity but separate variances for the tests ***
meqrlogit true sens spCT spMRI, nocons || study_id: seCT spCT, nocons cov(un) ///
|| study_id seMRI spMRI, nocons cov(un) binomial(n)  refineopts(iterate(3)) intpoints(5) variance nolr

estimates store F

/* Perform a likelihood ratio test comparing the model (F) that assumes the same sensitivity with the model (E) that allows for the effect of testtype. 
Use the stored values in F and E. */
lrtest F E

*** Fit model (G) assuming same specificity but separate variances for the tests ***
meqrlogit true seCT seMRI spec, nocons || study_id: seCT spCT, nocons cov(un) ///
|| study_id seMRI spMRI, nocons cov(un) binomial(n)  refineopts(iterate(3)) intpoints(5) variance nolr

estimates store G

/* Perform a likelihood ratio test comparing the model (G) that assumes the same specificity with the model (E) that allows for the effect of testtype. 
Use the stored values in G and E. */
lrtest G E

*** Meta-analysis with gllamm ***
eq eq1: sens
eq eq0: spec

gllamm true sens spec if testtype==1, nocons i(study_id) nrf(2) ///
eqs(eq1 eq0) family(binomial) link(logit) denom(n)  ip(g) nip(5) adapt



