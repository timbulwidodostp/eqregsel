{smcl}
{* *! version 1.2.1  25january2017}{...}
{viewerjumpto "Title" "eqregsel##title"}{...}
{viewerjumpto "Syntax" "eqregsel##syntax"}{...}
{viewerjumpto "Description" "eqregsel##description"}{...}
{viewerjumpto "Options" "eqregsel##options"}{...}
{viewerjumpto "Examples" "eqregsel##examples"}{...}
{viewerjumpto "Stored results" "eqregsel##storedresults"}{...}
{viewerjumpto "Version requirements" "eqregsel##version"}{...}
{viewerjumpto "Methods and formulas" "eqregsel##methods"}{...}
{viewerjumpto "Supporting package" "eqregsel##supporting"}{...}
{viewerjumpto "References" "eqregsel##references"}{...}
{viewerjumpto "Authors" "eqregsel##authors"}{...}
{viewerjumpto "Also see" "eqregsel##alsosee"}{...}
{cmd:help eqregsel}{right: ({browse "https://doi.org/10.1177/1536867X20930998":SJ20-2: st0598})}
{hline}

{marker title}{...}
{title:Title}

{p2colset 5 17 19 2}{...}
{p2col :{cmd:eqregsel} {hline 2}}Estimation method for sample-selection models
based on extremal quantile regressions{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:eqregsel}
{it:Y} {it:X}1 {it:X}2
{ifin}
[{cmd:,} {it:options}]

{synoptset 10}{...}
{synopthdr}
{synoptline}
{synopt:{opt hom(#)}}estimate the coefficients of the first {it:#} independent
variables; default is {cmd:hom(1)}{p_end}
{synopt:{opt subs(#)}}subsample size used for selecting quantile index (tau);
default is computed as a function of the size of the dataset{p_end}
{synopt:{opt grid(#)}}discretize the interval of tau into {it:#} subintervals;
default is {cmd:grid(40)}{p_end}
{synopt:{opt rep(#)}}perform {it:#} bootstrap and subsampling replications; default
is {cmd:rep(150)}{p_end}
{synopt:{opt small}}estimate the coefficients using the opposite values of all independent variables{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:eqregsel} estimates and provides inference for a semiparametric
sample-selection model.  This estimation method, which is based on quantile
regressions in the tails of the outcome (extremal quantile regression), can be
used in the absence of an instrument or a large support regressor.  The
estimation and inference methods are based on D'Haultf{c 0x9c}uille, Maurel,
and Zhang (2018).  The variables of interest are assumed to have homogeneous
effects on the outcome across its distribution.


{marker options}{...}
{title:Options}

{phang}
{opt hom(#)} specifies that the first {it:#} variables after the dependent
variable have homogeneous effects on the outcome across its distribution.  The
code then returns their estimated effects and standard errors.  The default is
{cmd:hom(1)}.

{phang}
{opt subs(#)} specifies the subsample size used when selecting the quantile
index tau following the data-driven procedure proposed in 
D'Haultf{c 0x9c}uille, Maurel, and Zhang (2018).  The quantile index tau is
computed as the minimizer of a criterion that captures the tradeoff between
bias and variance.  This criterion is computed with a subsample, whose size is
specified by {cmd:subs()}.  If no value (or a negative value) is specified, the
subsample size is calculated using the formula given in section 3.1 in
D'Haultf{c 0x9c}uille, Maurel, and Zhang (2018).

{phang}
{opt grid(#)} specifies the number of (equally spaced) grid points in the
admissible interval for searching for the quantile index tau.  The upper and
lower bound of the admissible interval are 0.3 and min(0.1, 80/b_n),
respectively, where b_n is specified by {cmd:subs()}.  The default is
{cmd:grid(40)}.

{phang}
{opt rep(#)} specifies the number of bootstrap and subsampling replications.
Subsampling is used when selecting tau, while bootstrap is used to compute the
standard errors of the estimators.  The default is {cmd:rep(150)}.

{phang}
{cmd:small} specifies that the opposite values of all independent variables
are used for the estimation.


{marker examples}{...}
{title:Examples}

{phang2}{cmd:. use bw_nlsy7997}{p_end}
{phang2}{cmd:. keep if cohort79==1}{p_end}
{phang2}{cmd:. generate afqt2=afqt^2}{p_end}

{pstd}
{cmd:black} is the only variable of interest{p_end}
{phang2}{cmd:. eqregsel log_wage black hispanic age afqt afqt2}

{pstd}
{cmd:black} and {cmd:hispanic} are the variables of interest{p_end}
{phang2}{cmd:. eqregsel log_wage black hispanic age afqt afqt2, hom(2)}

{pstd}
{cmd:black} is the only variable of interest; 500 bootstrap and subsampling
replications are computed{p_end}
{phang2}{cmd:. eqregsel log_wage black hispanic age afqt afqt2, rep(500)}


{marker storedresults}{...}
{title:Stored results}

{pstd}
{cmd:eqregsel} stores the following in {cmd:e()}:

{synoptset 22 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(tau0)}}quantile index{p_end}
{synopt:{cmd:e(specificationtest)}}p-value of the specification test{p_end}
{synopt:{cmd:e(subs)}}subsample size when selecting the quantile index{p_end}
{synopt:{cmd:e(homvar)}}number of variables with homogeneous effects on the outcome{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(beta_hom)}}a d_1 x 1 matrix containing the estimated
coefficient or coefficients of interest{p_end}
{synopt:{cmd:e(std_b)}}a d_1 x 1 matrix containing the standard error
of the estimator or estimators{p_end}


{marker version}{...}
{title:Version requirements}

{p 4 4 2}
This command requires Stata 12 or later.


{marker methods}{...}
{title:Methods and formulas}

{pstd}
See D'Haultf{c 0x9c}uille, Maurel, and Zhang (2018).


{marker supporting}{...}
{title:Supporting package}

{p 4 4 2}
This command requires the {cmd:moremata} package (Jann 2005).  Before calling
this command, please install {cmd:moremata} by typing

{phang2}
{cmd:. ssc install moremata}{p_end}


{marker references}{...}
{title:References}

{phang}
D'Haultf{c 0x9c}uille, X., A. Maurel, and Y. Zhang. 2018. Extremal quantile
regressions for selection models and the black-white wage gap.  
{it:Journal of Econometrics} 203: 129-142.
{browse "https://doi.org/10.1016/j.jeconom.2017.11.004"}.

{phang}
Jann, B. 2005. moremata: Stata module (Mata) to provide various functions.
Statistical Software Components S455001, Department of Economics, Boston
College.  {browse "https://ideas.repec.org/c/boc/bocode/s455001.html"}.


{marker authors}{...}
{title:Authors}

{pstd}
Xavier D'Haultf{c 0x9c}uille{break}
CREST-ENSAE{break}
Paris, France{break}
xavier.dhaultfoeuille@ensae.fr

{pstd}
Arnaud Maurel{break}
Duke University, NBER, and IZA{break}
Durham, NC{break}
apm16@duke.edu

{pstd}
Xiaoyun Qiu{break}
Northwestern University{break}
Evanston, IL{break}
xiaoyun.qiu@u.northwestern.edu

{pstd}
Yichong Zhang{break}
Singapore Management University{break}
Singapore, Singapore{break}
yczhang@smu.edu.sg


{marker alsosee}{...}
{title:Also see}

{p 4 14 2}
Article:  {it:Stata Journal}, volume 20, number 2: {browse "https://doi.org/10.1177/1536867X20930998":st0598}{p_end}
