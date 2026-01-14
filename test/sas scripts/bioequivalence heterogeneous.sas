/* Heterogeneous Bioequivalence Testing - MIXED vs GLIMMIX Comparison */
/* FDA Guidance model with heterogeneous residual variances */

DATA WORK.BE;
    INFILE '/home/u64165441/Data bioequivalence.csv' DLM=',' FIRSTOBS=2;
    INPUT subject formulation $ period sequence $ log_data;
RUN;

/* ============================================================== */
/* Section 1: FDA Guidance PROC MIXED - Heterogeneous Variance    */
/* ============================================================== */
/* Model Structure:
   - Random treatment effects by subject (FA0(2))
   - Heterogeneous residual variance by treatment group
   - NO pooled residual variance
*/

ODS OUTPUT
    Tests3=MIXED_Tests3
    SolutionF=MIXED_SolutionF
    CovParms=MIXED_CovParms
    Estimates=MIXED_Estimates;

PROC MIXED DATA=WORK.BE ASYCOV;
    CLASSES sequence subject period formulation;
    MODEL log_data = sequence period formulation / DDFM=SATTERTHWAITE SOLUTION;
    RANDOM formulation / TYPE=FA0(2) SUB=subject G;
    REPEATED / GRP=formulation SUB=subject;
    ESTIMATE 'T vs. R' formulation 1 -1 / CL ALPHA=0.1;
    ODS OUTPUT AsyCov=MIXED_AsyCov;
    TITLE1 'FDA Guidance Heterogeneous Model - PROC MIXED';
RUN;

/* Format and Export MIXED Results */
DATA MIXED_Tests3; SET MIXED_Tests3; FORMAT _NUMERIC_ BEST32.; RUN;
DATA MIXED_SolutionF; SET MIXED_SolutionF; FORMAT _NUMERIC_ BEST32.; RUN;
DATA MIXED_CovParms; SET MIXED_CovParms; FORMAT _NUMERIC_ BEST32.; RUN;
DATA MIXED_Estimates; SET MIXED_Estimates; FORMAT _NUMERIC_ BEST32.; RUN;
DATA MIXED_AsyCov; SET MIXED_AsyCov; FORMAT _NUMERIC_ BEST32.; RUN;

PROC EXPORT DATA=MIXED_Tests3 
    OUTFILE='/home/u64165441/Results bioequivalence heterogeneous mixed tests3.csv' 
    DBMS=CSV REPLACE; 
RUN;

PROC EXPORT DATA=MIXED_SolutionF 
    OUTFILE='/home/u64165441/Results bioequivalence heterogeneous mixed solution.csv' 
    DBMS=CSV REPLACE; 
RUN;

PROC EXPORT DATA=MIXED_CovParms 
    OUTFILE='/home/u64165441/Results bioequivalence heterogeneous mixed covparms.csv' 
    DBMS=CSV REPLACE; 
RUN;

PROC EXPORT DATA=MIXED_Estimates 
    OUTFILE='/home/u64165441/Results bioequivalence heterogeneous mixed estimates.csv' 
    DBMS=CSV REPLACE; 
RUN;

PROC EXPORT DATA=MIXED_AsyCov 
    OUTFILE='/home/u64165441/Results bioequivalence heterogeneous mixed asycov.csv' 
    DBMS=CSV REPLACE; 
RUN;

/* ============================================================== */
/* Section 2: PROC GLIMMIX - Equivalent Heterogeneous Model       */
/* ============================================================== */
/* Using GLIMMIX syntax to achieve same structure */

ODS OUTPUT
    Tests3=GLIMMIX_Tests3
    ParameterEstimates=GLIMMIX_PE
    CovParms=GLIMMIX_CovParms
    Estimates=GLIMMIX_Estimates;

PROC GLIMMIX DATA=WORK.BE ASYCOV;
    CLASSES sequence subject period formulation;
    MODEL log_data = sequence period formulation / DDFM=SATTERTHWAITE SOLUTION;
    RANDOM formulation / TYPE=FA0(2) SUB=subject G;
    RANDOM _RESIDUAL_ / GROUP=formulation;
    ESTIMATE 'T vs. R' formulation 1 -1 / CL ALPHA=0.1;
    ODS OUTPUT AsyCov=GLIMMIX_AsyCov;
    TITLE1 'FDA Guidance Heterogeneous Model - PROC GLIMMIX';
RUN;

/* Format and Export GLIMMIX Results */
DATA GLIMMIX_Tests3; SET GLIMMIX_Tests3; FORMAT _NUMERIC_ BEST32.; RUN;
DATA GLIMMIX_PE; SET GLIMMIX_PE; FORMAT _NUMERIC_ BEST32.; RUN;
DATA GLIMMIX_CovParms; SET GLIMMIX_CovParms; FORMAT _NUMERIC_ BEST32.; RUN;
DATA GLIMMIX_Estimates; SET GLIMMIX_Estimates; FORMAT _NUMERIC_ BEST32.; RUN;
DATA GLIMMIX_AsyCov; SET GLIMMIX_AsyCov; FORMAT _NUMERIC_ BEST32.; RUN;

PROC EXPORT DATA=GLIMMIX_Tests3 
    OUTFILE='/home/u64165441/Results bioequivalence heterogeneous glimmix tests3.csv' 
    DBMS=CSV REPLACE; 
RUN;

PROC EXPORT DATA=GLIMMIX_PE 
    OUTFILE='/home/u64165441/Results bioequivalence heterogeneous glimmix pe.csv' 
    DBMS=CSV REPLACE; 
RUN;

PROC EXPORT DATA=GLIMMIX_CovParms 
    OUTFILE='/home/u64165441/Results bioequivalence heterogeneous glimmix covparms.csv' 
    DBMS=CSV REPLACE; 
RUN;

PROC EXPORT DATA=GLIMMIX_Estimates 
    OUTFILE='/home/u64165441/Results bioequivalence heterogeneous glimmix estimates.csv' 
    DBMS=CSV REPLACE; 
RUN;

PROC EXPORT DATA=GLIMMIX_AsyCov 
    OUTFILE='/home/u64165441/Results bioequivalence heterogeneous glimmix asycov.csv' 
    DBMS=CSV REPLACE; 
RUN;

