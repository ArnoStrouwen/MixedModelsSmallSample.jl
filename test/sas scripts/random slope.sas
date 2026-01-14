%web_drop_table(WORK.SS);
FILENAME REFFILE '/home/u64165441/sleepstudy.csv';

PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV
	OUT=WORK.SS;
	GETNAMES=YES;
RUN;

DATA WORK.SS;
    SET WORK.SS;
RUN;

/* Model 1: Zero Correlation (Independent Random Effects) */
ODS OUTPUT
    ParameterEstimates=PE_Model1
    AsyCov=AsyCov_Model1;

PROC GLIMMIX DATA=WORK.SS ASYCOV;
	CLASS subj;
    MODEL reaction = days / DDFM=Kenwardroger SOLUTION;
    RANDOM INTERCEPT / SUBJECT=subj;
    RANDOM days / SUBJECT=subj;
RUN;

ODS OUTPUT ParameterEstimates=PE_Model1_SW;
PROC GLIMMIX DATA=WORK.SS;
	CLASS subj;
    MODEL reaction = days / DDFM=Satterthwaite SOLUTION;
    RANDOM INTERCEPT / SUBJECT=subj;
    RANDOM days / SUBJECT=subj;
RUN;

/* Model 2: Correlated (UN) */
ODS OUTPUT
    ParameterEstimates=PE_Model2
    AsyCov=AsyCov_Model2;

PROC GLIMMIX DATA=WORK.SS ASYCOV;
	CLASS subj;
    MODEL reaction = days / DDFM=Kenwardroger SOLUTION;
    RANDOM INTERCEPT days / type=un SUBJECT=subj;
RUN;

ODS OUTPUT ParameterEstimates=PE_Model2_SW;
PROC GLIMMIX DATA=WORK.SS;
	CLASS subj;
    MODEL reaction = days / DDFM=Satterthwaite SOLUTION;
    RANDOM INTERCEPT days / type=un SUBJECT=subj;
RUN;

/* Model 3: Quadratic Zero Correlation */
ODS OUTPUT
    ParameterEstimates=PE_Model3
    AsyCov=AsyCov_Model3;

PROC GLIMMIX DATA=WORK.SS ASYCOV;
	CLASS subj;
    MODEL reaction = days days*days / DDFM=Kenwardroger SOLUTION;
    RANDOM INTERCEPT / SUBJECT=subj;
    RANDOM days / SUBJECT=subj;
    RANDOM days*days / SUBJECT=subj;
RUN;

ODS OUTPUT ParameterEstimates=PE_Model3_SW;
PROC GLIMMIX DATA=WORK.SS;
	CLASS subj;
    MODEL reaction = days days*days / DDFM=Satterthwaite SOLUTION;
    RANDOM INTERCEPT / SUBJECT=subj;
    RANDOM days / SUBJECT=subj;
    RANDOM days*days / SUBJECT=subj;
RUN;

/* Model 4: Quadratic Correlated (UN) */
ODS OUTPUT
    ParameterEstimates=PE_Model4
    AsyCov=AsyCov_Model4;

PROC GLIMMIX DATA=WORK.SS ASYCOV;
	CLASS subj;
    MODEL reaction = days days*days / DDFM=Kenwardroger SOLUTION;
    RANDOM INTERCEPT days days*days / type=un SUBJECT=subj;
RUN;

ODS OUTPUT ParameterEstimates=PE_Model4_SW;
PROC GLIMMIX DATA=WORK.SS;
	CLASS subj;
    MODEL reaction = days days*days / DDFM=Satterthwaite SOLUTION;
    RANDOM INTERCEPT days days*days / type=un SUBJECT=subj;
RUN;

/* Model 5: Quadratic Some Correlated (UN blocks) */
ODS OUTPUT
    ParameterEstimates=PE_Model5
    AsyCov=AsyCov_Model5;

PROC GLIMMIX DATA=WORK.SS ASYCOV;
	CLASS subj;
    MODEL reaction = days days*days / DDFM=Kenwardroger SOLUTION;
    RANDOM INTERCEPT / type=un SUBJECT=subj;
    RANDOM days days*days / type=un SUBJECT=subj;
RUN;

ODS OUTPUT ParameterEstimates=PE_Model5_SW;
PROC GLIMMIX DATA=WORK.SS;
	CLASS subj;
    MODEL reaction = days days*days / DDFM=Satterthwaite SOLUTION;
    RANDOM INTERCEPT / type=un SUBJECT=subj;
    RANDOM days days*days / type=un SUBJECT=subj;
RUN;

/* Model 6: FA0(2) - Factor Analytic Cholesky (2 terms, equivalent to full Cholesky for 2-term random effect) */
ODS OUTPUT
    ParameterEstimates=PE_Model6
    AsyCov=AsyCov_Model6;

PROC GLIMMIX DATA=WORK.SS ASYCOV;
	CLASS subj;
    MODEL reaction = days / DDFM=Kenwardroger SOLUTION;
    RANDOM INTERCEPT days / type=FA0(2) SUBJECT=subj;
RUN;

ODS OUTPUT ParameterEstimates=PE_Model6_SW;
PROC GLIMMIX DATA=WORK.SS;
	CLASS subj;
    MODEL reaction = days / DDFM=Satterthwaite SOLUTION;
    RANDOM INTERCEPT days / type=FA0(2) SUBJECT=subj;
RUN;

/* Apply Formatting */
DATA PE_Model1; SET PE_Model1; FORMAT _NUMERIC_ BEST32.; RUN;
DATA AsyCov_Model1; SET AsyCov_Model1; FORMAT _NUMERIC_ BEST32.; RUN;
DATA PE_Model1_SW; SET PE_Model1_SW; FORMAT _NUMERIC_ BEST32.; RUN;

DATA PE_Model2; SET PE_Model2; FORMAT _NUMERIC_ BEST32.; RUN;
DATA AsyCov_Model2; SET AsyCov_Model2; FORMAT _NUMERIC_ BEST32.; RUN;
DATA PE_Model2_SW; SET PE_Model2_SW; FORMAT _NUMERIC_ BEST32.; RUN;

DATA PE_Model3; SET PE_Model3; FORMAT _NUMERIC_ BEST32.; RUN;
DATA AsyCov_Model3; SET AsyCov_Model3; FORMAT _NUMERIC_ BEST32.; RUN;
DATA PE_Model3_SW; SET PE_Model3_SW; FORMAT _NUMERIC_ BEST32.; RUN;

DATA PE_Model4; SET PE_Model4; FORMAT _NUMERIC_ BEST32.; RUN;
DATA AsyCov_Model4; SET AsyCov_Model4; FORMAT _NUMERIC_ BEST32.; RUN;
DATA PE_Model4_SW; SET PE_Model4_SW; FORMAT _NUMERIC_ BEST32.; RUN;

DATA PE_Model5; SET PE_Model5; FORMAT _NUMERIC_ BEST32.; RUN;
DATA AsyCov_Model5; SET AsyCov_Model5; FORMAT _NUMERIC_ BEST32.; RUN;
DATA PE_Model5_SW; SET PE_Model5_SW; FORMAT _NUMERIC_ BEST32.; RUN;

DATA PE_Model6; SET PE_Model6; FORMAT _NUMERIC_ BEST32.; RUN;
DATA AsyCov_Model6; SET AsyCov_Model6; FORMAT _NUMERIC_ BEST32.; RUN;
DATA PE_Model6_SW; SET PE_Model6_SW; FORMAT _NUMERIC_ BEST32.; RUN;

/* Export Results */

/* Model 1 */
PROC EXPORT DATA=PE_Model1 OUTFILE='/home/u64165441/Results sleep study sas kr.csv' DBMS=CSV REPLACE; RUN;
PROC EXPORT DATA=AsyCov_Model1 OUTFILE='/home/u64165441/Results sleep study sas asycov.csv' DBMS=CSV REPLACE; RUN;
PROC EXPORT DATA=PE_Model1_SW OUTFILE='/home/u64165441/Results sleep study sas sw.csv' DBMS=CSV REPLACE; RUN;

/* Model 2 */
PROC EXPORT DATA=PE_Model2 OUTFILE='/home/u64165441/Results sleep study corr sas kr.csv' DBMS=CSV REPLACE; RUN;
PROC EXPORT DATA=AsyCov_Model2 OUTFILE='/home/u64165441/Results sleep study corr sas asycov.csv' DBMS=CSV REPLACE; RUN;
PROC EXPORT DATA=PE_Model2_SW OUTFILE='/home/u64165441/Results sleep study corr sas sw.csv' DBMS=CSV REPLACE; RUN;

/* Model 3 */
PROC EXPORT DATA=PE_Model3 OUTFILE='/home/u64165441/Results sleep study quadratic sas kr.csv' DBMS=CSV REPLACE; RUN;
PROC EXPORT DATA=AsyCov_Model3 OUTFILE='/home/u64165441/Results sleep study quadratic sas asycov.csv' DBMS=CSV REPLACE; RUN;
PROC EXPORT DATA=PE_Model3_SW OUTFILE='/home/u64165441/Results sleep study quadratic sas sw.csv' DBMS=CSV REPLACE; RUN;

/* Model 4 */
PROC EXPORT DATA=PE_Model4 OUTFILE='/home/u64165441/Results sleep study corr quadratic sas kr.csv' DBMS=CSV REPLACE; RUN;
PROC EXPORT DATA=AsyCov_Model4 OUTFILE='/home/u64165441/Results sleep study corr quadratic sas asycov.csv' DBMS=CSV REPLACE; RUN;
PROC EXPORT DATA=PE_Model4_SW OUTFILE='/home/u64165441/Results sleep study corr quadratic sas sw.csv' DBMS=CSV REPLACE; RUN;

/* Model 5 */
PROC EXPORT DATA=PE_Model5 OUTFILE='/home/u64165441/Results sleep study some corr quadratic sas kr.csv' DBMS=CSV REPLACE; RUN;
PROC EXPORT DATA=AsyCov_Model5 OUTFILE='/home/u64165441/Results sleep study some corr quadratic sas asycov.csv' DBMS=CSV REPLACE; RUN;
PROC EXPORT DATA=PE_Model5_SW OUTFILE='/home/u64165441/Results sleep study some corr quadratic sas sw.csv' DBMS=CSV REPLACE; RUN;

/* Model 6: FA0(2) */
PROC EXPORT DATA=PE_Model6 OUTFILE='/home/u64165441/Results sleep study fa0 sas kr.csv' DBMS=CSV REPLACE; RUN;
PROC EXPORT DATA=AsyCov_Model6 OUTFILE='/home/u64165441/Results sleep study fa0 sas asycov.csv' DBMS=CSV REPLACE; RUN;
PROC EXPORT DATA=PE_Model6_SW OUTFILE='/home/u64165441/Results sleep study fa0 sas sw.csv' DBMS=CSV REPLACE; RUN;
