/* Intermediate Bioequivalence Model
   Heterogeneous random effects (FA0(2)) + Homogeneous residual variance
   This is a middle ground between:
   - Homogeneous: RANDOM INTERCEPT / SUBJECT=subject
   - Heterogeneous: RANDOM formulation / TYPE=FA0(2) + REPEATED / GRP=formulation
*/

DATA WORK.BE;
    INFILE '/home/u64165441/Data bioequivalence.csv' DLM=',' FIRSTOBS=2;
    INPUT subject formulation $ period sequence $ log_data;
RUN;

/* Model with FA0(2) random effects and common residual - Kenward-Roger */
ODS OUTPUT
    Tests3=Tests3_KR
    ParameterEstimates=PE_KR
    AsyCov=AsyCov_KR
    CovParms=CovParms_KR;

PROC GLIMMIX DATA=WORK.BE ASYCOV;
    CLASS subject sequence period formulation;
    MODEL log_data = period formulation sequence / DDFM=KENWARDROGER SOLUTION COVB(DETAILS);
    RANDOM formulation / TYPE=FA0(2) SUBJECT=subject;
RUN;

DATA PE_KR; SET PE_KR; FORMAT _NUMERIC_ BEST32.; RUN;
DATA Tests3_KR; SET Tests3_KR; FORMAT _NUMERIC_ BEST32.; RUN;
DATA AsyCov_KR; SET AsyCov_KR; FORMAT _NUMERIC_ BEST32.; RUN;
DATA CovParms_KR; SET CovParms_KR; FORMAT _NUMERIC_ BEST32.; RUN;

PROC EXPORT DATA=PE_KR OUTFILE='/home/u64165441/Results bioequivalence intermediate sas kr.csv' DBMS=CSV REPLACE; RUN;
PROC EXPORT DATA=Tests3_KR OUTFILE='/home/u64165441/Results bioequivalence intermediate sas tests3 kr.csv' DBMS=CSV REPLACE; RUN;
PROC EXPORT DATA=AsyCov_KR OUTFILE='/home/u64165441/Results bioequivalence intermediate sas asycov.csv' DBMS=CSV REPLACE; RUN;
PROC EXPORT DATA=CovParms_KR OUTFILE='/home/u64165441/Results bioequivalence intermediate sas covparms.csv' DBMS=CSV REPLACE; RUN;

/* Same model - Satterthwaite */
ODS OUTPUT
    Tests3=Tests3_SW
    ParameterEstimates=PE_SW;

PROC GLIMMIX DATA=WORK.BE;
    CLASS subject sequence period formulation;
    MODEL log_data = period formulation sequence / DDFM=SATTERTHWAITE SOLUTION;
    RANDOM formulation / TYPE=FA0(2) SUBJECT=subject;
RUN;

DATA PE_SW; SET PE_SW; FORMAT _NUMERIC_ BEST32.; RUN;
DATA Tests3_SW; SET Tests3_SW; FORMAT _NUMERIC_ BEST32.; RUN;

PROC EXPORT DATA=PE_SW OUTFILE='/home/u64165441/Results bioequivalence intermediate sas sw.csv' DBMS=CSV REPLACE; RUN;
PROC EXPORT DATA=Tests3_SW OUTFILE='/home/u64165441/Results bioequivalence intermediate sas tests3 sw.csv' DBMS=CSV REPLACE; RUN;
