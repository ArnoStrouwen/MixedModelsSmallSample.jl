DATA WORK.WT;
    INFILE '/home/u64165441/Data bioequivalence.csv' DLM=',' FIRSTOBS=2;
    INPUT subject formulation $ period sequence $ log_data;
RUN;

ODS OUTPUT
    Tests3=Tests3_SW_ML
    ParameterEstimates=PE_SW_ML
    AsyCov=AsyCov_ML;

/*
GLIMMIX does not accept METHOD=ML; for Gaussian LMMs:
  - METHOD=MSPL corresponds to ML
  - METHOD=RSPL corresponds to REML
*/
PROC GLIMMIX DATA=WORK.WT ASYCOV METHOD=MSPL;
    CLASS subject sequence period formulation;
    MODEL log_data = period formulation sequence / DDFM=SATTERTHWAITE SOLUTION;
    RANDOM INTERCEPT / SUBJECT=subject;
RUN;

DATA PE_SW_ML;
    SET PE_SW_ML;
    FORMAT _NUMERIC_ BEST32.;
RUN;

PROC EXPORT DATA=PE_SW_ML
    OUTFILE='/home/u64165441/Results bioequivalence homogeneous sas sw ml.csv'
    DBMS=CSV REPLACE;
RUN;

DATA Tests3_SW_ML;
    SET Tests3_SW_ML;
    FORMAT _NUMERIC_ BEST32.;
RUN;

PROC EXPORT DATA=Tests3_SW_ML
    OUTFILE='/home/u64165441/Results bioequivalence homogeneous sas tests3 sw ml.csv'
    DBMS=CSV REPLACE;
RUN;

DATA AsyCov_ML;
    SET AsyCov_ML;
    FORMAT _NUMERIC_ BEST32.;
RUN;

PROC EXPORT DATA=AsyCov_ML
    OUTFILE='/home/u64165441/Results bioequivalence homogeneous sas asycov ml.csv'
    DBMS=CSV REPLACE;
RUN;
