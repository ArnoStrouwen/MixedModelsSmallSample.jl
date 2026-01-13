DATA WORK.WT;
    INFILE '/home/u64165441/Data bioequivalence.csv' DLM=',' FIRSTOBS=2;
    INPUT subject formulation $ period sequence $ log_data;
RUN;

ODS OUTPUT
    Tests3=Tests3_KR
    ParameterEstimates=PE_KR
    AsyCov=AsyCov;

PROC GLIMMIX DATA=WORK.WT ASYCOV;
    CLASS subject sequence period formulation;
    MODEL log_data = period formulation sequence / DDFM=KENWARDROGER SOLUTION;
    RANDOM INTERCEPT / SUBJECT=subject;
RUN;

DATA PE_KR;
    SET PE_KR;
    FORMAT _NUMERIC_ BEST32.;
RUN;

PROC EXPORT DATA=PE_KR
    OUTFILE='/home/u64165441/Results bioequivalence homogeneous sas kr.csv'
    DBMS=CSV REPLACE;
RUN;

DATA Tests3_KR;
    SET Tests3_KR;
    FORMAT _NUMERIC_ BEST32.;
RUN;

PROC EXPORT DATA=Tests3_KR
    OUTFILE='/home/u64165441/Results bioequivalence homogeneous sas tests3 kr.csv'
    DBMS=CSV REPLACE;
RUN;

DATA AsyCov;
    SET AsyCov;
    FORMAT _NUMERIC_ BEST32.;
RUN;

PROC EXPORT DATA=AsyCov
    OUTFILE='/home/u64165441/Results bioequivalence homogeneous sas asycov.csv'
    DBMS=CSV REPLACE;
RUN;

ODS OUTPUT
    Tests3=Tests3_SW
    ParameterEstimates=PE_SW;

PROC GLIMMIX DATA=WORK.WT;
    CLASS subject sequence period formulation;
    MODEL log_data = period formulation sequence / DDFM=SATTERTHWAITE SOLUTION;
    RANDOM INTERCEPT / SUBJECT=subject;
RUN;

DATA PE_SW;
    SET PE_SW;
    FORMAT _NUMERIC_ BEST32.;
RUN;

PROC EXPORT DATA=PE_SW
    OUTFILE='/home/u64165441/Results bioequivalence homogeneous sas sw.csv'
    DBMS=CSV REPLACE;
RUN;

DATA Tests3_SW;
    SET Tests3_SW;
    FORMAT _NUMERIC_ BEST32.;
RUN;

PROC EXPORT DATA=Tests3_SW
    OUTFILE='/home/u64165441/Results bioequivalence homogeneous sas tests3 sw.csv'
    DBMS=CSV REPLACE;
RUN;


