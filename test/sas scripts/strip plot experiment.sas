%web_drop_table(WORK.BC);
FILENAME REFFILE '/home/u64165441/Data battery cell Chapter 11.csv';

PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV
	OUT=WORK.BC;
	GETNAMES=YES;
RUN;

DATA WORK.BC;
    SET WORK.BC;
    RENAME "Whole Plots"n = WP;
    RENAME "SubPlots"n = SP;
RUN;

ODS OUTPUT
    ParameterEstimates=PE_KR
    AsyCov=MyAsyCov;

PROC GLIMMIX DATA=WORK.BC ASYCOV;
	CLASS WP SP;
    MODEL Y = 
        X1 X2 X3 X4 X5 X6
        X1*X2 X1*X3 X1*X4 X1*X5 X1*X6 
        X2*X3 X2*X4 X2*X5 X2*X6
        X3*X4 X3*X5 X3*X6
        X4*X5 X4*X6 
        X5*X6 / DDFM=Kenwardroger SOLUTION;
    RANDOM INTERCEPT / SUBJECT=WP;
    RANDOM INTERCEPT / SUBJECT=SP;
RUN;

ODS OUTPUT
    ParameterEstimates=PE_SW;

PROC GLIMMIX DATA=WORK.BC;
	CLASS WP SP;
    MODEL Y = 
        X1 X2 X3 X4 X5 X6
        X1*X2 X1*X3 X1*X4 X1*X5 X1*X6 
        X2*X3 X2*X4 X2*X5 X2*X6
        X3*X4 X3*X5 X3*X6
        X4*X5 X4*X6 
        X5*X6 / DDFM=Satterthwaite SOLUTION;
    RANDOM INTERCEPT / SUBJECT=WP;
    RANDOM INTERCEPT / SUBJECT=SP;
RUN;

DATA PE_KR;
    SET PE_KR;
    FORMAT _NUMERIC_ BEST32.;
RUN;

DATA MyAsyCov;
    SET MyAsyCov;
    FORMAT _NUMERIC_ BEST32.;
RUN;

DATA PE_SW;
    SET PE_SW;
    FORMAT _NUMERIC_ BEST32.;
RUN;

PROC EXPORT DATA=PE_KR
    OUTFILE='/home/u64165441/Results battery cell sas kr.csv'
    DBMS=CSV
    REPLACE;
RUN;

PROC EXPORT DATA=MyAsyCov
    OUTFILE='/home/u64165441/Results battery cell sas asycov.csv'
    DBMS=CSV
    REPLACE;
RUN;

PROC EXPORT DATA=PE_SW
    OUTFILE='/home/u64165441/Results battery cell sas sw.csv'
    DBMS=CSV
    REPLACE;
RUN;