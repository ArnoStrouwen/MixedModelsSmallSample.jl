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

PROC CONTENTS DATA=WORK.BC; RUN;
%web_open_table(WORK.BC);

ODS OUTPUT
    ParameterEstimates=PE;
    
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

DATA PE;
    SET PE;
    FORMAT Estimate StdErr DF tValue Probt BEST12.10;
RUN;

PROC EXPORT DATA=PE
    OUTFILE='/home/u64165441/Results battery cell sas.csv'
    DBMS=CSV
    REPLACE;
RUN;