%web_drop_table(WORK.WT);
FILENAME REFFILE '/home/u64165441/Data wind tunnel Chapter 10.csv';

PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV
	OUT=WORK.WT;
	GETNAMES=YES;
RUN;

DATA WORK.WT;
    SET WORK.WT;
    RENAME "Whole Plots"n = WP;
RUN;

PROC CONTENTS DATA=WORK.WT; RUN;
%web_open_table(WORK.WT);

ODS OUTPUT
    ParameterEstimates=PE;
    
PROC GLIMMIX DATA=WORK.WT ASYCOV;
	CLASS WP;
    MODEL EFFICIENCY = 
        FRH RRH YA GC
        FRH*RRH FRH*YA FRH*GC
        RRH*YA RRH*GC
        YA*GC
        FRH*FRH RRH*RRH YA*YA GC*GC / DDFM=Kenwardroger SOLUTION;

    RANDOM INTERCEPT / SUBJECT=WP;
RUN;

DATA PE;
    SET PE;
    FORMAT Estimate StdErr DF tValue Probt BEST12.10;
RUN;

PROC EXPORT DATA=PE
    OUTFILE='/home/u64165441/Results wind tunnel sas kr.csv'
    DBMS=CSV
    REPLACE;
RUN;

ODS OUTPUT
    ParameterEstimates=PE;
    
PROC GLIMMIX DATA=WORK.WT ASYCOV;
	CLASS WP;
    MODEL EFFICIENCY = 
        FRH RRH YA GC
        FRH*RRH FRH*YA FRH*GC
        RRH*YA RRH*GC
        YA*GC
        FRH*FRH RRH*RRH YA*YA GC*GC / DDFM=Satterthwaite SOLUTION;

    RANDOM INTERCEPT / SUBJECT=WP;
RUN;

DATA PE;
    SET PE;
    FORMAT Estimate StdErr DF tValue Probt BEST12.10;
RUN;

PROC EXPORT DATA=PE
    OUTFILE='/home/u64165441/Results wind tunnel sas sw.csv'
    DBMS=CSV
    REPLACE;
RUN;