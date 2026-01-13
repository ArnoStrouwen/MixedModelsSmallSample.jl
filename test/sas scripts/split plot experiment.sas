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

ODS OUTPUT
    ParameterEstimates=PE_KR
    AsyCov=MyAsyCov;
    
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

ODS OUTPUT
    ParameterEstimates=PE_SW;

PROC GLIMMIX DATA=WORK.WT;
	CLASS WP;
    MODEL EFFICIENCY = 
        FRH RRH YA GC
        FRH*RRH FRH*YA FRH*GC
        RRH*YA RRH*GC
        YA*GC
        FRH*FRH RRH*RRH YA*YA GC*GC / DDFM=Satterthwaite SOLUTION;

    RANDOM INTERCEPT / SUBJECT=WP;
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
    OUTFILE='/home/u64165441/Results wind tunnel sas kr.csv'
    DBMS=CSV
    REPLACE;
RUN;

PROC EXPORT DATA=MyAsyCov
    OUTFILE='/home/u64165441/Results wind tunnel sas asycov.csv'
    DBMS=CSV
    REPLACE;
RUN;

PROC EXPORT DATA=PE_SW
    OUTFILE='/home/u64165441/Results wind tunnel sas sw.csv'
    DBMS=CSV
    REPLACE;
RUN;