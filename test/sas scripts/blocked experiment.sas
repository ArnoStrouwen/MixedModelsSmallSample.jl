%web_drop_table(WORK.PD);
FILENAME REFFILE '/home/u64165441/Data Pastry Dough Experiment Chapter 7.csv';

PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV
	OUT=WORK.PD;
	GETNAMES=YES;
RUN;

DATA WORK.PD;
    SET WORK.PD;
    RENAME "Flow Rate"n = FR;
    RENAME "Moisture Content"n = MC;
    RENAME "Screw Speed"n = SS;
    RENAME "Longitudinal Expansion Index"n = LEI;
RUN;

PROC STDIZE DATA=WORK.PD OUT=WORK.PD METHOD=RANGE MULT=2 ADD=-1;
	var FR MC SS;
RUN;
PROC CONTENTS DATA=WORK.PD; RUN;
%web_open_table(WORK.PD);

ODS OUTPUT
    ParameterEstimates=PE;
    
PROC GLIMMIX DATA=WORK.PD ASYCOV;
	CLASS Day;
    MODEL LEI = 
        FR MC SS
        FR*MC FR*SS
        MC*SS
        FR*FR MC*MC SS*SS / DDFM=Kenwardroger SOLUTION;
    RANDOM INTERCEPT / SUBJECT=Day;

RUN;

DATA PE;
    SET PE;
    FORMAT Estimate StdErr DF tValue Probt BEST12.10;
RUN;

PROC EXPORT DATA=PE
    OUTFILE='/home/u64165441/Results pastry dough sas kr.csv'
    DBMS=CSV
    REPLACE;
RUN;

ODS OUTPUT
    ParameterEstimates=PE;
    
PROC GLIMMIX DATA=WORK.PD ASYCOV;
	CLASS Day;
    MODEL LEI = 
        FR MC SS
        FR*MC FR*SS
        MC*SS
        FR*FR MC*MC SS*SS / DDFM=Satterthwaite SOLUTION;
    RANDOM INTERCEPT / SUBJECT=Day;

RUN;

DATA PE;
    SET PE;
    FORMAT Estimate StdErr DF tValue Probt BEST12.10;
RUN;

PROC EXPORT DATA=PE
    OUTFILE='/home/u64165441/Results pastry dough sas sw.csv'
    DBMS=CSV
    REPLACE;
RUN;