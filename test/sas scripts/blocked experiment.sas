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

ODS OUTPUT 
    AsyCov=MyAsyCov
    ParameterEstimates=PE_KR;

PROC GLIMMIX DATA=WORK.PD ASYCOV;
	CLASS Day;
    MODEL LEI = 
        FR MC SS
        FR*MC FR*SS
        MC*SS
        FR*FR MC*MC SS*SS / DDFM=Kenwardroger SOLUTION;
    RANDOM INTERCEPT / SUBJECT=Day;
RUN;

ODS OUTPUT 
    ParameterEstimates=PE_SW;

PROC GLIMMIX DATA=WORK.PD;
	CLASS Day;
    MODEL LEI = 
        FR MC SS
        FR*MC FR*SS
        MC*SS
        FR*FR MC*MC SS*SS / DDFM=Satterthwaite SOLUTION;
    RANDOM INTERCEPT / SUBJECT=Day;
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

PROC EXPORT DATA=MyAsyCov
    OUTFILE="/home/u64165441/Results pastry dough sas asycov.csv"
    DBMS=CSV
    REPLACE;
RUN;

PROC EXPORT DATA=PE_KR
    OUTFILE="/home/u64165441/Results pastry dough sas kr.csv"
    DBMS=CSV
    REPLACE;
RUN;

PROC EXPORT DATA=PE_SW
    OUTFILE="/home/u64165441/Results pastry dough sas sw.csv"
    DBMS=CSV
    REPLACE;
RUN;