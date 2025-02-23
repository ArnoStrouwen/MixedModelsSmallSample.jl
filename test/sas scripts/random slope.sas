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

PROC CONTENTS DATA=WORK.SS; RUN;
%web_open_table(WORK.SS);

ODS OUTPUT
    ParameterEstimates=PE;
    
PROC GLIMMIX DATA=WORK.SS ASYCOV;
	CLASS subj;
    MODEL reaction = days / DDFM=Kenwardroger SOLUTION;
    RANDOM INTERCEPT / SUBJECT=subj;
    RANDOM days / SUBJECT=subj;
RUN;

DATA PE;
    SET PE;
    FORMAT Estimate StdErr DF tValue Probt BEST12.10;
RUN;

PROC EXPORT DATA=PE
    OUTFILE='/home/u64165441/Results sleep study sas kr.csv'
    DBMS=CSV
    REPLACE;
RUN;

ODS OUTPUT
    ParameterEstimates=PE;
    
PROC GLIMMIX DATA=WORK.SS ASYCOV;
	CLASS subj;
    MODEL reaction = days / DDFM=Satterthwaite SOLUTION;
    RANDOM INTERCEPT / SUBJECT=subj;
    RANDOM days / SUBJECT=subj;
RUN;

DATA PE;
    SET PE;
    FORMAT Estimate StdErr DF tValue Probt BEST12.10;
RUN;

PROC EXPORT DATA=PE
    OUTFILE='/home/u64165441/Results sleep study sas sw.csv'
    DBMS=CSV
    REPLACE;
RUN;