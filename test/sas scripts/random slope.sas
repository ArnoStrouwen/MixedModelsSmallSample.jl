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

ODS OUTPUT
    ParameterEstimates=PE;
    
PROC GLIMMIX DATA=WORK.SS ASYCOV;
	CLASS subj;
    MODEL reaction = days / DDFM=Kenwardroger SOLUTION;
    RANDOM INTERCEPT days / type=un SUBJECT=subj;
RUN;

DATA PE;
    SET PE;
    FORMAT Estimate StdErr DF tValue Probt BEST12.10;
RUN;

PROC EXPORT DATA=PE
    OUTFILE='/home/u64165441/Results sleep study corr sas kr.csv'
    DBMS=CSV
    REPLACE;
RUN;

ODS OUTPUT
    ParameterEstimates=PE;
    
PROC GLIMMIX DATA=WORK.SS ASYCOV;
	CLASS subj;
    MODEL reaction = days / DDFM=Satterthwaite SOLUTION;
    RANDOM INTERCEPT days / type=un SUBJECT=subj;
RUN;

DATA PE;
    SET PE;
    FORMAT Estimate StdErr DF tValue Probt BEST12.10;
RUN;

PROC EXPORT DATA=PE
    OUTFILE='/home/u64165441/Results sleep study corr sas sw.csv'
    DBMS=CSV
    REPLACE;
RUN;

ODS OUTPUT
    ParameterEstimates=PE;
    
PROC GLIMMIX DATA=WORK.SS ASYCOV;
	CLASS subj;
    MODEL reaction = days days*days / DDFM=Kenwardroger SOLUTION;
    RANDOM INTERCEPT / SUBJECT=subj;
    RANDOM days / SUBJECT=subj;
    RANDOM days*days / SUBJECT=subj;
RUN;

DATA PE;
    SET PE;
    FORMAT Estimate StdErr DF tValue Probt BEST12.10;
RUN;

PROC EXPORT DATA=PE
    OUTFILE='/home/u64165441/Results sleep study quadratic sas kr.csv'
    DBMS=CSV
    REPLACE;
RUN;

ODS OUTPUT
    ParameterEstimates=PE;
    
PROC GLIMMIX DATA=WORK.SS ASYCOV;
	CLASS subj;
    MODEL reaction = days days*days/ DDFM=Satterthwaite SOLUTION;
    RANDOM INTERCEPT / SUBJECT=subj;
    RANDOM days / SUBJECT=subj;
    RANDOM days*days / SUBJECT=subj;
RUN;

DATA PE;
    SET PE;
    FORMAT Estimate StdErr DF tValue Probt BEST12.10;
RUN;

PROC EXPORT DATA=PE
    OUTFILE='/home/u64165441/Results sleep study quadratic sas sw.csv'
    DBMS=CSV
    REPLACE;
RUN;

ODS OUTPUT
    ParameterEstimates=PE;
    
PROC GLIMMIX DATA=WORK.SS ASYCOV;
	CLASS subj;
    MODEL reaction = days days*days / DDFM=Kenwardroger SOLUTION;
    RANDOM INTERCEPT days days*days / type=un SUBJECT=subj;
RUN;

DATA PE;
    SET PE;
    FORMAT Estimate StdErr DF tValue Probt BEST12.10;
RUN;

PROC EXPORT DATA=PE
    OUTFILE='/home/u64165441/Results sleep study corr quadratic sas kr.csv'
    DBMS=CSV
    REPLACE;
RUN;

ODS OUTPUT
    ParameterEstimates=PE;
    
PROC GLIMMIX DATA=WORK.SS ASYCOV;
	CLASS subj;
    MODEL reaction = days days*days / DDFM=Satterthwaite SOLUTION;
    RANDOM INTERCEPT days days*days / type=un SUBJECT=subj;
RUN;

DATA PE;
    SET PE;
    FORMAT Estimate StdErr DF tValue Probt BEST12.10;
RUN;

PROC EXPORT DATA=PE
    OUTFILE='/home/u64165441/Results sleep study corr quadratic sas sw.csv'
    DBMS=CSV
    REPLACE;
RUN;

ODS OUTPUT
    ParameterEstimates=PE;
    
PROC GLIMMIX DATA=WORK.SS ASYCOV;
	CLASS subj;
    MODEL reaction = days days*days / DDFM=Kenwardroger SOLUTION;
    RANDOM INTERCEPT / type=un SUBJECT=subj;
    RANDOM days days*days / type=un SUBJECT=subj;
RUN;

DATA PE;
    SET PE;
    FORMAT Estimate StdErr DF tValue Probt BEST12.10;
RUN;

PROC EXPORT DATA=PE
    OUTFILE='/home/u64165441/Results sleep study some corr quadratic sas kr.csv'
    DBMS=CSV
    REPLACE;
RUN;

ODS OUTPUT
    ParameterEstimates=PE;
    
PROC GLIMMIX DATA=WORK.SS ASYCOV;
	CLASS subj;
    MODEL reaction = days days*days / DDFM=Satterthwaite SOLUTION;
    RANDOM INTERCEPT / type=un SUBJECT=subj;
    RANDOM days days*days / type=un SUBJECT=subj;
RUN;

DATA PE;
    SET PE;
    FORMAT Estimate StdErr DF tValue Probt BEST12.10;
RUN;

PROC EXPORT DATA=PE
    OUTFILE='/home/u64165441/Results sleep study some corr quadratic sas sw.csv'
    DBMS=CSV
    REPLACE;
RUN;
