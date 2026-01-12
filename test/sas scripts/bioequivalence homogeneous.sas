%web_drop_table(WORK.WT);
FILENAME REFFILE '/home/u64165441/Data bioequivalence.csv';

PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV
	OUT=WORK.WT;
	GETNAMES=YES;
RUN;

DATA WORK.WT;
    SET WORK.WT;
RUN;

PROC CONTENTS DATA=WORK.WT; RUN;
%web_open_table(WORK.WT);

ODS OUTPUT
    Tests3=Tests3;
PROC GLIMMIX DATA=WORK.WT ASYCOV;
	CLASS subject sequence period formulation ;
    MODEL log_data = period formulation sequence / DDFM=kenwardroger SOLUTION;
    RANDOM INTERCEPT / SUBJECT=subject;
RUN;

DATA Tests3;
    SET Tests3;
    FORMAT DenDF FValue ProbF BEST12.10;
RUN;

proc contents data=Tests3;
run;

PROC EXPORT DATA=Tests3
    OUTFILE='/home/u64165441/Results bioequivalence homogeneous sas kr.csv'
    DBMS=CSV
    REPLACE;
RUN;

ODS OUTPUT
    Tests3=Tests3;
PROC GLIMMIX DATA=WORK.WT ASYCOV;
	CLASS subject sequence period formulation ;
    MODEL log_data = period formulation sequence / DDFM=satterthwaite SOLUTION;
    RANDOM INTERCEPT / SUBJECT=subject;
RUN;

DATA Tests3;
    SET Tests3;
    FORMAT DenDF FValue ProbF BEST12.10;
RUN;

proc contents data=Tests3;
run;

PROC EXPORT DATA=Tests3
    OUTFILE='/home/u64165441/Results bioequivalence homogeneous sas sw.csv'
    DBMS=CSV
    REPLACE;
RUN;

PROC MIXED DATA=WORK.WT ASYCOV;
	CLASS subject sequence(ref="ABAB") period(ref="1") formulation(ref="R");
	MODEL log_data = sequence period formulation/ DDFM=SATTERTH SOLUTION;
	RANDOM formulation/TYPE=FA0(2) SUB=subject G;
	REPEATED/R=1,2,3,4,5,6,7,8, GRP=formulation;
RUN;

PROC GLIMMIX DATA=WORK.WT ASYCOV outdesign(z)=zmatrix;
	CLASS subject sequence(ref="ABAB") period(ref="1") formulation(ref="R");
	MODEL log_data = sequence period formulation/ DDFM=SATTERTH SOLUTION;
	RANDOM formulation/TYPE=FA0(2) SUB=subject G;
RUN;

proc print data=zmatrix;
run;

