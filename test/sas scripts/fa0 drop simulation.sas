/*
FA0(2) boundary/drop simulation

Goal
- Fit the same random intercept + random slope model across many simulated datasets
- Increase correlation between intercept and slope across datasets
- Export AsyCov and CovParms for PROC MIXED and PROC GLIMMIX separately

How to run
- Set &INFILE and &OUTDIR to your paths
- Run in SAS; outputs are written as CSV to &OUTDIR
*/

/* Match other repo SAS scripts: hard-coded input/output paths */

DATA WORK.SIMDATA;
    INFILE '/home/u64165441/Data fa0 drop simulation.csv' DLM=',' DSD FIRSTOBS=2 TRUNCOVER;
    INPUT sim_id rho rep subject time y;
RUN;

PROC SORT DATA=WORK.SIMDATA;
    BY sim_id subject time;
RUN;

/* ============================================================ */
/* PROC MIXED                                                    */
/* ============================================================ */

ODS OUTPUT
    CovParms=Mixed_CovParms
    AsyCov=Mixed_AsyCov;

PROC MIXED DATA=WORK.SIMDATA ASYCOV;
    BY sim_id;
    CLASS subject;
    MODEL y = time / DDFM=KENWARDROGER SOLUTION;
    RANDOM INTERCEPT time / TYPE=FA0(2) SUBJECT=subject;
RUN;

DATA Mixed_CovParms; SET Mixed_CovParms; FORMAT _NUMERIC_ BEST32.; RUN;
DATA Mixed_AsyCov; SET Mixed_AsyCov; FORMAT _NUMERIC_ BEST32.; RUN;

PROC EXPORT DATA=Mixed_CovParms
    OUTFILE='/home/u64165441/Results fa0 drop simulation mixed covparms.csv'
    DBMS=CSV
    REPLACE;
RUN;

PROC EXPORT DATA=Mixed_AsyCov
    OUTFILE='/home/u64165441/Results fa0 drop simulation mixed asycov.csv'
    DBMS=CSV
    REPLACE;
RUN;

/* ============================================================ */
/* PROC GLIMMIX                                                  */
/* ============================================================ */

ODS OUTPUT
    CovParms=Glimmix_CovParms
    AsyCov=Glimmix_AsyCov;

PROC GLIMMIX DATA=WORK.SIMDATA ASYCOV;
    BY sim_id;
    CLASS subject;
    MODEL y = time / DDFM=KENWARDROGER SOLUTION;
    RANDOM INTERCEPT time / TYPE=FA0(2) SUBJECT=subject;
RUN;

DATA Glimmix_CovParms; SET Glimmix_CovParms; FORMAT _NUMERIC_ BEST32.; RUN;
DATA Glimmix_AsyCov; SET Glimmix_AsyCov; FORMAT _NUMERIC_ BEST32.; RUN;

PROC EXPORT DATA=Glimmix_CovParms
    OUTFILE='/home/u64165441/Results fa0 drop simulation glimmix covparms.csv'
    DBMS=CSV
    REPLACE;
RUN;

PROC EXPORT DATA=Glimmix_AsyCov
    OUTFILE='/home/u64165441/Results fa0 drop simulation glimmix asycov.csv'
    DBMS=CSV
    REPLACE;
RUN;
