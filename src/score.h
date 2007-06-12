/**
   This software is distributed under the terms of the GNU GENERAL
   PUBLIC LICENSE Version 2, June 1991.  The terms of this license
   are in a file called COPYING which you should have received with
   this software.
*/
#include <R.h>
#include <Rdefines.h>
#include <stdio.h>
#include <ctype.h>   	/* character handling*/
#include <stdlib.h>     /* def of RAND_MAX */
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <limits.h>


extern SEXP showArgs(SEXP args);

extern void revstring(char *ch);
/*extern void merror(char *);*/		/** error handling */
extern int char2AA(char ch);
extern int AA2char(int x);
extern int maximum(int n,...);
extern double maximumdouble(int n,...);

extern SEXP alignSEXP(
	    SEXP scc1 /*AAstring 1*/
	    ,SEXP scc2 /*AAstring 2*/
	    ,SEXP ssimm /*simm*/
	    ,SEXP smdelta /*gap delta*/
	    ,SEXP sgapext /*gap extension costs*/
	    ,SEXP stype /*type of alignment*/
	    );

extern void globalB
(
 const char *cc1,
 const char *cc2,
 double *simm, /*similarity matrix*/
 int mdim,  /*dimension of similarity matrix*/
 double mdelta, /*gap penalty*/
 double gapext, /*gap extendsion costs*/
 const char *type, /*type of alignment*/
 char *al1, /*aligment for string 2*/
 char *al2, /*alignment fo r string 2*/
 int *all, /*alignment length*/
 double *score, /*score of alignment*/
 double *selfscore1, /*score of selfalignment*/
 double *selfscore2, /*score of selfalignment*/
 int *identity, /*nr of identities*/
 int *alignsimilarity, /*the similarity of the alignment*/
 char *errormsg
 );

extern void pozitiveScore
(
 const char *cc1
 ,const char *cc2
 ,double *simm /*similarity matrix*/
 ,int mdim  /*dimension of similarity matrix*/
 ,double mdelta /*gap penalty*/
 ,double gapext /*gap extendsion costs*/
 ,const char *align_type /*type of alignment*/
 ,double *score
 ,char *errormsg
 ,const char *score_t
 );
/*test*/

extern void globalScore
(
 const char *cc1
 ,const char *cc2
 ,double *simm /*similarity matrix*/
 ,int mdim  /*dimension of similarity matrix*/
 ,double mdelta /*gap penalty*/
 ,double gapext /*gap extendsion costs*/
 ,const char *alig_type /*type of alignment*/
 ,double *score /*score of alignment*/
 ,char *errormsg
 ,const char *score_t
 );

extern void identSimilarScore
(
 const char *cc1
 ,const char *cc2
 ,double *simm /*similarity matrix*/
 ,int mdim  /*dimension of similarity matrix*/
 ,double mdelta /*gap penalty*/
 ,double gapext /*gap extendsion costs*/
 ,const char *align_type /*type of alignment*/
 ,double *score
 ,char *errormsg
 ,const char *score_t
 );


/*Functions in basic h.*/
/**
* translates the sequence string to a index array.
*/
extern int seq2index(int *nc, const char *c, int lc);
/**
copies the weight from simm to dimm
*/
extern void initsim(double **sim, int dim, double *simm);
extern void printsim(double **sim, int test);
extern void printmat(double **sim, int n, int m);
extern void allocmatrixint(int ***mat, int x, int y);
extern void allocmatrixdouble(double ***mat, int x, int y);
extern void freematrixint(int ***mat,int x, int y);
extern void freematrixdouble(double ***mat,int x, int y);
/*extern int getselfalign(int *nc , int lc);*/

extern SEXP alignScoreSEXP(
	    SEXP scc1 /*AAstring 1*/
	    ,SEXP scc2 /*AAstring 2*/
	    ,SEXP ssimm /*simm*/
	    ,SEXP smdelta /*gap delta*/
	    ,SEXP sgapext /*gap extension costs*/
	    ,SEXP stype /*type of alignment*/
	    ,SEXP scoretype /*scoretype*/
	    );
