#include "score.h"

SEXP alignScoreSEXP(
	    SEXP scc1 /*AAstring 1*/
	    ,SEXP scc2 /*AAstring 2*/
	    ,SEXP ssimm /*simm*/
	    ,SEXP smdelta /*gap delta*/
	    ,SEXP sgapext /*gap extension costs*/
	    ,SEXP salign_type /*type of alignment*/
	    ,SEXP scoretype /*scoretype*/
	    )
{
  double *simm = REAL(ssimm);
  const char *align_type = CHAR(STRING_ELT(salign_type,0)); /*alignment type*/
  const char *score_t = CHAR(STRING_ELT(scoretype,0)); /*scoretype*/
  int mdim = (int) sqrt((double)length(ssimm));
  double delta = REAL(smdelta)[0];
  double gapext = REAL(sgapext)[0];
  /** getting the char */
  const char *cc1,*cc2;  
  SEXP score;
  char *errormsg;
  double score_d=0.0;
  errormsg="                               ";
  if(!isString(scc1))
    error("is not a string\n");
  PROTECT(ssimm = AS_NUMERIC(ssimm));
  PROTECT(smdelta = AS_NUMERIC(smdelta));
  PROTECT(sgapext = AS_NUMERIC(sgapext));
  PROTECT(salign_type= AS_CHARACTER(salign_type));
  PROTECT(scoretype = AS_CHARACTER(scoretype));
  cc1 = CHAR(STRING_ELT(scc1,0));
  cc2 = CHAR(STRING_ELT(scc2,0));
  PROTECT(score = NEW_NUMERIC(1));
  if(strcmp(score_t,"score")==0 || strcmp(score_t,"scoreN")==0)
    {
      globalScore(
		  cc1
		  ,cc2
		  ,simm /*similarity matrix*/
		  ,mdim  /*dimension of similarity matrix*/
		  ,delta /*gap penalty*/
		  ,gapext /*gap extendsion costs*/
		  ,align_type /*alignment type*/
		  ,&score_d /*score result*/
		  ,errormsg
		  ,score_t
		  );
    }
  else if(strcmp(score_t,"similarity")==0 || strcmp(score_t,"identity")==0)
    {
      identSimilarScore(
			cc1
			,cc2
			,simm /*similarity matrix*/
			,mdim  /*dimension of similarity matrix*/
			,delta /*gap penalty*/
			,gapext /*gap extendsion costs*/
			,align_type /*type of alignment*/
			,&score_d
			,errormsg
			,score_t
			);
    }
  else if(strcmp(score_t,"pozitive")==0)
    {
      pozitiveScore(
		    cc1
		    ,cc2
		    ,simm /*similarity matrix*/
		    ,mdim  /*dimension of similarity matrix*/
		    ,delta /*gap penalty*/
		    ,gapext /*gap extendsion costs*/
		    ,align_type /*type of alignment*/
		    ,&score_d
		    ,errormsg
		    ,score_t
		    );
    }
  else
    {
      printf("No such score\n");
    }

  REAL(score)[0]=score_d; /*score of alignment*/
  UNPROTECT(6);
  /*setting the answer*/
  return(score);
  /*  return(testMemalloc());*/
}



SEXP alignSEXP(
	    SEXP scc1 /*AAstring 1*/
	    ,SEXP scc2 /*AAstring 2*/
	    ,SEXP ssimm /*simm*/
	    ,SEXP smdelta /*gap delta*/
	    ,SEXP sgapext /*gap extension costs*/
	    ,SEXP stype /*type of alignment*/
	    )
{
  double *simm;

  int mdim= (int) sqrt((double)length(ssimm));
  double delta = REAL(smdelta)[0];
  double gapext = REAL(sgapext)[0];
  const char *type = CHAR(STRING_ELT(stype,0));
  /** getting the char */
  const char *cc1,*cc2;  
  char *al1;
  char *al2;
  SEXP ral1;
  SEXP ral2;
  SEXP all;
  SEXP score;
  SEXP selfscore1;
  SEXP selfscore2;
  SEXP identity;
  SEXP alignsimilarity;
  SEXP err;
  char *errormsg;
  int ident_i=0,simil_i=0,all_i=0;
  double score_d=0.0,self1_d = 0.0,self2_d=0.0;
  SEXP ans,ansnames;

  errormsg="                               ";
  simm = REAL(ssimm);
  if(!isString(scc1))
    error("is not a string\n");

  PROTECT(ssimm = AS_NUMERIC(ssimm));
  PROTECT(smdelta = AS_NUMERIC(smdelta));
  PROTECT(sgapext = AS_NUMERIC(sgapext));
  PROTECT(stype= AS_CHARACTER(stype));

  cc1 = CHAR(STRING_ELT(scc1,0));
  cc2 = CHAR(STRING_ELT(scc2,0));
  al1 = R_alloc((strlen(cc1) + strlen(cc2)),sizeof(char));
  al2 = R_alloc((strlen(cc1) + strlen(cc2)),sizeof(char));
    
  PROTECT(ral1 = NEW_CHARACTER(1));
  PROTECT(ral2 = NEW_CHARACTER(1));
  PROTECT(all = NEW_INTEGER(1));
  PROTECT(score = NEW_NUMERIC(1));
  PROTECT(selfscore1 = NEW_NUMERIC(1));
  PROTECT(selfscore2 = NEW_NUMERIC(1));
  PROTECT(identity = NEW_INTEGER(1));
  PROTECT(alignsimilarity = NEW_INTEGER(1));
  PROTECT(err = NEW_CHARACTER(1));
  
  globalB(
        cc1
        ,cc2
        ,simm /*similarity matrix*/
        ,mdim  /*dimension of similarity matrix*/
        ,delta /*gap penalty*/
        ,gapext /*gap extendsion costs*/
        ,type /*type of alignment 1,2,3,4*/
        ,al1 /*aligment for string 2*/
        ,al2 /*alignment for string 2*/
	,&all_i
	,&score_d
	,&self1_d
	,&self2_d
	,&ident_i
	,&simil_i
	,errormsg
    );
  /*
    printf("al1:%s\n",al1 );
    printf("al2:%s\n",al2 );
    printf("length:%d\n",all_i );
    printf("score:%f\n",score_d );
    printf("self1:%f\n",self1_d );
    printf("self2:%f\n",self2_d );
    printf("ident:%i\n",ident_i);
    printf("simil:%i\n",simil_i);
  */

  INTEGER(all)[0]=all_i; /*alignment length*/
  REAL(score)[0]=score_d; /*score of alignment*/
  REAL(selfscore1)[0] = self1_d; /*score of selfalignment*/
  REAL(selfscore2)[0] = self2_d; /*score of selfalignment*/
  INTEGER(identity)[0] = ident_i;
  INTEGER(alignsimilarity)[0] = simil_i; /*the similarity of the alignment*/
  


  /*setting the answer*/
  
  PROTECT(ans = allocVector(VECSXP, 9));
  PROTECT(ansnames = allocVector(STRSXP, 9));
  
  SET_STRING_ELT(ral1,0,COPY_TO_USER_STRING(al1));
  SET_VECTOR_ELT(ans,0,ral1);
  SET_STRING_ELT(ansnames, 0, COPY_TO_USER_STRING("alig1"));
  
  SET_STRING_ELT(ral2,0,COPY_TO_USER_STRING(al2));
  SET_VECTOR_ELT(ans,1,ral2);
  SET_STRING_ELT(ansnames, 1, COPY_TO_USER_STRING("alig2"));
  
  SET_VECTOR_ELT(ans,2,all);
  SET_STRING_ELT(ansnames, 2, COPY_TO_USER_STRING("all"));
  
  SET_VECTOR_ELT(ans,3,score);
  SET_STRING_ELT(ansnames, 3, COPY_TO_USER_STRING("score"));
  
  SET_VECTOR_ELT(ans,4,selfscore1);
  SET_STRING_ELT(ansnames, 4, COPY_TO_USER_STRING("selfscore1"));
  
  SET_VECTOR_ELT(ans,5,selfscore2);
  SET_STRING_ELT(ansnames, 5, COPY_TO_USER_STRING("selfscore2"));
  
  SET_VECTOR_ELT(ans,6,identity);
  SET_STRING_ELT(ansnames, 6, COPY_TO_USER_STRING("identity"));
  
  SET_VECTOR_ELT(ans,7,alignsimilarity);
  SET_STRING_ELT(ansnames, 7, COPY_TO_USER_STRING("alignsimilarity"));
  
  SET_STRING_ELT(err,0,COPY_TO_USER_STRING(errormsg));
  SET_VECTOR_ELT(ans,8,err);
  SET_STRING_ELT(ansnames, 8, COPY_TO_USER_STRING("errmsg"));
  
  setAttrib(ans, R_NamesSymbol, ansnames);    
  UNPROTECT(15);
  
  return(ans);
  /*  return(testMemalloc());*/
}
