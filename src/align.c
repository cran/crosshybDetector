#include "score.h"

const double MIN_DOUBLE = -1000000;

typedef struct
{
  double **h; /*scoring matrix*/
  double **Ix; /*matrix necessary for affince gap penalty*/
  double **Iy; /*matrix necessary for affine gap penalty*/
  int **itrace; /*tracback matrix*/
  int **jtrace; /*tracback matrix*/
  int lc1; /*size 1 matrix*/
  int lc2; /*size 2 matrix*/
  double **sim; /*similarity matrix*/
  int dim; /*dim of matrix*/
  const char *c1; /*string1*/
  const char *c2; /*string2*/
  char *ac1; /*alignment string*/
  char *ac2; /*alignment string*/
  int *nc1; /*index representation of c1*/
  int *nc2; /*index representation of c2*/
  double delta; /*linear gap penalty*/
  double gapext; /*gap extension costs*/
  double score; /*the resulting score*/
  int si; /*starting i  postion for backtrace*/
  int sj; /*starting j postion for backtrace*/
  /*for POZITIVE SCORE*/
  int *ung1;
  int *ung2;
  int ungl;
  short pozitiv;
} DYNAM;

void printstrut(DYNAM *test)
{
  printf("----------------------------------\n");
   printf("in printstrut:\n");
   printf("STRING 1 : %s\n" , (*test).c1 );
   printf("STRING 2 : %s\n" , (*test).c2 );
   printf("SIZE 1 : %d\n" , test->lc1 );
   printf("SIZE 2 : %d\n" , test->lc2 );
   printf("Gap Penalty: %f\n" , test->delta );
   printf("Extension : %f\n" , test->gapext );
   printf("----------------------------------\n");
}


void allocDYNAM(DYNAM *p_dynam)
{
  /*DYNAM tmp = *p_dynam;*/
    /*printf("dims %d %d\n", p_dynam->lc1 , p_dynam->lc2);*/
    allocmatrixdouble(&(p_dynam->h) , p_dynam->lc1 , p_dynam->lc2);
    allocmatrixdouble(&(p_dynam->Ix), p_dynam->lc1 , p_dynam->lc2);
    allocmatrixdouble(&(p_dynam->Iy), p_dynam->lc1 , p_dynam->lc2);
    allocmatrixdouble(&(p_dynam->sim), p_dynam->dim , p_dynam->dim );
    allocmatrixint(&(p_dynam->itrace), p_dynam->lc1, p_dynam->lc2);
    allocmatrixint(&(p_dynam->jtrace), p_dynam->lc1, p_dynam->lc2);
    p_dynam->nc1=malloc((p_dynam->lc1)*sizeof(int));
    p_dynam->nc2=malloc((p_dynam->lc2)*sizeof(int));
}


void freeDYNAM(DYNAM *p_dynam)
{
    freematrixdouble(&(p_dynam->h) , p_dynam->lc1 , p_dynam->lc2);
    freematrixdouble(&(p_dynam->Ix), p_dynam->lc1 , p_dynam->lc2);
    freematrixdouble(&(p_dynam->Iy), p_dynam->lc1 , p_dynam->lc2);
    freematrixdouble(&(p_dynam->sim), p_dynam->dim , p_dynam->dim );
    freematrixint(&(p_dynam->itrace), p_dynam->lc1, p_dynam->lc2);
    freematrixint(&(p_dynam->jtrace), p_dynam->lc1, p_dynam->lc2);
    free((void *) p_dynam->nc1);
    free((void *) p_dynam->nc2);
}

void initglobalA(DYNAM *p_dynam)
{
  int i;
  p_dynam->h[0][0] = 0.0;
  p_dynam->Ix[0][0] = MIN_DOUBLE;
  p_dynam->Iy[0][0] = MIN_DOUBLE;
  for(i=1;i < p_dynam->lc1; i++ )
    {
      p_dynam->Iy[i][0] = (p_dynam->delta) + (i-1)*(p_dynam->gapext);
      p_dynam->h[i][0] = MIN_DOUBLE;
      p_dynam->Ix[i][0] = MIN_DOUBLE;
      p_dynam->itrace[i][0] = i-1;
    }
  for(i=1; i< p_dynam->lc2;i++)
    {
      p_dynam->Ix[0][i] = (p_dynam->delta) + (i-1)*(p_dynam->gapext);
      p_dynam->h[0][i] = MIN_DOUBLE;
      p_dynam->Iy[0][i] = MIN_DOUBLE;
      p_dynam->jtrace[0][i] = i-1;
    }
}

double hglob(DYNAM *p_dynam, int i, int j)
{
  double sim,m,x,y;
  sim = p_dynam->sim[p_dynam->nc1[i]][p_dynam->nc2[j]];
  /*  printf("%c %c %f\n",AA2char(p_dynam->nc1[i]),AA2char(p_dynam->nc2[j]),sim);*/
  m = p_dynam->h[i-1][j-1] + sim;
  x = p_dynam->Ix[i-1][j-1] + sim;
  y = p_dynam->Iy[i-1][j-1] + sim;
  return(maximumdouble(3,m,x,y));
}

double ixglob(DYNAM *p_dynam, int i, int j)
{
  double m,x,y;
  m=p_dynam->h[i][j-1] +  p_dynam->delta;
  y=p_dynam->Iy[i][j-1] + p_dynam->delta;  
  x=p_dynam->Ix[i][j-1] + p_dynam->gapext;
  return(maximumdouble(3,m,x,y));
}

double iyglob(DYNAM *p_dynam, int i, int j)
{
  double m,y,x;
  m=p_dynam->h[i-1][j] +   p_dynam->delta;
  x=p_dynam->Ix[i-1][j] +  p_dynam->delta;
  y=p_dynam->Iy[i-1][j] + p_dynam->gapext;
  return(maximumdouble(3,m,y,x));
}


void settrace(DYNAM *p_dynam, int i, int j)
{
  
	  double max = maximumdouble(3,p_dynam->h[i][j],p_dynam->Ix[i][j],p_dynam->Iy[i][j]);
	  if(fabs(max - p_dynam->h[i][j])<1e-14)
	    {
	      p_dynam->itrace[i][j] = i-1;
	      p_dynam->jtrace[i][j] = j-1;
	      return;
	    }
	  else 
	    if(fabs(max - p_dynam->Ix[i][j])<1e-14)
	    {
	      p_dynam->itrace[i][j] = i;
	      p_dynam->jtrace[i][j] = j-1;
	      return;
	    }
	  else 
	    if(fabs(max - p_dynam->Iy[i][j])<1e-14)
	    {
	      p_dynam->itrace[i][j] = i-1;
	      p_dynam->jtrace[i][j] = j;
	      return;
	    }
	    else{printf("Something spoiled\n");}
}

void globalAlign(DYNAM *p_dynam)
{
  int i,j;
  int lc1,lc2;
  double max, h, Ix, Iy;

  initglobalA(p_dynam); 
  lc1=p_dynam->lc1;
  lc2=p_dynam->lc2;
  p_dynam->si=lc1-1;
  p_dynam->sj=lc2-1;
  for(i=1; i <lc1; i++)
    {
      for(j=1; j <lc2; j++)
	{
	  p_dynam->h[i][j]= hglob(p_dynam,i,j);
	  p_dynam->Ix[i][j]= ixglob(p_dynam,i,j);
	  p_dynam->Iy[i][j]= iyglob(p_dynam,i,j);
	  settrace(p_dynam,i,j);
	}
    }
  h = p_dynam->h[lc1-1][lc2-1];
  Ix = p_dynam->Ix[lc1-1][lc2-1];
  Iy = p_dynam->Iy[lc1-1][lc2-1];
  max=maximumdouble(3,h,Ix,Iy);
  p_dynam->score=max;
}

/*backtracer if the POZITIVE score are computed*/
void tracerPOZITIV(DYNAM *p_dynam,int *i, int *j, int *a, int *b)
{
  if((p_dynam->itrace[*i][*j] == (*i)-1) && (p_dynam->jtrace[*i][*j]== (*j)-1) )
    {
      p_dynam->ung1[p_dynam->ungl] = p_dynam->nc1[(*i)];
      p_dynam->ung2[p_dynam->ungl] = p_dynam->nc2[(*j)];
      p_dynam->ungl++;
      (*i)--;
      (*j)--;
      return;
    }
  else if((p_dynam->itrace[*i][*j]== *i) && (p_dynam->jtrace[*i][*j]== (*j)-1))
    {
      (*j)--;
      return;
    }
  else if((p_dynam->itrace[*i][*j] == *i-1 )  && ( p_dynam->jtrace[*i][*j] == *j))
    {
      (*i)--;
      return;
    }
}

/*backtracer if the alignment are computed.*/
void tracer(DYNAM *p_dynam,int *i, int *j, int *a, int *b)
{
  if((p_dynam->itrace[*i][*j] == (*i)-1) && (p_dynam->jtrace[*i][*j]== (*j)-1) )
    {
      (*i)--;
      (*j)--;
      p_dynam->ac1[(*a)++] = p_dynam->c1[(*i)];
      p_dynam->ac2[(*b)++] = p_dynam->c2[(*j)];
      return;
    }
  if((p_dynam->itrace[*i][*j]== *i) && (p_dynam->jtrace[*i][*j]== (*j)-1))
    {
      (*j)--;
      p_dynam->ac1[(*a)++] = '-';
      p_dynam->ac2[(*b)++] = p_dynam->c2[*j];
      return;
    }
  if((p_dynam->itrace[*i][*j] == *i-1 )  && ( p_dynam->jtrace[*i][*j] == *j))
    {
      (*i)--;
      p_dynam->ac1[(*a)++] = p_dynam->c1[*i];
      p_dynam->ac2[(*b)++] = '-';
      return;
    }
}

void tracebackPOZITIVOverlap(DYNAM *p_dynam)
{
  /**
     finding traceback for global alignment
   */
  int i,j,a=0,b=0;
  i=p_dynam->si;
  j=p_dynam->sj;
  /*dp the backtrack*/
  do{
    tracerPOZITIV(p_dynam,&i,&j,&b,&a);
  }while(!(i==0 || j==0)); /*break if one of the values are true*/
}


void tracebackOverlap(DYNAM *p_dynam)
{
  /**
     finding traceback for global alignment
   */
  int i,j,a=0,b=0;
  i=p_dynam->si;
  j=p_dynam->sj;
  /*dp the backtrack*/
  do{
    tracer(p_dynam,&i,&j,&b,&a);
  }while(!(i==0 || j==0)); /*break if one of the values are true*/
  p_dynam->ac1[a++]=0;
  p_dynam->ac2[b++]=0;
  revstring(p_dynam->ac1);
  revstring(p_dynam->ac2);
}


void tracebackPOZITIVGlobal(DYNAM *p_dynam)
{
  /**
     finding traceback for global alignment
   */
  int i,j,a=0,b=0;
  i=p_dynam->si; /*get starting postions*/
  j=p_dynam->sj; /*get starting postions*/
  do{
    tracerPOZITIV(p_dynam,&i,&j,&b,&a);
  }while(!(i==0 && j==0));
}

void tracebackGlobal(DYNAM *p_dynam)
{
  /**
     finding traceback for global alignment
   */
  int i,j,a=0,b=0;
  i=p_dynam->si; /*get starting postions*/
  j=p_dynam->sj; /*get starting postions*/
  do{
    tracer(p_dynam,&i,&j,&b,&a);
  }while(!(i==0 && j==0));
  p_dynam->ac1[a++]=0;
  p_dynam->ac2[b++]=0;
  revstring(p_dynam->ac1);
  revstring(p_dynam->ac2);
}

/*local traceback for POZITIVE SCORE*/
void tracebackPOZITIVLocal(DYNAM *p_dynam)
{
  int i,j,a=0,b=0;
  i=p_dynam->si; /*get starting postions*/
  j=p_dynam->sj; /*get starting postions*/
  do{
    tracerPOZITIV(p_dynam,&i,&j,&b,&a);
  }while(fabs(maximumdouble(3,p_dynam->h[i][j],p_dynam->Ix[i][j],p_dynam->Iy[i][j])) > 1e-12);
}

void tracebackLocal(DYNAM *p_dynam)
{
  /**
     finding traceback for local alignment
   */
  int i,j,a=0,b=0;
  i=p_dynam->si; /*get starting postions*/
  j=p_dynam->sj; /*get starting postions*/
  do{
    tracer(p_dynam,&i,&j,&b,&a);
  }while( fabs(maximumdouble(3,p_dynam->h[i][j],p_dynam->Ix[i][j],p_dynam->Iy[i][j])) >1e-12);
  p_dynam->ac1[a++]=0;
  p_dynam->ac2[b++]=0;
  revstring(p_dynam->ac1);
  revstring(p_dynam->ac2);
}

void initLocalA(DYNAM *p_dynam)
{
  int i,j;
  for(i=0; i< p_dynam->lc1;i++)
    {
      for(j=0; j< p_dynam->lc2;j++)
	{
	  p_dynam->h[i][j]=0;
	  p_dynam->Ix[i][j]=MIN_DOUBLE;
	  p_dynam->Iy[i][j]=MIN_DOUBLE;
	}
    }
}

void localAlign(DYNAM *p_dynam)
{
  int i,j;
  int lc1,lc2;
  double score = MIN_DOUBLE;
  double tmp;

  initLocalA(p_dynam);
  lc1=p_dynam->lc1;
  lc2=p_dynam->lc2;
  for(i = 1; i< lc1;i++)
    {
      for(j =1; j<lc2;j++)
	{
	  /*The same conditions like for global, except for the h matrix.*/
	  p_dynam->h[i][j]= maximumdouble(2, 0.0 ,hglob(p_dynam,i,j));
	  p_dynam->Ix[i][j]= ixglob(p_dynam,i,j);
	  p_dynam->Iy[i][j]= iyglob(p_dynam,i,j);
	  settrace(p_dynam, i, j);
	}
    }

  /*find maximal score in matrices*/
  for(i = 1 ; i< lc1; i++)
    {
      for(j = 1 ; j < lc2 ; j++)
	{
	  tmp = maximumdouble(3,  p_dynam->h[i][j] , p_dynam->Ix[i][j] , p_dynam->Iy[i][j] );
	  if(tmp>score)
	    {
	      p_dynam->si=i; /*set start of traceback*/
	      p_dynam->sj=j; /*set start of traceback*/
	      score=tmp;
	    }
	}
    }
  p_dynam->score = score;
}

void overlapAlign(DYNAM *p_dynam)
{
  int i , j , lc1 , lc2;
  double score=MIN_DOUBLE;
  double tmp;
  initLocalA(p_dynam);
  lc1=p_dynam->lc1;
  lc2=p_dynam->lc2;
  for(i = 1; i< lc1;i++)
    {
      for(j =1; j<lc2;j++)
	{
	  /*The same conditions like for global but different initialization*/
	  p_dynam->h[i][j]= hglob(p_dynam,i,j);
	  p_dynam->Ix[i][j]= ixglob(p_dynam,i,j);
	  p_dynam->Iy[i][j]= iyglob(p_dynam,i,j);
	  settrace(p_dynam,i,j);
	}
    }

  /*find maximal score in last row or last column*/
  for(i = 1 ; i< lc1; i++)
    {
      tmp = maximumdouble(3 , p_dynam->h[i][lc2-1] , p_dynam->Ix[i][lc2-1] , p_dynam->Iy[i][lc2-1]);
      if(tmp > score)
	{
	  /*set traceback start*/
	  p_dynam->si=i;
	  p_dynam->sj=lc2-1;
	  score=tmp;
	}
    }
  for(i=1 ; i< lc2; i++)
    {
      tmp = maximumdouble(3 ,  p_dynam->h[lc1-1][i] , p_dynam->Ix[lc1-1][i] , p_dynam->Iy[lc1-1][i]);
      if(tmp>score)
	{
	  /*set traceback start*/
	  p_dynam->si=lc1-1;
	  p_dynam->sj=i;
	  score=score;
	}
    }
  p_dynam->score = score;
}

/*Computes the score of the aligned sequence*/

double getselfalign(double **sim, const char *ac)
{

  int i;
  double score=0;
  int lc = strlen(ac);
  for(i=0;i<lc;i++)
    {
      if(ac[i]!='-') /*do not count matches*/
	{
	    score += sim[char2AA(ac[i])][char2AA(ac[i])];
	}
    }
  return score;
}



/*
Counts how many amino acids.
This is a count of the number of positions over the length of the alignment where >= 51% of the residues or bases at that position are similar. 
Any two residues or bases are defined as similar when they have positive comparisons (as defined by the comparison matrix being used in the alignment algorithm). 
*/
int getalignsimilarity(double **sim, char *ac1, char *ac2)
{
  int i; 
  int score = 0 ;
  int lc1 = strlen(ac1) ;
  int lc2 = strlen(ac2) ;
  if(lc1!=lc2)
    return -1; /** error handling */
  for(i=0; i<lc1; i++)
    {
      if(ac1[i]!='-' && ac2[i] != '-') /**if one is gap they cant be similar **/
	{
	  if(sim[char2AA(ac1[i])][char2AA(ac2[i])]>0)
	    {
	      score++;
	    }
	}
    }
  return score;
}

int getalignidentity(char *ac1, char *ac2)
{
  int i, score=0;
  int lc1 = strlen(ac1);
  int lc2 = strlen(ac2);
  if(lc1!=lc2)
    return -1;
  for(i=0; i<lc1; i++)
    {
      if(ac1[i]!='-' && ac2[i] != '-') /**if one is gap they cant be identical **/
	{
	  if(ac1[i]==ac2[i])
	    {
	      score++;
	    }
	}
    }
  return score;
}


void identSimilarScore(
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
    )
{
/*  printf("all %d",*all);
    printf("score %f\n",*score);
    printf("mdim %d\n",*mdim);
*/
  DYNAM test;
  /*set chars*/
  int malloccount=0;
  test.c1 = cc1;
  test.c2 = cc2;
  test.ac1 = malloc((strlen(cc1) + strlen(cc2))*sizeof(char));
  malloccount++;
  test.ac2 =  malloc((strlen(cc1) + strlen(cc2))*sizeof(char));
  malloccount++;
  /*set dimensions of the dynamic programming matrices*/
  test.lc1 = strlen(test.c1)+1; /*length of the numerical representation of the string.*/
  test.lc2 = strlen(test.c2)+1;
  /*set dimension of the similarity matrix*/
  test.dim = mdim;
  /*set gap opening cost and gap extend*/
  test.delta = mdelta;
  test.gapext = gapext;
  /*alocate memmory*/
  allocDYNAM(&test);
    
    /* translate char to similartity matrix indices */
    if(-1==seq2index(test.nc1, test.c1 , test.lc1))
      {
	strcpy(errormsg,"Nonstandard amino acid");
	return;
      }
    if(-1==seq2index(test.nc2, test.c2 , test.lc2))
      {
	strcpy(errormsg,"Nonstandard amino acid");
	return;
      };

    /*copy stuff into similarity matrix*/
    initsim(test.sim,test.dim,simm);
    /*        printstrut(&test);*/
    if(strcmp(align_type,"global")==0)
      {
	globalAlign(&test);
	tracebackGlobal(&test);
      }
    else if(strcmp(align_type,"local")==0)
      {
	localAlign(&test);
	tracebackLocal(&test);
      }
    else if(strcmp(align_type,"overlap")==0)
      {
	overlapAlign(&test);
	tracebackOverlap(&test);
      }
    else
      {
	printf("No such type of alignment : [ %s ] \n",align_type);
	return;
      }
    if(strcmp(score_t,"identity")==0)
      {
	if(test.lc1<=test.lc2)
	  {
	    *score = ((double)getalignidentity(test.ac1,test.ac2))/((double)(strlen(test.c1)));
	  }
	else
	  *score = ((double)getalignidentity(test.ac1,test.ac2))/((double)(strlen(test.c2)));
      }
    else
      {
	if(test.lc1<=test.lc2)
	  {
	    *score = ((double)getalignsimilarity(test.sim, test.ac1 , test.ac2 ))/((double)(strlen(test.c1)));
	  }
	else
	  *score = ((double)getalignsimilarity(test.sim, test.ac1 , test.ac2 ))/((double)(strlen(test.c2)));
      }
    free(test.ac1);
    malloccount--;
    free(test.ac2);
    malloccount--;
    freeDYNAM(&test);
}

void globalB(
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
    )
{
/*  printf("all %d",*all);
    printf("score %f\n",*score);
    printf("mdim %d\n",*mdim);
*/
  DYNAM test;
  /*set chars*/
    test.c1 = cc1;
    test.c2 = cc2;
    test.ac1 = al1;
    test.ac2 = al2;
    /*set dimensions of the dynamic programming matrices*/
    test.lc1 = strlen(test.c1)+1; /*length of the numerical representation of the string.*/
    test.lc2 = strlen(test.c2)+1;
    /*set dimension of the similarity matrix*/
    test.dim = mdim;
    /*set gap opening cost and gap extend*/
    test.delta = mdelta;
    test.gapext = gapext;
    /*alocate memmory*/
    allocDYNAM(&test);
    /* translate char to similartity matrix indices */
    if(-1==seq2index(test.nc1, test.c1 , test.lc1))
      {
	strcpy(errormsg,"Nonstandard amino acid");
	return;
      }
    if(-1==seq2index(test.nc2, test.c2 , test.lc2))
      {
	strcpy(errormsg,"Nonstandard amino acid");
	return;
      };

    /*copy stuff into similarity matrix*/
    initsim(test.sim,test.dim,simm);
    /*        printstrut(&test);*/
    if(strcmp(type,"global")==0)
      {
	globalAlign(&test);
	/*	printmat(test.h,test.lc1,test.lc2);*/
	/*	printmat(test.Ix,test.lc1,test.lc2);*/
	/*	printmat(test.Iy,test.lc1,test.lc2);*/
	tracebackGlobal(&test);
      }
    else if(strcmp(type,"local")==0)
      {
	localAlign(&test);
	tracebackLocal(&test);
      }
    else if(strcmp(type,"overlap")==0)
      {
	overlapAlign(&test);
	tracebackOverlap(&test);
      }
    else
      {
	printf("No such type of alignment : [ %s ] \n",type);
	return;
      }
    *score = test.score;
    al1=(test.ac1);
    al2=(test.ac2);
    *selfscore1 =  getselfalign(test.sim , test.c1);
    *selfscore2 =  getselfalign(test.sim , test.c2);
    *all = strlen(test.ac1);
    *alignsimilarity =  getalignsimilarity(test.sim, test.ac1 , test.ac2 );
    *identity = getalignidentity(test.ac1,test.ac2);
    freeDYNAM(&test);
}

/*Returns Smith Watermann Score*/
void globalScore
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
)
{
  DYNAM test;
  double selfscore1,selfscore2;
  /*set chars*/
  test.c1=cc1;
  test.c2=cc2;
  /*set dimensions of the dynamic programming matrices*/
  test.lc1 = strlen(test.c1)+1; /*length of the numerical representation of the string.*/
  test.lc2 = strlen(test.c2)+1;
  /*set dimension of the similarity matrix*/
  test.dim = mdim;
  /*set gap opening cost and gap extend*/
  test.delta = mdelta;
  test.gapext = gapext;
  /*alocate memmory*/
  allocDYNAM(&test);
  /* translate char to similartity matrix indices */
  if(-1==seq2index(test.nc1, test.c1 , test.lc1)){
    strcpy(errormsg,"Nonstandard amino acid");
    return;
  }
  if(-1==seq2index(test.nc2, test.c2 , test.lc2)){
    strcpy(errormsg,"Nonstandard amino acid");
    return;
  };
  /*copy stuff into similarity matrix*/
  initsim(test.sim,test.dim,simm);
  if(strcmp(alig_type,"global")==0)
    {
      globalAlign(&test);
    }
  else if(strcmp(alig_type,"local")==0)
    {
      localAlign(&test);
    }
  else if(strcmp(alig_type,"overlap")==0)
    {
      overlapAlign(&test);
    }
  else
    {
      printf("No such type of alignment : %s !\n",alig_type);
    }

  if(strcmp(score_t,"scoreN")==0)
    {
      selfscore1 =  getselfalign(test.sim , test.c1);
      selfscore2 =  getselfalign(test.sim , test.c2);
      if(test.score<=0)
	{
	  *score=0;
	}
      else
	{
	  if(selfscore1<selfscore2)
	    {
	      *score = test.score/selfscore1;
	    }
	  else
	    {
	      *score = test.score/selfscore2;
	    }
	}
    }
  else
    {
      *score = test.score;
    }
  freeDYNAM(&test);

}

double pozitive(DYNAM *dynam_p)
{

  /*now its time to generate 2 strings without gaps.
  compute mu.*/
  int il, jl; 
  double mu = 0; /*sum of s_ij */
  double mu_2=0; /*sum of squared s_ij */
  double s_ij; /*s_ij entry*/
  double mu_i = 0; /*row sum*/
  double mu_i_2 = 0; /*squared row sum*/
  double mu_j = 0; /*column sum*/
  double mu_j_2 = 0; /*squared col sum*/
  double sigma = 0;
  double zscore = 0;
  for(il = 0; il < dynam_p->ungl; il++)
    {
      for(jl = 0; jl< dynam_p->ungl; jl++)
	{
	  s_ij = dynam_p->sim[dynam_p->ung1[il]][dynam_p->ung2[jl]];
	  mu +=  s_ij;
	  mu_2 +=  (s_ij * s_ij);
	  mu_j +=  s_ij;
	}
      mu_i_2 +=  (mu_j * mu_j);
      mu_j = 0;
    }

  for(jl=0; jl< dynam_p->ungl; jl++)
    {
      for(il=0; il< dynam_p->ungl; il++)
	{
	  s_ij = dynam_p->sim[dynam_p->ung1[il]][dynam_p->ung2[jl]];
	  mu_i += s_ij;
	}
      mu_j_2 += (mu_i*mu_i);
      mu_i = 0;
    }

  /*  printf("mu_paths = %f \n",mu/test.ungl);*/
  /*compute sigma.*/

  sigma = mu_2/dynam_p->ungl 
    + (mu*mu - mu_i_2 - mu_j_2 + mu_2)/(dynam_p->ungl*(dynam_p->ungl-1)) 
    - (mu/dynam_p->ungl)*(mu/dynam_p->ungl);

  zscore = (dynam_p->score-(mu/dynam_p->ungl))/sigma;
  /*  printf("SA score = %f\n",test.score);*/
  /*  printf("zscore = %f\n",zscore);*/
  return(zscore);
}


void pozitiveScore(
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
    )
{
  DYNAM test;
  int malloccount=0;
  malloccount++;
  test.c1 = cc1;
  test.c2 = cc2;
  /*set dimensions of the dynamic programming matrices*/
  test.lc1 = strlen(test.c1)+1; /*length of the numerical representation of the string.*/
  test.lc2 = strlen(test.c2)+1;
  /*set dimension of the similarity matrix*/
  test.dim = mdim;
  /*set gap opening cost and gap extend*/
  test.delta = mdelta;
  test.gapext = gapext;
  /*alocate memmory*/
  allocDYNAM(&test);
  /*pozitive*/
  test.pozitiv=1;
  test.ungl=0;
  /*allocate space for inidices of ungapped strings.*/
  if(test.lc1<test.lc2)
    {
      test.ung1 = calloc(test.lc1 , sizeof(int));
      malloccount++;
      test.ung2 = calloc(test.lc1 , sizeof(int));
      malloccount++;
    }
  else
    {
      test.ung1 = calloc(test.lc2 , sizeof(int));
      malloccount++;
      test.ung2 = calloc(test.lc2 , sizeof(int));
      malloccount++;
    }
  /* translate char to similartity matrix indices */
  if(-1==seq2index(test.nc1, test.c1 , test.lc1))
    {
      strcpy(errormsg,"Nonstandard amino acid");
      return;
    }
  if(-1==seq2index(test.nc2, test.c2 , test.lc2))
    {
      strcpy(errormsg,"Nonstandard amino acid");
      return;
    };
  
  /*copy stuff into similarity matrix*/
  initsim(test.sim,test.dim,simm);
  
  /*printstrut(&test);*/
  if(strcmp(align_type,"global")==0)
    {
      globalAlign(&test);
      tracebackPOZITIVGlobal(&test);
    }
  else if(strcmp(align_type,"local")==0)
    {
      localAlign(&test);
      tracebackPOZITIVLocal(&test);
    }
  else if(strcmp(align_type,"overlap")==0)
    {
      overlapAlign(&test);
      tracebackPOZITIVOverlap(&test);
    }
  else
    {
      printf("No such type of alignment : [ %s ] \n",align_type);
      return;
    }
  *score=pozitive(&test);
  free(test.ung1);
  malloccount--;
  free(test.ung2);
  malloccount--;
  if(malloccount==0)
    error("The Devil are in your Code!\n");
  freeDYNAM(&test);
}
