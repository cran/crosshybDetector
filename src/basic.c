/**
* translates the sequence string to a index array.
*/
#include "score.h"

int  seq2index(int *nc, const char *c, int lc)
{
    int i;
    int tmp;
    nc[0]=lc;
    for(i = 1; i<lc; i++)
    {
      tmp=char2AA(c[i-1]);
      if(tmp ==-1)
	{
	  return -1;
	}else
	{
	  nc[i] = tmp ;	// copy to int array and cast
	}
    }
    return 1;
}
    
/**
copies the weight from simm to dimm
*/

void initsim(double **sim, int dim, double *simm)
{
    int i,j;
    for(i =0; i<dim; i++)
    {
        for(j=0; j<dim; j++)
            {
                        sim[i][j] = simm[j+i*dim];
            }
    }
}

void printsim(double **sim, int test)
{
    int j,i;
    for(i=0;i<test; i++)
    {
        for(j=0;j<test; j++)
            {
                printf("\t%f",sim[i][j]);    
            }
            printf("\n");
    }
}

void printmat(double **sim, int n, int m)
{
    int j,i;
    printf("\n");

    for(i=0;i<n; i++)
    {
        for(j=0;j<m; j++)
            {
                printf(" %3f",sim[i][j]);    
            }
            printf("\n");
    }
    printf("\n\n");
}

/**allocates storage for a matrix and initializes with zeros*/

void allocmatrixint(int ***mat, int x, int y)
{
  int i;
  int **mati;
  if(x<=0){printf(" Refusing to claim no (or even less memmory ... !");exit(1);}
  if(y<=0){printf(" Refusing to claim no (or even less memmory ... !");exit(1);}
  
  mati = malloc(x * sizeof(int*));
  for(i = 0; i< x ; i++)
    {
      mati[i]=  calloc( y, sizeof( int ) );
    }
  *mat = mati;
}


void allocmatrixdouble(double ***mat, int x, int y)
{
    int i;
    double **mati;
    if(x<=0){printf(" Refusing to claim no (or even less memmory ... !");exit(1);}
    if(y<=0){printf(" Refusing to claim no (or even less memmory ... !");exit(1);}

    mati = malloc(x * sizeof(double*));
    for(i = 0; i< x ; i++)
    {
      mati[i]=  calloc( y ,sizeof( double ) );
    }
    *mat = mati;
}


/*free the memmory allocated by the matrix allocmatrixint*/
void freematrixint(int ***mat,int x, int y)
{
  int i;
  for(i=0; i<x ; i++)
    {
      free((void *)(*mat)[i]);
    }
  free((void *)(*mat));
}

void freematrixdouble(double ***mat,int x, int y)
{
  int i;
  for(i=0; i<x ; i++)
    {
      free((void *)(*mat)[i]);
    }
  free((void *)(*mat));
}

/**allocates storage for a matrix and initializes with zeros*/

void allocmatrixshort(short ***mat, int x, int y)
{
    int i,j;
    short **mati;
    if(x<=0){printf(" Refusing to claim no (or even less memmory ... !");exit(1);}
    if(y<=0){printf(" Refusing to claim no (or even less memmory ... !");exit(1);}
    mati = malloc(x * sizeof(short*));
    for(i = 0; i< x ; i++)
    {
        mati[i]=  malloc( y * sizeof( short ) );
    }
    /*initialize matrix with zeros.*/
    for(i=0;i < x;i++)
    {
        for(j=0;j < y;j++)
        {
                mati[i][j]=0;
        }
    }
    *mat = mati;
}


/**
void merror(char * s) {
	fprintf(stderr,"%s\n",s);
	exit(1);
}
*/
int char2AA(char ch){
  switch (ch) {
  case 'A' : return 0;
  case 'R' : return 1;
  case 'N' : return 2;
  case 'D' : return 3;
  case 'C' : return 4;
  case 'Q' : return 5;
  case 'E' : return 6;
  case 'G' : return 7;
  case 'H' : return 8;
  case 'I' : return 9;
  case 'L' : return 10;
  case 'K' : return 11;
  case 'M' : return 12;
  case 'F' : return 13;
  case 'P' : return 14;
  case 'S' : return 15;
  case 'T' : return 16;
  case 'W' : return 17;
  case 'Y' : return 18;
  case 'V' : return 19;
  case 'B' : return 20;
  case 'Z' : return 21;
  case 'X' : return 22;
  case '*' : return 23;
  default  : 
    return -1;
  }
  return -1;      
}


int AA2char(int x){
  switch (x) {
  case 0 : return 'A';
  case 1 : return 'R';
  case 2 : return 'N';
  case 3 : return 'D';
  case 4 : return 'C';
  case 5 : return 'Q';
  case 6 : return 'E';
  case 7 : return 'G';
  case 8 : return 'H';
  case 9 : return 'I';
  case 10: return 'L';
  case 11: return 'K';
  case 12: return 'M';
  case 13: return 'F';
  case 14: return 'P';
  case 15: return 'S';
  case 16: return 'T';
  case 17: return 'W';
  case 18: return 'Y';
  case 19: return 'V';
  case 20: return 'B';
  case 21: return 'Z';
  case 22: return 'X';
  case 23: return '*';
  default: return -1;
  }
  return -1;
}


double maximumdouble3(double a, double b, double c )
{
  double res;
  if(a>b)
    res=a;
  else
    res=b;
  if(c>res)
    res=c;
  return(res);
}
double maximumdouble(int n,...)
{
    va_list ap;
    double ival, max = -100000.0;
    int i;
    va_start(ap,n);
    for(i = 0; i < n; i++)
    { 
        ival = va_arg(ap, double);
        if (ival > max) max = ival;
    }
    va_end(ap);
    return max;
}


int maximum(int n,...)
{
    va_list ap;
    int ival, max = -100000;
    int i;
    va_start(ap,n);
    for(i = 0; i < n; i++)
    { 
        ival = va_arg(ap, int);
        if (ival > max) max = ival;
    }
    va_end(ap);
    return max;
}

/*reverses string*/
void revstring(char *ch)
{
  int i,j,strl,lp;
  char c;
  strl = strlen(ch);
  lp=strl-1;
  j = ((int) lp/2)+1;
  for(i = 0 ; i<j ; i++)
    {
      c=(ch)[i];
      ch[i] = ch[lp-i];
      ch[lp-i] = c;

      /*printf("%d %c\n",i,(*ch)[i]);*/
    }
}
