/**********************************************************************
 * This program reports the following parameters calculated from DDA: *
 *   Cabs - absorption cross section (m^2)			      *
 *   Csca - scattering cross section (m^2)			      *
 *   Cbsc - backscattering cross section (m^2)			      *
 *   g    - asymmetry parameter	(0 - 1)				      *
 *   P(Q) - phase function from 0 (foreward) to 180 (backward) degree *
 *          for every 5 degrees. Normalized, so that integral of      *
 *         1/2*integral[P(Q)dcosQ] = 1				      *
 * For randomally oriented ice particles of the following shapes:     *
 *  SHAPE-ID		SHAPE#		SHAPE			      * 
 *   hexl		  0		long hexagonal column l/d=4   *
 *   hexs		  1		short hexagonal col   l/d=2   *
 *   hexb		  2             block hex col	      l/d=1   *
 *   hexf		  3             thick hex plate	      l/d=0.2 *
 *   hexp		  4		thin  hex plate	     l/d=0.05 *
 *   ros3		  5		3-bullet rosette	      *
 *   ros4		  6 		4-bullet rosette	      *
 *   ros5		  7		5-bullet rosette	      *
 *   ros6		  8		6-bullet rosette	      *
 *   sstr		  9 		sector-like snowflake	      *
 *   sden		 10		dendrite snowflake	      *
 * For frequency range:  3   - 340 GHz				      *
 * For temperature range: 0 - -40 C				      *
 * References:							      *
 * Liu,2004: JAS, 15, 2441-2456.				      *
 * Liu, 2008: BAMS, 89, 1563-1570.				      *
 *								      *
 * call by fortran:						      *
 *  call scatdb(f,t,nshape,dmax,cabs,csca,cbsc,g,p,re,iret,is_loaded) *
 * call by c:							      *
 *  iret = scatdb(f,t,nshape,dmax,&cabs,&csca,&cbsc,&g,&p,&re,        *
 *  			&is_loaded)        			      *
 * Input:							      *
 *  f - frequency in GHz	real/float			      *
 *            Range: 3    - 340 GHz				      *
 *  t - temperature in K	real/float			      *
 *            Range: 273.15 - 233.15 K			      	      *
 *  nshape - shape#	(see above)	integer/int		      *
 *            Range: 0 - 10					      *
 *  dmax - particle's maximum dimension (micron) real/float	      *
 *            Range is shape dependent as follows		      *
 *             SHAPE-ID	SHAPE#	dmax-Range(um)	re-Range(um)          *
 *		hexl	0	121 - 4835	25 - 1000	      *
 *              hexs	1	83  - 3304	25 - 1000	      *
 * 		hexb	2	66  - 2632	25 - 1000	      *
 *  		hexf	3	81  - 3246	25 - 1000	      *
 *		hexp	4	127 - 5059	25 - 1000	      *	
 *		ros3	5	50  - 10000	19 - 1086	      *
 *		ros4	6	50  - 10000	19 - 984	      *
 *		ros5	7	50  - 10000	21 - 1058	      *
 *		ros6	8	50  - 10000	21 - 1123	      *
 *		sstr	9	50  - 10000	25 - 672	      *
 *              sden	10	75  - 12454     33 - 838	      *
 *     Note: While the program will return values when dmax is out of *
 *           range, the resultant values are probably wrong because   *
 *           they are EXTRAPOLATED, particularly for phase function   *
 *           you may see negative values, obviously wrong)	      *	
 * Output:							      *
 *  cabs - absorption cross section m^2  real/float		      *		
 *  csca - scattering cross section m^2  real/float		      *
 *  cbsc - backscattering cross section m^2 real/float		      *
 *  g - asymmetry parameter		real/float		      *
 *  p - phase function, 37 elements (0 to 180deg every 5deg)          *
 *  		real/float					      *
 *  re - sphere-equivalent radius (assume density=0.916)              *
 *  		real/float					      *
 *  iret: return status	integer/int				      *
 *  	0	- 	success					      *
 *     xx1(2)   - 	frequency higher (lower) than 340 ( 3  ) GHz  *
 *     x1(2)x   -	temperature warmer (colder) than 273.15       *
 *  					(233.15) K		      *
 *     1(2)xx  -	dmax larger (smaller) than that in database   *
 *   Extrapolation will be used for 10>iret>0                         *
 *      1000	- 	can not find suitable ice shape		      *
 *      	        No valid return values			      * 
 *      2000    -       cannot find scat_db2.dda file                  *
 *								      *
 * In&Out:							      *
 *  is_loaded - indicator whether dda_database is loaded to memory,   *
 *              assign 0 at the beginning of the program and Never    *
 *              change it again. integer/int                          *
 *   NOTE:							      *
 *   The scat_db2.dda file contains the data computed by DDA. It needs*
 *   to be located in the directory indicated by the SCATDB_DATA     *
 *   environmental variable. If the variable is not set, the program *
 *   assumes the data file is in the current (run-time) directory.   *
 *   To set SCATDB_DATA, under csh/tcsh:                             *
 *          setenv SCATDB_DATA (the-directory)			      *
 *   under sh/bash:						      *
 *           export SCATDB_DATA=(the-directory)                       *
 *                                                                    *
 * ********************************************************************/

/* Ver2.1  2010.3.22 add 5 GHz*/
/* Ver2.2  2010.7.9  add 3, 9, 24.1*/
/* Ver2.3  correct minor error near 1000 um*/
/* Ver2.4  2011.3.22 add 10,15,19,60,70,80,90 GHz*/
/* Ver2.5  2011.4.9 add 50 GHz*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define	TRUE	1
#define FALSE	0

#define NSHAP		11   /* number of shapes */
#define NTEMP		5    /* number of temps  */
#define NFREQ		22   /* number of frequencies */
#define NSIZE		20   /* max number of sizes */
#define NQ		37   /* number of anges in PF */

typedef char bool;
int little_endian();
void reverse(void*,int);
void scatdb_(float*,float*,int*,float*,float*,float*,float*,float*,float*,float*,int*,int*);
int scatdb(float,float,int,float,float*,float*,float*,float*,float*,float*,int*);
float linear_interp(float,float,float,float,float);

void scatdb_(float *f,float *t,int *nshape,float *dmax,float *cabs,float *csca,float *cbsc,float *g,float *p,float *re,int *iret, int *is_loaded) {
/* wrapper for fortran interface */
	float f1, t1, dmax1;
	int n1;
	f1 = *f;
	t1 = *t;
	n1 = *nshape;
	dmax1 = *dmax;
	*iret = scatdb(f1,t1,n1,dmax1,cabs,csca,cbsc,g,p,re,is_loaded);
}

float linear_interp(float x,float x1,float x2,float y1,float y2) {
/* linear interpolation */
	float y;
	y=y1+(y2-y1)/(x2-x1)*(x-x1);
	return y;
}

int little_endian() {
	const int endian_test = 1;
	bool check_endian;
	memcpy(&check_endian, &endian_test, 1);
	if(check_endian) return TRUE;
	else return FALSE;
}

void reverse(void *v, int n) {
/* Reverse n bytes */
    char c;
    char *p1, *p2;

    for (p1 = v, p2 = p1 + n - 1; p1 < p2; c = *p1, *p1++ = *p2, *p2-- = c)
        ;
}

int scatdb(float f, float t, int nshape, float dmax, float *cabs, float *csca, float *cbsc, float *g, float *p, float *re, int *is_loaded) {
       	char scat_db_dir[180]={"."};
	const char scat_db_fn[]={"scat_db2.dda"};
	FILE *fp;
	static float fs[NFREQ],ts[NTEMP],szs[NSIZE][NSHAP],abss[NFREQ][NTEMP][NSHAP][NSIZE],scas[NFREQ][NTEMP][NSHAP][NSIZE],
		     bscs[NFREQ][NTEMP][NSHAP][NSIZE],gs[NFREQ][NTEMP][NSHAP][NSIZE],reff[NSIZE][NSHAP],pqs[NFREQ][NTEMP][NSHAP][NSIZE][NQ];
	static int shs[NSHAP],mf,mt,msh,msz[NSHAP];
	int iret=0,it1=-1,it2=-1,if1=-1,if2=-1,ir1=-1,ir2=-1,big=TRUE;
	float x1,x2,y1,y2,a1,b1,c1,d1;
	int i,j,k,m,n;
	char *db_dir, wk[180];
	if(*is_loaded != TRUE) {
		if((db_dir=getenv("SCATDB_DATA")))
		  strcpy(scat_db_dir,db_dir);
		sprintf(wk,"%s/%s",scat_db_dir,scat_db_fn);
		if(!(fp=fopen(wk,"rb")))	{
			fprintf(stderr,"Cannot find scat_db2.dda file. %s\n",wk);
			return(iret=2000);
		}
	
		if(little_endian() == TRUE) big = FALSE;	

		fread(&msh,4,1,fp); 
		if(big) reverse(&msh,4);
		fread(&mf,4,1,fp); 
		if(big) reverse(&mf,4);
		fread(&mt,4,1,fp); 
		if(big) reverse(&mt,4);

		for(i=0;i<msh;i++) {
		        fread(&shs[i],4,1,fp); 
			if(big) reverse(&shs[i],4);	
			fread(&msz[i],4,1,fp); 
			if(big) reverse(&msz[i],4); 
				for(j=0;j<msz[i];j++) {
					fread(&szs[j][shs[i]],4,1,fp);
					if(big) reverse(&szs[j][shs[i]],4);
					fread(&reff[j][shs[i]],4,1,fp);
					if(big) reverse(&reff[j][shs[i]],4);
					for(m=0;m<mf;m++) {
						fread(&fs[m],4,1,fp);
						if(big) reverse(&fs[m],4);
						for(n=0;n<mt;n++) {
							fread(&ts[n],4,1,fp);
							if(big) reverse(&ts[n],4);
							fread(&abss[m][n][shs[i]][j],4,1,fp);
							if(big) reverse(&abss[m][n][shs[i]][j],4);
							fread(&scas[m][n][shs[i]][j],4,1,fp);
							if(big) reverse(&scas[m][n][shs[i]][j],4);
							fread(&bscs[m][n][shs[i]][j],4,1,fp);
							if(big) reverse(&bscs[m][n][shs[i]][j],4);
 							fread(&gs[m][n][shs[i]][j],4,1,fp);
							if(big) reverse(&gs[m][n][shs[i]][j],4);
							for(k=0;k<NQ;k++) {
								fread(&pqs[m][n][shs[i]][j][k],4,1,fp);
								if(big) reverse(&pqs[m][n][shs[i]][j][k],4);
							}
						}
					}
				}
		}
		fclose(fp);
		*is_loaded = TRUE;
	}

	if((nshape<0) || (nshape>msh-1)) return(iret=1000);

	for(i=0;i<mf-1;i++) 
		if((f>=fs[i]) && (f<=fs[i+1])) {
			if1=i;
			if2=i+1;
			break;
		}
	if(if1 == -1) {
		if(f<fs[0]) {if1=0;if2=1;iret += 2;}
		else {if1=mf-2;if2=mf-1;iret += 1;}
	}

	for(i=0;i<mt-1;i++) 
		if((t<=ts[i]) && (t>=ts[i+1])) {
			it1=i;
			it2=i+1;
			break;
		}

	if(it1 == -1) {
		if(t>ts[0]) {it1=0;it2=1;iret += 10;}
		else {it1=mt-2;it2=mt-1;iret += 20;}
	}

	for(i=0;i<msz[nshape]-1;i++) 
		if((dmax>=szs[i][nshape]) && (dmax<=szs[i+1][nshape])) {
			ir1=i;
			ir2=i+1;
			break;
		}

	if(ir1 == -1) { 
		if(dmax<szs[0][nshape]) {ir1=0;ir2=1;iret += 200;}
		else {ir1=msz[nshape]-2;ir2=msz[nshape]-1;iret += 100;}
	}

/* sphere equavalent radius */
	x1=szs[ir1][nshape];y1=reff[ir1][nshape];x2=szs[ir2][nshape];y2=reff[ir2][nshape];
        *re=linear_interp(dmax,x1,x2,y1,y2);

/* absorption */
	x1=fs[if1]; x2=fs[if2];
	y1=abss[if1][it1][nshape][ir1]; y2=abss[if2][it1][nshape][ir1];
	a1=linear_interp(f,x1,x2,y1,y2);
	y1=abss[if1][it1][nshape][ir2]; y2=abss[if2][it1][nshape][ir2];
	b1=linear_interp(f,x1,x2,y1,y2);
	x1=szs[ir1][nshape]; x2=szs[ir2][nshape]; y1=a1; y2=b1;
	c1=linear_interp(dmax,x1,x2,y1,y2);
	x1=fs[if1]; x2=fs[if2];
	y1=abss[if1][it2][nshape][ir1]; y2=abss[if2][it2][nshape][ir1];
	a1=linear_interp(f,x1,x2,y1,y2);
	y1=abss[if1][it2][nshape][ir2]; y2=abss[if2][it2][nshape][ir2];
	b1=linear_interp(f,x1,x2,y1,y2);
	x1=szs[ir1][nshape]; x2=szs[ir2][nshape]; y1=a1; y2=b1;
	d1=linear_interp(dmax,x1,x2,y1,y2);
	x1=ts[it1]; x2=ts[it2]; y1=c1; y2=d1;
	*cabs=linear_interp(t,x1,x2,y1,y2);

/* scattering */
	x1=fs[if1]; x2=fs[if2];
	y1=scas[if1][it1][nshape][ir1]; y2=scas[if2][it1][nshape][ir1];
	a1=linear_interp(f,x1,x2,y1,y2);
	y1=scas[if1][it1][nshape][ir2]; y2=scas[if2][it1][nshape][ir2];
	b1=linear_interp(f,x1,x2,y1,y2);
	x1=szs[ir1][nshape]; x2=szs[ir2][nshape]; y1=a1; y2=b1;
	c1=linear_interp(dmax,x1,x2,y1,y2);
	x1=fs[if1]; x2=fs[if2];
	y1=scas[if1][it2][nshape][ir1]; y2=scas[if2][it2][nshape][ir1];
	a1=linear_interp(f,x1,x2,y1,y2);
	y1=scas[if1][it2][nshape][ir2]; y2=scas[if2][it2][nshape][ir2];
	b1=linear_interp(f,x1,x2,y1,y2);
	x1=szs[ir1][nshape]; x2=szs[ir2][nshape]; y1=a1; y2=b1;
	d1=linear_interp(dmax,x1,x2,y1,y2);
	x1=ts[it1]; x2=ts[it2]; y1=c1; y2=d1;
	*csca=linear_interp(t,x1,x2,y1,y2);	

/* backscattering */
	x1=fs[if1]; x2=fs[if2];
	y1=bscs[if1][it1][nshape][ir1]; y2=bscs[if2][it1][nshape][ir1];
	a1=linear_interp(f,x1,x2,y1,y2);
	y1=bscs[if1][it1][nshape][ir2]; y2=bscs[if2][it1][nshape][ir2];
	b1=linear_interp(f,x1,x2,y1,y2);
	x1=szs[ir1][nshape]; x2=szs[ir2][nshape]; y1=a1; y2=b1;
	c1=linear_interp(dmax,x1,x2,y1,y2);
	x1=fs[if1]; x2=fs[if2];
	y1=bscs[if1][it2][nshape][ir1]; y2=bscs[if2][it2][nshape][ir1];
	a1=linear_interp(f,x1,x2,y1,y2);
	y1=bscs[if1][it2][nshape][ir2]; y2=bscs[if2][it2][nshape][ir2];
	b1=linear_interp(f,x1,x2,y1,y2);
	x1=szs[ir1][nshape]; x2=szs[ir2][nshape]; y1=a1; y2=b1;
	d1=linear_interp(dmax,x1,x2,y1,y2);
	x1=ts[it1]; x2=ts[it2]; y1=c1; y2=d1;
	*cbsc=linear_interp(t,x1,x2,y1,y2);

/* asymmetry parameter */
	x1=fs[if1]; x2=fs[if2];
	y1=gs[if1][it1][nshape][ir1]; y2=gs[if2][it1][nshape][ir1];
	a1=linear_interp(f,x1,x2,y1,y2);
	y1=gs[if1][it1][nshape][ir2]; y2=gs[if2][it1][nshape][ir2];
	b1=linear_interp(f,x1,x2,y1,y2);
	x1=szs[ir1][nshape]; x2=szs[ir2][nshape]; y1=a1; y2=b1;
	c1=linear_interp(dmax,x1,x2,y1,y2);
	x1=fs[if1]; x2=fs[if2];
	y1=gs[if1][it2][nshape][ir1]; y2=gs[if2][it2][nshape][ir1];
	a1=linear_interp(f,x1,x2,y1,y2);
	y1=gs[if1][it2][nshape][ir2]; y2=gs[if2][it2][nshape][ir2];
	b1=linear_interp(f,x1,x2,y1,y2);
	x1=szs[ir1][nshape]; x2=szs[ir2][nshape]; y1=a1; y2=b1;
	d1=linear_interp(dmax,x1,x2,y1,y2);
	x1=ts[it1]; x2=ts[it2]; y1=c1; y2=d1;
	*g=linear_interp(t,x1,x2,y1,y2);

/* phase function */
	for (i=0;i<NQ;i++) {
	x1=fs[if1]; x2=fs[if2];
	y1=pqs[if1][it1][nshape][ir1][i]; y2=pqs[if2][it1][nshape][ir1][i];
	a1=linear_interp(f,x1,x2,y1,y2);
	y1=pqs[if1][it1][nshape][ir2][i]; y2=pqs[if2][it1][nshape][ir2][i];
	b1=linear_interp(f,x1,x2,y1,y2);
	x1=szs[ir1][nshape]; x2=szs[ir2][nshape]; y1=a1; y2=b1;
	c1=linear_interp(dmax,x1,x2,y1,y2);
	x1=fs[if1]; x2=fs[if2];
	y1=pqs[if1][it2][nshape][ir1][i]; y2=pqs[if2][it2][nshape][ir1][i];
	a1=linear_interp(f,x1,x2,y1,y2);
	y1=pqs[if1][it2][nshape][ir2][i]; y2=pqs[if2][it2][nshape][ir2][i];
	b1=linear_interp(f,x1,x2,y1,y2);
	x1=szs[ir1][nshape]; x2=szs[ir2][nshape]; y1=a1; y2=b1;
	d1=linear_interp(dmax,x1,x2,y1,y2);
	x1=ts[it1]; x2=ts[it2]; y1=c1; y2=d1;
	p[i]=linear_interp(t,x1,x2,y1,y2);
	}

	return iret;
}
