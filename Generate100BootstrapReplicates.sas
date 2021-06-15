/*********** Investigator: Rebecca Brown ***********/
/*********** Statistician: Grisell Diaz-Ramirez *********** /
/*********** Date created: 2021.06.14 ***********/
/*********** Purpose: Create 100 boostrap samples with 32,486  Ids/sample (same as original data) ***********/

/* Check system options specified at SAS invocation */
proc options option=work; run;
proc options option=utilloc; run;
proc options group=(memory performance); run;

libname in 'path';
libname out 'path'; 
options nosource nonotes; /*nosource: suppress the listing of the SAS statements to the log, causes only source lines that contain errors to be written to the log*/
options errorabend; /*so the process stops as soon as an error is encountered*/

/********************************************************************************************************************************************************************************/
/*** Create input and multimodel data ***/

data input; set in.data1obPerId (keep=newid normwgt_198int RAESTRAT RAEHSAMP) ; run;
proc sort data=input; by RAESTRAT RAEHSAMP; run;
/*  32,486 observations and 4 variables */

data multimodel (keep=newid wavenum weightbl RAESTRAT RAEHSAMP ageint RAGENDER bathdif beddif dressdif eatdif toiletdif walkdif phonedif moneydif meddif shopdif mealdif);
 set in.MultiModel;
proc sort; by newid wavenum; run; /* 232218 observations and 18 variables.*/


/********************************************************************************************************************************************************************************/
/*** Create bootstrap samples using PROC SURVEYSELECT ***/

*Create sample size input data set to provide the stratum sample sizes in the bootstrapping; 
proc freq data=input noprint;
  tables RAESTRAT*RAEHSAMP/out=nsize(rename=(count=_nsize_));
run;
proc print ; run;

*Run SURVEYSELECT to generate data with no replicated Ids, I will use variable NumberHits to replicate these Ids later;
proc surveyselect data=input out=bsample method=urs sampsize=nsize reps=100 seed=1953;
strata RAESTRAT RAEHSAMP;
run;
/*Total Sample size is : 32,486*100=32,486,000  */


/********************************************************************************************************************************************************************************/
/*** Merge each of the bootstrap samples with data that has repeated observations/ Id ***/

%macro bsmulti (B=); 

 %do i=1 %to &B; /*for each bootstrap sample B*/

	*Each Id  needs to be repeated as many times as their normwgt_198int*NumberHits;
	data bsample2;
	 set bsample;
	 where replicate=&i;
	 countId=normwgt_198int*NumberHits;
	proc sort; by newid; run;

	*Merge with replicated data;
	data bsample3;
	 merge multimodel bsample2(in=A keep=newid countId NumberHits SamplingWeight replicate);
	 by newid;
	 if A;
	proc sort; by newid wavenum; run;

	* Repeat each newid and iwdate up to their countid;
	data bsample4(drop=i count check);
	  set bsample3;
	  if _N_=1 then count=1; /* start count */
	  by newid wavenum; /* for each newid and iwdate */
	  do i=1 to countId; /* replicate each newid and iwdate up to their countId */
	   newid2=compbl(catx(".", newid, count)); /* create newid2 variable */
	   /*catx: concatenate (join) two or more character strings, stripping both leading and trailing blanks and inserting "." separator between the strings*/
	   count+1; /* for each replicate of newid and iwdate count will increase 1 */
	   check=count-countId;
	   if check=1 then count=1; /*check=1 the first time count>coundtId, so that means that we reach the total number of replicates needed for Id and its iw, so we need to reset count to 1*/
	   output;
	  end;
	run;

   * Create permanent dataset and csv file so we can use in R to fit Multistate model;
	data out.bs&i; 
	 set bsample4;
	 label newid2="Fake respondent Id for bootstrap samples"
	       countId="Number of times R is selected=normwgt_198int*NumberHits";
	run;

   proc delete data=bsample2 bsample3 bsample4; run; quit;

  %end;	
%mend bsmulti;
%bsmulti(B=100);





