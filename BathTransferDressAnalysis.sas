/*********** Investigator: Rebecca Brown ***********/
/*********** Statistician: Grisell Diaz-Ramirez *********** /
/*********** Date created: 2021.06.14 ***********/
/*********** Purpose: Simultaneous analysis of bathing, transferring, and dressing difficulty ***********/


libname rap 'path';
libname rand 'path';
proc format cntlin=rand.sasfmts; run;


*Get rap.cohort_adl_iadl2 dataset that has all ADLs/IADLs;
proc contents data=rap.cohort_adl_iadl2; run;
proc sort data=rap.cohort_adl_iadl2 out=all; by newid; run;
/*32486 observations and 369 variables.*/

proc freq data=all; tables bathdifbl bathpredif dressdifbl dresspredif beddifbl bedpredif; run;
proc print data=all; 
 var newid bathdifbl bathpredif dressdifbl dresspredif beddifbl bedpredif;
 where bathdifbl=. or dressdifbl=. or beddifbl=.;
run; /*6 Ids with missing bathdifbl/beddifbl*/

*Create new variables to record ADLdiff baseline or before;
data all2;
 set all;
 bathdifprebl=max(bathdifbl, bathpredif);
 dressdifprebl=max(dressdifbl, dresspredif);
 beddifprebl=max(beddifbl, bedpredif);

 if bathdifbl=. or dressdifbl=. or beddifbl=. then delete;

 label bathdifprebl='Bathing diff at baseline or before. 0.no, 1.yes'
       dressdifprebl='Dressing diff at baseline or before. 0.no, 1.yes' 
       beddifprebl='Transferring diff at baseline or before. 0.no, 1.yes' ;
run;
/* 32486-6=32480 observations and 369+2=372 variables. */


/* QC */
proc freq data=all2; tables bathdifprebl*bathdifbl*bathpredif dressdifprebl*dressdifbl*dresspredif beddifprebl*beddifbl*bedpredif /list missing; run;
proc freq data=all2; tables bathdifprebl*dressdifprebl*beddifprebl /list missing; run;
proc freq data=all2; tables agegroupbl; run;

/* Difficulty at baseline or before */

*Gral;
proc freq data=all2; tables bathdifprebl*dressdifprebl*beddifprebl /list missing; run;
proc surveyfreq data=all2;
 weight weightbl;
 stratum raestrat; 
 cluster raehsamp; 
 tables bathdifprebl*dressdifprebl*beddifprebl /  cl ;
run;

*By agegroup at baseline;
proc freq data=all2; tables bathdifprebl*dressdifprebl*beddifprebl /list; where agegroupbl=1; run; *agegroupbl:50-64;
proc freq data=all2; tables bathdifprebl*dressdifprebl*beddifprebl /list; where agegroupbl=2; run; *agegroupbl:65-74;
proc freq data=all2; tables bathdifprebl*dressdifprebl*beddifprebl /list; where agegroupbl=3; run; *agegroupbl:75-84;
proc freq data=all2; tables bathdifprebl*dressdifprebl*beddifprebl /list; where agegroupbl=4; run; *agegroupbl:85+;

proc surveyfreq data=all2;
 weight weightbl;
 stratum raestrat; 
 cluster raehsamp; 
 tables agegroupbl*bathdifprebl*dressdifprebl*beddifprebl / cl ;
run;

/* Difficulty at follow-up for those with bathdifprebl=0 and dressdifprebl=0 and beddifprebl=0 */
data all3;
 set all2;
 where bathdifprebl=0 and dressdifprebl=0 and beddifprebl=0;
run;
/* 28659 observations and 372 variables. */

proc freq data=all3; tables bathdifprebl dressdifprebl beddifprebl; run;

/********************************************************** Macro to compute difficulty for each ADL *****************************************************************************/

%macro cohort (outputdata=, diff=, randvars=);

data &outputdata (drop=wave1-wave13);
 set all3;

 &diff=0;

 array dif[13]     &randvars; 
 array agecont[13]  age1 	age2     age3    	age4     age5     age6     age7     age8     age9     age10     age11		age12		age13;
 array ageint[13]	ageint1 ageint2  ageint3 	ageint4  ageint5  ageint6  ageint7  ageint8  ageint9  ageint10  ageint11	ageint12	ageint13;
 array IWDT[13] 	iwdate1 iwdate2  iwdate3 	iwdate4  iwdate5  iwdate6  iwdate7	iwdate8  iwdate9  iwdate10  iwdate11	iwdate12	iwdate13;
 array ITYPE[13] 	R1IWSTAT R2IWSTAT R3IWSTAT  R4IWSTAT R5IWSTAT R6IWSTAT R7IWSTAT R8IWSTAT R9IWSTAT R10IWSTAT R11IWSTAT   R12IWSTAT	R13IWSTAT;
 array wave[13] 	(1,2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,13); 

 do i=2 to 13; /*start in wave 2 since in wave 1 all Ids with wave1stelig=1 should be adl/iadl free*/
  if wave[i]>=wave1stelig and &diff=0 and dif[i]=1 then do;
    &diff=1; /*make the indicator variable=1 for the first time. This is the first time a participant has diff*/
    &diff.wave=wave[i];

    if IWDT[i-1]>0 and IWDT[i]>0 and ITYPE[i-1]=1 and ITYPE[i]=1 then do; 
	 time2&diff=yrdif(iwdate, IWDT[i-1],'actual')+((yrdif(IWDT[i-1],IWDT[i],'actual'))/2);  
	 /*Time to event: (get the time between entering the study and free-disease time)plus(get the time between when the event happened and the previous time/2: interval censoring )*/
	 age&diff=agecont[i-1]+((yrdif(IWDT[i-1],IWDT[i],'ACTUAL'))/2);
    end;
	else if IWDT[i-1]=. or ITYPE[i-1] ne 1 then do; 
	 time2&diff=yrdif(iwdate, IWDT[i],'actual');
	 age&diff=agecont[i];
	end;

  end;
 end; drop i;

 label &diff="R had any &diff diff during follow-up time. 0.no, 1.yes"
	   time2&diff="R time to 1st &diff diff or follow-up time in yrs"
	   age&diff="R age at 1st &diff diff or at last follow-up interview"
       &diff.wave="R 1st wave of &diff diff";

  /*Get last date when interviewed with diff=0*/
 &diff.enddt=.;
 do i=13 to 1 by -1;
	if IWDT[i]^=. and ITYPE[i]=1 and dif[i]=0 and &diff.enddt=. then do; &diff.enddt=IWDT[i]; &diff.difendage=agecont[i]; end;
 end; drop i;
 format &diff.enddt mmddyy10.;
 label &diff.enddt="R last date with core interview and &diff.dif=0"
       &diff.difendage="R age at last date with core interview and &diff.dif=0";

 /*Competing-risk outcome*/
 status_&diff=&diff; time_&diff=time2&diff; age_&diff=age&diff;

 if &diff=0 and death=1 then do; status_&diff=2; time_&diff=timetodeath; age_&diff=agedeathcont; end;
 else if &diff=0 and death=0 then do; time_&diff=yrdif(iwdate, &diff.enddt,'actual'); age_&diff=&diff.difendage; end;

 label status_&diff="R had any &diff diff during follow-up time. 0.no, 1.yes, 2.died"
	   time_&diff="R time to 1st &diff diff, death, or follow-up time in yrs"
	   age_&diff="R age at 1st &diff diff, death, or at last follow-up interview";

 /*Define age-group cohort base on diff at follow-up*/
 if 50<=age_&diff<65 then age_&diff._group=1;
 else if 65<=age_&diff<75 then age_&diff._group=2;
 else if 75<=age_&diff<85 then age_&diff._group=3;
 else if age_&diff>=85 then age_&diff._group=4;
 label age_&diff._group="R age at 1st dif. 1.50-64, 2.65-74, 3.75-84, 4.85+";

/* if time_difficulty=0 then delete;*/

run;

%mend cohort;

/* Bath */
%cohort(outputdata=bath, diff=bath, randvars=R1BATHW R2BATHA  R3BATHA R4BATHA  R5BATHA R6BATHA R7BATHA R8BATHA  R9BATHA R10BATHA R11BATHA R12BATHA R13BATHA);
/* 28659 observations and 372+10-13=369 variables. */

/* QC bath */
proc contents data=bath position; run;
proc freq data=bath; tables bath*death*status_bath /missing list; run;
proc freq data=bath; tables bath*bathwave /missing list; run;
proc freq data=bath; tables status_bath*bathwave /missing list; run;
proc freq data=bath; tables age_bath_group; run;
proc means data=bath n nmiss min max; var age_bath; class age_bath_group; run;
proc freq data=bath; tables agegroupbl*age_bath_group /list ; run;
proc means data=bath n nmiss min max; var time2bath agebath time_bath age_bath bathenddt bathdifendage timetodeath; run;


*Check some Ids;
proc print data=bath (obs=2); var newid RACOHBYR wave1stelig iwdate iwdate1-iwdate13 ;
where agegroupbl=2 and  status_bath=1; run;
proc print data=bath (obs=2); var newid R1IWSTAT R2IWSTAT R3IWSTAT  R4IWSTAT R5IWSTAT R6IWSTAT R7IWSTAT R8IWSTAT R9IWSTAT R10IWSTAT R11IWSTAT R12IWSTAT	R13IWSTAT;
where agegroupbl=2 and  status_bath=1; run;
proc print data=bath (obs=2); var newid dod agecontbl agegroupbl status_bath age_bath age_bath_group time_bath bathdifendage bathenddt date65 date75 date85;
where agegroupbl=2 and status_bath=1; run;
proc print data=bath (obs=2); var newid bathwave R1BATHW R2BATHA  R3BATHA R4BATHA  R5BATHA R6BATHA R7BATHA R8BATHA  R9BATHA R10BATHA R11BATHA R12BATHA R13BATHA;
where agegroupbl=2 and status_bath=1; run;
proc print data=bath (obs=2); var newid age1 	age2     age3    	age4     age5     age6     age7     age8     age9     age10     age11		age12		age13;
where agegroupbl=2 and status_bath=1; run;


/* Dress */
%cohort(outputdata=dress, diff=dress, randvars=R1DRESSW R2DRESSA R3DRESSA R4DRESSA R5DRESSA R6DRESSA R7DRESSA  R8DRESSA R9DRESSA R10DRESSA R11DRESSA R12DRESSA R13DRESSA);
/* 28659 observations and 372+10-13=369 variables. */

/* QC dress */
proc contents data=dress position; run;
proc freq data=dress; tables dress*death*status_dress /missing list; run;
proc freq data=dress; tables dress*dresswave /missing list; run;
proc freq data=dress; tables status_dress*dresswave /missing list; run;
proc freq data=dress; tables age_dress_group; run;
proc means data=dress n nmiss min max; var age_dress; class age_dress_group; run;
proc freq data=dress; tables agegroupbl*age_dress_group /list ; run;
proc means data=dress n nmiss min max; var time2dress agedress time_dress age_dress dressenddt dressdifendage timetodeath; run;


*Check some Ids;
proc print data=dress (obs=2); var newid RACOHBYR wave1stelig iwdate iwdate1-iwdate13 ;
where agegroupbl=3 and  status_dress=1; run;
proc print data=dress (obs=2); var newid R1IWSTAT R2IWSTAT R3IWSTAT  R4IWSTAT R5IWSTAT R6IWSTAT R7IWSTAT R8IWSTAT R9IWSTAT R10IWSTAT R11IWSTAT R12IWSTAT	R13IWSTAT;
where agegroupbl=3 and  status_dress=1; run;
proc print data=dress (obs=2); var newid dod agecontbl agegroupbl status_dress age_dress age_dress_group time_dress dressdifendage dressenddt date65 date75 date85;
where agegroupbl=3 and status_dress=1; run;
proc print data=dress (obs=2); var newid dresswave R1DRESSW R2DRESSA R3DRESSA R4DRESSA R5DRESSA R6DRESSA R7DRESSA  R8DRESSA R9DRESSA R10DRESSA R11DRESSA R12DRESSA R13DRESSA;
where agegroupbl=3 and status_dress=1; run;
proc print data=dress (obs=2); var newid age1 	age2     age3    	age4     age5     age6     age7     age8     age9     age10     age11		age12		age13;
where agegroupbl=3 and status_dress=1; run;

/* bed */
%cohort(outputdata=bed, diff=bed, randvars=R1BEDW R2BEDA R3BEDA R4BEDA R5BEDA R6BEDA R7BEDA R8BEDA R9BEDA R10BEDA R11BEDA R12BEDA R13BEDA);
/* 28659 observations and 372+10-13=369 variables. */

/* QC bed */
proc contents data=bed position; run;
proc freq data=bed; tables bed*death*status_bed /missing list; run;
proc freq data=bed; tables bed*bedwave /missing list; run;
proc freq data=bed; tables status_bed*bedwave /missing list; run;
proc freq data=bed; tables age_bed_group; run;
proc means data=bed n nmiss min max; var age_bed; class age_bed_group; run;
proc freq data=bed; tables agegroupbl*age_bed_group /list ; run;
proc means data=bed n nmiss min max; var time2bed agebed time_bed age_bed bedenddt beddifendage timetodeath; run;


*Check some Ids;
proc print data=bed (obs=2); var newid RACOHBYR wave1stelig iwdate iwdate1-iwdate13 ;
where agegroupbl=4 and  status_bed=1; run;
proc print data=bed (obs=2); var newid R1IWSTAT R2IWSTAT R3IWSTAT  R4IWSTAT R5IWSTAT R6IWSTAT R7IWSTAT R8IWSTAT R9IWSTAT R10IWSTAT R11IWSTAT R12IWSTAT	R13IWSTAT;
where agegroupbl=4 and  status_bed=1; run;
proc print data=bed (obs=2); var newid dod agecontbl agegroupbl status_bed age_bed age_bed_group time_bed beddifendage bedenddt date65 date75 date85;
where agegroupbl=4 and status_bed=1; run;
proc print data=bed (obs=2); var newid bedwave R1BEDW R2BEDA R3BEDA R4BEDA R5BEDA R6BEDA R7BEDA R8BEDA R9BEDA R10BEDA R11BEDA R12BEDA R13BEDA;
where agegroupbl=4 and status_bed=1; run;
proc print data=bed (obs=2); var newid age1 	age2     age3    	age4     age5     age6     age7     age8     age9     age10     age11		age12		age13;
where agegroupbl=4 and status_bed=1; run;

proc sort data=bath; by newid; run;
proc sort data=dress; by newid; run;
proc sort data=bed; by newid; run;

*Merge bath, dress, and bed datasets;
data bath_dress_bed;
 merge bath dress(keep=newid dress dresswave time2dress agedress dressenddt dressdifendage status_dress time_dress age_dress age_dress_group)
            bed  (keep=newid bed bedwave time2bed agebed bedenddt beddifendage status_bed time_bed age_bed age_bed_group);
 by newid;
 length order $30;

 adlsum_baddressbed=sum(bath,dress,bed);

 if adlsum_baddressbed=0 then order="No_3_ADLs";

 else if adlsum_baddressbed=1 and bath=1 and dress=0 and bed=0 then order="bath";
 else if adlsum_baddressbed=1 and bath=0 and dress=1 and bed=0 then order="dress";
 else if adlsum_baddressbed=1 and bath=0 and dress=0 and bed=1 then order="bed";

 else if adlsum_baddressbed=2 and bath=1 and dress=1 and bed=0 then do;
  if bathwave<dresswave then order="bath1st_dress2nd";
  else if bathwave>dresswave then order="dress1st_bath2nd";
  else if bathwave=dresswave then order="bath_dress";
 end;
 else if adlsum_baddressbed=2 and bath=1 and dress=0 and bed=1 then do;
  if bathwave<bedwave then order="bath1st_bed2nd";
  else if bathwave>bedwave then order="bed1st_bath2nd";
  else if bathwave=bedwave then order="bath_bed";
 end;
 else if adlsum_baddressbed=2 and bath=0 and dress=1 and bed=1 then do;
  if dresswave<bedwave then order="dress1st_bed2nd";
  else if dresswave>bedwave then order="bed1st_dress2nd";
  else if dresswave=bedwave then order="dress_bed";
 end;

 else if adlsum_baddressbed=3 then do;
  if bathwave=dresswave and dresswave=bedwave then order="bath_dress_bed";

  else if bathwave=dresswave and dresswave<bedwave then order="bathdress1st_bed2nd";
  else if bathwave=dresswave and dresswave>bedwave then order="bed1st_bathdress2nd";

  else if bathwave=bedwave and bedwave<dresswave then order="bathbed1st_dress2nd";
  else if bathwave=bedwave and bedwave>dresswave then order="dress1st_bathbed2nd";

  else if dresswave=bedwave and bedwave<bathwave then order="dressbed1st_bath2nd";
  else if dresswave=bedwave and bedwave>bathwave then order="bath1st_dressbed2nd";

  else if bathwave<dresswave and dresswave<bedwave then order="bath1st_dress2nd_bed3rd";
  else if bathwave<bedwave and bedwave<dresswave then order="bath1st_bed2nd_dress3rd";

  else if dresswave<bathwave and bathwave<bedwave then order="dress1st_bath2nd_bed3rd";
  else if dresswave<bedwave and bedwave<bathwave then order="dress1st_bed2nd_bath3rd";

  else if bedwave<bathwave and bathwave<dresswave then order="bed1st_bath2nd_dress3rd";
  else if bedwave<dresswave and dresswave<bathwave then order="bed1st_dress2nd_bath3rd";
 end;

 label adlsum_baddressbed="R total number of ADLs during follow-up (0-3)"
       order="Order of ADLs among bath, dress, bed";
 
run;
/* 28659 observations and 369+20+2=391 variables. */

/* QC */
proc freq data=bath_dress_bed; tables bath*dress*bed /list missing; run;
proc freq data=bath_dress_bed; tables bathwave / missing; where bath=1; run;
proc freq data=bath_dress_bed; tables dresswave / missing; where dress=1; run;
proc freq data=bath_dress_bed; tables bedwave / missing; where bed=1; run;

*1 ADL;
proc freq data=bath_dress_bed; tables order*bathwave*dresswave*bedwave /list missing; where adlsum_baddressbed=1; run;

*2 ADLs;
proc freq data=bath_dress_bed; tables order*bathwave*dresswave /list missing; where adlsum_baddressbed=2 and bath=1 and dress=1 and bed=0; run;
proc freq data=bath_dress_bed; tables order*bathwave*bedwave /list missing; where adlsum_baddressbed=2 and bath=1 and dress=0 and bed=1; run;
proc freq data=bath_dress_bed; tables order*dresswave*bedwave /list missing; where adlsum_baddressbed=2 and bath=0 and dress=1 and bed=1; run;
proc freq data=bath_dress_bed; tables order; run;
proc freq data=bath_dress_bed; tables adlsum_baddressbed*order /list missing; run;

*3 ADLs;
proc freq data=bath_dress_bed; tables order*bathwave*dresswave*bedwave /list missing; where adlsum_baddressbed=3; run;


*Check some Ids with adlsum_baddressbed=0;
proc print data=bath_dress_bed (obs=2); var newid RACOHBYR wave1stelig iwdate iwdate1-iwdate13 ;
where adlsum_baddressbed=0; run;
proc print data=bath_dress_bed (obs=2); var newid R1IWSTAT R2IWSTAT R3IWSTAT  R4IWSTAT R5IWSTAT R6IWSTAT R7IWSTAT R8IWSTAT R9IWSTAT R10IWSTAT R11IWSTAT R12IWSTAT R13IWSTAT;
where adlsum_baddressbed=0; run;
proc print data=bath_dress_bed (obs=2); var newid bathwave R1BATHW R2BATHA  R3BATHA R4BATHA  R5BATHA R6BATHA R7BATHA R8BATHA  R9BATHA R10BATHA R11BATHA R12BATHA R13BATHA;
where adlsum_baddressbed=0; run;
proc print data=bath_dress_bed (obs=2); var newid dresswave R1DRESSW R2DRESSA R3DRESSA R4DRESSA R5DRESSA R6DRESSA R7DRESSA  R8DRESSA R9DRESSA R10DRESSA R11DRESSA R12DRESSA R13DRESSA;
where adlsum_baddressbed=0; run;
proc print data=bath_dress_bed (obs=2); var newid bedwave R1BEDW R2BEDA R3BEDA R4BEDA R5BEDA R6BEDA R7BEDA R8BEDA R9BEDA R10BEDA R11BEDA R12BEDA R13BEDA;
where adlsum_baddressbed=0; run;


*Check some Ids with adlsum_baddressbed=1;
proc print data=bath_dress_bed (obs=2); var newid RACOHBYR wave1stelig order iwdate iwdate1-iwdate13 ;
where adlsum_baddressbed=1; run;
proc print data=bath_dress_bed (obs=2); var newid R1IWSTAT R2IWSTAT R3IWSTAT  R4IWSTAT R5IWSTAT R6IWSTAT R7IWSTAT R8IWSTAT R9IWSTAT R10IWSTAT R11IWSTAT R12IWSTAT R13IWSTAT;
where adlsum_baddressbed=1; run;
proc print data=bath_dress_bed (obs=2); var newid bathwave R1BATHW R2BATHA  R3BATHA R4BATHA  R5BATHA R6BATHA R7BATHA R8BATHA  R9BATHA R10BATHA R11BATHA R12BATHA R13BATHA;
where adlsum_baddressbed=1; run;
proc print data=bath_dress_bed (obs=2); var newid dresswave R1DRESSW R2DRESSA R3DRESSA R4DRESSA R5DRESSA R6DRESSA R7DRESSA  R8DRESSA R9DRESSA R10DRESSA R11DRESSA R12DRESSA R13DRESSA;
where adlsum_baddressbed=1; run;
proc print data=bath_dress_bed (obs=2); var newid bedwave R1BEDW R2BEDA R3BEDA R4BEDA R5BEDA R6BEDA R7BEDA R8BEDA R9BEDA R10BEDA R11BEDA R12BEDA R13BEDA;
where adlsum_baddressbed=1; run;

*Check some Ids with adlsum_baddressbed=2;
proc print data=bath_dress_bed (obs=5); var newid RACOHBYR wave1stelig order iwdate iwdate1-iwdate13 ;
where adlsum_baddressbed=2; run;
proc print data=bath_dress_bed (obs=5); var newid R1IWSTAT R2IWSTAT R3IWSTAT  R4IWSTAT R5IWSTAT R6IWSTAT R7IWSTAT R8IWSTAT R9IWSTAT R10IWSTAT R11IWSTAT R12IWSTAT R13IWSTAT;
where adlsum_baddressbed=2; run;
proc print data=bath_dress_bed (obs=5); var newid bathwave R1BATHW R2BATHA  R3BATHA R4BATHA  R5BATHA R6BATHA R7BATHA R8BATHA  R9BATHA R10BATHA R11BATHA R12BATHA R13BATHA;
where adlsum_baddressbed=2; run;
proc print data=bath_dress_bed (obs=5); var newid dresswave R1DRESSW R2DRESSA R3DRESSA R4DRESSA R5DRESSA R6DRESSA R7DRESSA  R8DRESSA R9DRESSA R10DRESSA R11DRESSA R12DRESSA R13DRESSA;
where adlsum_baddressbed=2; run;
proc print data=bath_dress_bed (obs=5); var newid bedwave R1BEDW R2BEDA R3BEDA R4BEDA R5BEDA R6BEDA R7BEDA R8BEDA R9BEDA R10BEDA R11BEDA R12BEDA R13BEDA;
where adlsum_baddressbed=2; run;

*Check some Ids with adlsum_baddressbed=3;
proc print data=bath_dress_bed (obs=5); var newid RACOHBYR wave1stelig order iwdate iwdate1-iwdate13 ;
where adlsum_baddressbed=3; run;
proc print data=bath_dress_bed (obs=5); var newid R1IWSTAT R2IWSTAT R3IWSTAT  R4IWSTAT R5IWSTAT R6IWSTAT R7IWSTAT R8IWSTAT R9IWSTAT R10IWSTAT R11IWSTAT R12IWSTAT R13IWSTAT;
where adlsum_baddressbed=3; run;
proc print data=bath_dress_bed (obs=5); var newid bathwave R1BATHW R2BATHA  R3BATHA R4BATHA  R5BATHA R6BATHA R7BATHA R8BATHA  R9BATHA R10BATHA R11BATHA R12BATHA R13BATHA;
where adlsum_baddressbed=3; run;
proc print data=bath_dress_bed (obs=5); var newid dresswave R1DRESSW R2DRESSA R3DRESSA R4DRESSA R5DRESSA R6DRESSA R7DRESSA  R8DRESSA R9DRESSA R10DRESSA R11DRESSA R12DRESSA R13DRESSA;
where adlsum_baddressbed=3; run;
proc print data=bath_dress_bed (obs=5); var newid bedwave R1BEDW R2BEDA R3BEDA R4BEDA R5BEDA R6BEDA R7BEDA R8BEDA R9BEDA R10BEDA R11BEDA R12BEDA R13BEDA;
where adlsum_baddressbed=3; run;


/* Stats Gral */
proc freq data=bath_dress_bed; tables adlsum_baddressbed*bath*dress*bed /list missing; run;
proc surveyfreq data=bath_dress_bed;
 weight weightbl;
 stratum raestrat; 
 cluster raehsamp; 
 tables bath*dress*bed / cl;
run;

*Order;
proc freq data=bath_dress_bed; tables order; run;
proc freq data=bath_dress_bed; tables adlsum_baddressbed*order /list missing; run;

proc surveyfreq data=bath_dress_bed;
 weight weightbl;
 stratum raestrat; 
 cluster raehsamp; 
 tables order / cl;
run;

*Order by agegroupbl;
proc freq data=bath_dress_bed; tables agegroupbl; run;
proc freq data=bath_dress_bed; tables adlsum_baddressbed*order / list missing; where agegroupbl=1; run;
proc freq data=bath_dress_bed; tables adlsum_baddressbed*order / list missing; where agegroupbl=2; run;
proc freq data=bath_dress_bed; tables adlsum_baddressbed*order / list missing; where agegroupbl=3; run;
proc freq data=bath_dress_bed; tables adlsum_baddressbed*order / list missing; where agegroupbl=4; run;

proc surveyfreq data=bath_dress_bed;
 weight weightbl;
 stratum raestrat; 
 cluster raehsamp; 
 tables agegroupbl*order / domain=row cl ; /* for domain analysis see: https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.3/statug/statug_surveyfreq_details06.htm#:~:text=You%20can%20perform%20domain%20analysis,for%20the%20entire%20study%20population.*/
run;


*Create permanent dataset: bath_dress_bed;
data rap.bath_dress_bed;
 set bath_dress_bed;
proc sort; by newid; run;
/*28659 observations and 391 variables*/

proc contents data=rap.bath_dress_bed position; run;

