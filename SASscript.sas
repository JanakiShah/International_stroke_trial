/*********************************************
Step 1: Load the IST dataset
*********************************************/
proc import 
  datafile="/home/u62275281/ILE/DataIST.csv" out=dIST   dbms=dlm    replace;
delimiter=',';
 getnames=yes;
run;

/*********************************************
Step 2: Data Cleaning and Preparation
*********************************************/
/* Check for missing values and handle as appropriate */

/* Create binary variable for stroke type (1 = ischemic, 2 = hemorrhagic) */
proc format;
value stroke_typef 1='ischemic' 2='hemorrhagic';
value FDEADCf 2='ischemic' 3='hemorrhagic';
value  DDEADCf 2='ischemic' 3='hemorrhagic';
run;


data ist_new;
  set dIST;
  if DDIAGISC='Y' then stroke_type= '1';
 else  if DDIAGHA='Y' then  stroke_type= '2';
  else if stroke_type= '';
  run;
  
  data istclean;
  Set ist_new;
  if stroke_type='' then delete;
  if DDEAD ='Y' then deadforteen=1;	
else if DDEAD ='N' then deadforteen=0;	
if FDEAD ='Y' then dead=1;	
else if FDEAD ='N' then dead=0;	
format  stroke_type  stroke_typef. FDEADCf FDEADCf.  DDEADC DDEADCf.;
 label DDEADD='Follow-up time (days)' FDEADD='Follow-up time (days)';
  run;
   
/*********************************************
Step 3: Descriptive Statistics
*********************************************/
/* Generate descriptive statistics for patient characteristics */
proc means data=istclean maxdec=2;
class stroke_type;
  var age RSBP ;
  run;
  proc sort data=istclean;
  by stroke_type;
  run;
  proc freq data=istclean;
  by stroke_type;
  tables sex RSLEEP RCT RATRIAL RVISINF  RHEP24 RCONSC RASP3 RDEF1 RDEF2 RDEF3 RDEF4
  RDEF5 RDEF6 RDEF7 RDEF8 STYPE DRSISC DRSH  RXASP RXHEP/ oddsratio cl;
run;



proc freq data=istclean;
tables DRSISC DRSH;
run;

proc freq data=istclean;
tables DRSISC*FDEADC;
run;
proc freq data=istclean;
tables DRSH*FDEADC;
run;


/*********************************************
Step 4: Survival Analysis - Kaplan-Meier Curves
*********************************************/


*kaplan mier curve for comparision of Ischemic and haemorrhagic Strokes at 6months ;

proc lifetest data=istclean plots=survival(atrisk) notable;
   where FDEADC in (2, 3);
   time FDEADD*DEAD(0);
   strata FDEADC;
   format FDEADC FDEADCf.;
   title1 'Kaplan-Meier Survival at 6months';
   title2 'by Ischemic and haemorrhagic Strokes';
run;

proc lifetest data=istclean plots=survival(atrisk) notable;
   where DDEADC in (2, 3);
   time DDEADD*deadforteen(0);
   strata DDEADC;
   format DDEADC DDEADCf.;
   title1 'Kaplan-Meier Survival at 14 days';
   title2 'by Ischemic and haemorrhagic Strokes';
run;



proc lifetest data=istclean ;
   where FDEADC in (2, 3);
   time FDEADD*DEAD(0);
   strata FDEADC;
   format FDEADC FDEADCf.;
   title1 'Kaplan-Meier Survival at 6months';
   title2 'by Ischemic and haemorrhagic Strokes';
run;
proc lifetest data=istclean ;
   where DDEADC in (2, 3);
   time DDEADD*deadforteen(0);
   strata DDEADC;
   format DDEADC DDEADCf.;
   title1 'Kaplan-Meier Survival at 14 days';
   title2 'by Ischemic and haemorrhagic Strokes';
run;
/*********************************************
Step 5: Cox Proportional Hazards Regression
*********************************************/


  proc phreg data=istclean;
   where FDEADC in (2, 3);
class RXASP(ref='Y')  RSLEEP(ref='N') RCT(ref='Y') RATRIAL(ref='N') RVISINF(ref='N') RCONSC(ref='F') SEX(ref='F') stroke_type(ref='2');
  model FDEADD*DEAD(0) = stroke_type  RXASP  RATRIAL RVISINF RCONSC RSLEEP SEX AGE/ ties=exact risklimits= wald;
strata FDEADC;
hazardratio stroke_type;
run;


proc phreg data=istclean;
 where DDEADC in (2, 3);
class RXASP RXHEP RSLEEP(ref='N') RCT(ref='N') RATRIAL(ref='N') RVISINF(ref='N') RCONSC(ref='F') SEX(ref='M') stroke_type(ref='1');
 model DDEADD*deadforteen(0) =  RXASP  RATRIAL RVISINF RCONSC RSLEEP SEX AGE/ ties=exact risklimits= wald;
strata DDEADC;
run;


proc sort data=istclean;
by stroke_type;
run;

/*********************************************
Step 6: Additional Subgroup Analyses
*********************************************/

data ist_new2;
  set dIST;
   if DDIAGISC='Y' then stroke_type= '1';
  else if stroke_type= '';
  run;
  
  data istclean2;
  Set ist_new;
  if stroke_type='' then delete;
  if DDEAD ='Y' then deadforteen=1;	
else if DDEAD ='N' then deadforteen=0;	
if FDEAD ='Y' then dead=1;	
else if FDEAD ='N' then dead=0;	
  run;
  

  proc phreg data=istclean2;
class RXASP RXHEP RSLEEP(ref='N') RCT(ref='N') RATRIAL(ref='N') RVISINF(ref='N') RCONSC(ref='F') SEX(ref='M') stroke_type(ref='1');
  model FDEADD*DEAD(0) =  RXASP  RATRIAL RVISINF RCONSC RSLEEP SEX AGE/ ties=exact risklimits= wald;
run;

proc phreg data=istclean2;
class RXASP RXHEP RSLEEP(ref='N') RCT(ref='N') RATRIAL(ref='N') RVISINF(ref='N') RCONSC(ref='F') SEX(ref='M') stroke_type(ref='1');
  model DDEADD*deadforteen(0)  =  RXASP  RATRIAL RVISINF RCONSC RSLEEP SEX AGE/ ties=exact risklimits= wald;
run;

/* Perform Cox proportional hazard analysis */
proc phreg data=istclean2;
  class RSLEEP(ref='N') RCT(ref='N') RATRIAL(ref='N') RVISINF(ref='N') RCONSC(ref='F') SEX(ref='M') stroke_type(ref='1');
  model FDEADD*DEAD(0) = RATRIAL RVISINF RCONSC RSLEEP SEX AGE / ties=exact risklimits=wald ;
  hazardratio  RATRIAL RVISINF RCONSC RSLEEP SEX AGE / CL=WALD;
  
run;


