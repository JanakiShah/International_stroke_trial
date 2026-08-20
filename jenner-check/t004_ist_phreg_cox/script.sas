/*********************************************
Step 1: Load the IST dataset
(adapted from the repo's SASscript.sas PROC IMPORT step: original
reads /home/u62275281/ILE/DataIST.csv; here we inline a small
representative sample of IST-shaped rows via DATALINES so the
bundle is self-contained)
*********************************************/
data dIST;
  input ID AGE SEX$ RSBP RCONSC$ RSLEEP$ RCT$ RATRIAL$ RVISINF$ RHEP24$ RASP3$ RDEF1$ RDEF2$ RDEF3$ RDEF4$ RDEF5$ RDEF6$ RDEF7$ RDEF8$ STYPE$ RXASP$ RXHEP$ DDIAGISC$ DDIAGHA$ DDEAD$ DDEADD DDEADC FDEAD$ FDEADD FDEADC DRSISC DRSH;
  datalines;
1 70 M 183 F Y Y N Y Y Y Y N N Y Y Y N Y LACS Y Y Y N N 14 2 Y 57 2 0 0
2 63 M 153 F Y N Y Y Y N Y Y Y Y N N N N POCS N M N Y N 14 3 Y 179 3 0 1
3 78 F 163 D N N Y Y N Y N Y N N Y Y N N TACS N L Y N Y 14 2 Y 70 2 3 0
4 48 M 193 U N N N N N Y N N Y Y N Y Y N PACS Y L N Y N 14 3 N 180 3 0 0
5 80 F 135 F N N N N N Y Y Y Y Y Y Y Y N PACS N M Y N Y 3 2 N 180 2 4 0
6 53 F 188 U Y N N N N N Y N N Y Y Y Y N PACS Y M Y N Y 2 2 Y 146 2 1 0
7 84 F 103 F Y N Y N N N N Y Y N N N N N LACS Y Y Y N N 14 2 N 180 2 2 0
8 78 M 102 F N Y Y N Y N N Y N Y N Y Y Y POCS Y N Y N N 14 2 N 180 2 5 0
9 62 M 160 D Y N N N N Y Y Y Y N Y N Y N LACS N M Y N Y 14 2 N 180 2 0 0
10 75 M 122 D N Y N N N Y Y Y Y Y Y N Y N TACS Y N N Y Y 1 3 N 180 3 0 5
11 72 M 124 F Y N Y N Y N N N Y Y N N N Y PACS Y L Y N N 14 2 N 180 2 0 0
12 56 M 118 D Y Y N N Y Y Y Y N Y Y N Y Y POCS N N N Y N 14 3 N 180 3 0 4
13 77 F 131 U N Y N Y N Y N N N Y Y N Y Y TACS Y N Y N N 14 2 Y 65 2 1 0
14 51 M 150 D Y Y Y N N N N Y N N Y N Y N POCS N Y N Y N 14 3 N 180 3 0 4
15 52 M 200 F Y Y N N Y Y N Y N N N Y N N LACS N Y N Y N 14 3 N 180 3 0 0
16 85 M 111 D Y Y Y N Y N Y N N N Y Y Y Y PACS N Y Y N N 14 2 Y 80 2 5 0
17 63 M 157 U Y N N Y N Y Y Y Y N Y N Y N POCS N M Y N N 14 2 Y 88 2 1 0
18 70 M 144 F Y Y Y N N Y Y Y N N Y N Y N PACS Y M N Y N 14 3 Y 68 3 0 2
19 60 F 104 D Y N Y Y N N Y N N Y Y Y Y N LACS Y L N Y Y 7 3 Y 77 3 0 2
20 82 M 167 F N N N Y N Y Y N Y Y Y Y Y Y PACS N Y N Y N 14 3 N 180 3 0 4
21 85 M 168 U Y N N Y N Y Y Y N N Y N Y Y PACS N L Y N N 14 2 Y 123 2 5 0
22 84 M 180 U Y Y Y N N N Y Y N Y N N Y Y POCS N M Y N N 14 2 N 180 2 3 0
23 64 M 110 D Y N N Y N N N Y Y Y Y Y N N PACS N Y N Y N 14 3 Y 128 3 0 3
24 45 M 162 U N N N Y N N N N Y N Y N N N LACS Y Y Y N N 14 2 N 180 2 2 0
;
run;

/*********************************************
Step 2: Data Cleaning and Preparation
*********************************************/
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
Step 3: Cox Proportional Hazards Regression
(same PROC PHREG model as the repo's SASscript.sas Step 5, first block;
 ties=efron here in place of the source's ties=exact, since Jenner
 does not yet support TIES=EXACT)
*********************************************/

  proc phreg data=istclean;
   where FDEADC in (2, 3);
class RXASP(ref='Y')  RSLEEP(ref='N') RCT(ref='Y') RATRIAL(ref='N') RVISINF(ref='N') RCONSC(ref='F') SEX(ref='F') stroke_type(ref='2');
  model FDEADD*DEAD(0) = stroke_type  RXASP  RATRIAL RVISINF RCONSC RSLEEP SEX AGE/ ties=efron risklimits= wald;
strata FDEADC;
hazardratio stroke_type;
run;
