*set working dictionary;
%let path = C:\Users\jiangy12\Documents\Group project /*path to the working folder*/;

proc import out = df 
	datafile = "&path.\colon_tidy.csv"
	dbms=csv replace;
	getnames=yes;
run;

libname myFiles "&path.\";

proc print data=df;run;


proc format;
	value rx_f 1 = "Observation" 
			   2 = "Levamisole" 
			   3 = "Levamisole+5-FU";


%macro stats_sur(class,var,label,order);
proc lifetest data=df;
	time tt_&var*&var(0);
	strata &class;
	ods output ProductLimitEstimates=ple_&var;
run;


data ple_&var;
	set ple_&var;
	if survival = . then delete;
run;


proc sql;
	CREATE TABLE out_&var AS
	SELECT failed, failure
	FROM (
    	SELECT MAX(failed) AS failed, MAX(failure) AS failure
    	FROM ple_&var
    	GROUP BY stratum) as f;
quit;

data out_&var;
	set out_&var;
	failed1 = put(failed,6.0);
	failure1 = cats("(",put(failure*100,6.1),"%)");
	value = catx(" ", failed1,failure1);
	drop failed failure failed1 failure1;
run;

proc phreg data=df;
	class &class(ref = '1');
	model tt_&var*&var(0)=&class;
	hazardratio &class / alpha=0.05 cl=wald diff=ref;
    ods output ParameterEstimates=estimate_&var HazardRatios=hr_&var;
run; 

* add hazard ratio;
data hr_&var;
	set hr_&var (keep = HazardRatio WaldLower WaldUpper firstobs =1 obs=2);
	hr = put(HazardRatio,5.2);
	ci = "("||trim(left(put(WaldLower,5.2)))||" - "||trim(left(put(WaldUpper,5.2)))||")";
	hr_ci = hr||" "||ci;
	row_num = 2;
	drop hr ci HazardRatio WaldLower WaldUpper;
run;

proc sql;
	insert into hr_&var
		values("Ref.", 1);
quit;

proc sort data=hr_&var;
	by row_num;
run;
	
data out_&var;
	merge out_&var hr_&var (keep = hr_ci);
run;

* add p-value;
data estimate_&var;
	set estimate_&var (keep =ProbChiSq firstobs =1 obs=2);
	pvalue = put(ProbChiSq,pvalue5.3); 
	row_num = 2;
	drop ProbChiSq;
run;

proc sql;
	insert into estimate_&var
		values(" ", 1);
quit;

proc sort data=estimate_&var;
	by row_num;
run;
	
data out_&var;
	merge out_&var estimate_&var (keep = pvalue);
	num = _N_;
run;

* add label;
proc sql;
	insert into out_&var
		values(" "," "," ",0);
quit;

proc sort data=out_&var;
	by num;
run;

data out_&var;
	set out_&var;
	length label $50;
	if num = 0 then label = &label;
	else label = put(num, rx_f.);
	drop num;
run;

data out&order;
	set out_&var;
run;
%mend;



%stats_sur(rx, death,"Death",1);
%stats_sur(rx, recur,"Recurrence",2);
%stats_sur(rx, comp,"Death, or Recurrence",3);



data alldata;
	format label value hr_ci pvalue;
	set out1 out2 out3;
run;



ods rtf file="&path.\Table3. Overall associations between three treatment groups and adverse events.rtf";
ods escapechar="^";
ods rtf text="^S={width=100% just=c fontweight=bold fontsize=12pt} Table 3. Overall associations between three treatment groups and adverse events";
options nodate nonumber center; title; 
proc report data=alldata style=[background = white outputwidth=80%] style(report)={outputwidth=80% frame = hsides rules = group};
	column label value hr_ci pvalue;
	define label/" " style(header)=[verticalalign=center just=left background=white];
	define value/"Event(%)" style(column)=[cellwidth=1in just=center] style(header)=[verticalalign=center background=white];
	define hr_ci/"Hazard ratio ^n (95% CI)" style(column)=[cellwidth=1in just=center] style(header)=[verticalalign=center background=white];
	define pvalue/"p-value" style(header)=[verticalalign=center background=white] style(column)=[just=center];
	compute label;
   		if label="Observation" or label="Levamisole" or label="Levamisole+5-FU" then do;
    			call define(_COL_,'style','style={leftmargin=.20in}');
   		end;
		if label="Death" or label="Recurrence" or label="Death, or Recurrence" then do;
    			call define(_COL_,'style','style={font_weight=bold}');
   		end;
	endcomp;
run;
ods rtf close;
