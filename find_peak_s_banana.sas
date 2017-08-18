
/*Finds peak s for banana dataset*/


libname dataset "U:\dev\FISVDD\data";


data trainSet;
	set banana;
run;

/*SVDD WITH LOOP FOR MULTIPLE S VALUES;*/

/*Setup CAS*/
%let sysparm=runcas:ls+pbox:nd+host:rdcgrd033;
%cassetup; 

data sascas1.trainds;
set trainSet;
run;

%let sStart = 0.5;
%let sby=0.01;

/*Run PROC SVDD on different s*/
%macro _call;

/*Prepare the summary table*/

proc sql noprint;
   create table work.summary
       (
        s num,
        sv num,
        radius num,
		last_term num);
quit;

/*Run PROC SVDD on different s*/

%do i=1 %to 70;
    %let s=%sysevalf (&sStart.+&i.*&sby.-&sby.);
    ods listing close;
    ods output trainingresults=tr;
    proc svdd data=sascas1.trainds outsv=sascas1.dk; 
    	input x1 x2; 
        kernel rbf / bw=&s.;
		savestate rstore=sascas1.state_s;
    run;

%local nsv;
%local radius;
%local last_term;

/*Obtain information from each run*/

proc sql noprint;
	select value into :nsv
    from tr
    where Description = "Number of Support Vectors";

    select value into :radius
    from tr
    where Description = "Threshold R^2 Value";

	select value into :last_term
    from tr
    where Description = "Constant (C_r) Value";


    insert into work.summary
    values(&s, &nsv, &radius, &last_term);
quit;
%end;
%mend _call;

%_call;

/*Close CAS*/
%casclear;










/*FIND PEAK S USING SECOND DERIV OF OPTIMAL OBJECTIVE FUNCTION*/
%macro a(inds=,sby=);

/* Take the second derivative of last_term */
data dk;
    set &inds;
    l1_last_term=lag(last_term);
    l2_last_term=lag2(last_term);
    if _n_ > 1 then d_last_term = (last_term - l1_last_term)/(&sby);
    if _n_ > 2 then d2_last_term = (last_term - 2 * l1_last_term + l2_last_term) / ( &sby * &sby ) ;
run;


data dk;
	set dk;
	if not missing(d2_last_term);
run;


/* Define custom appearance. */
%let ds=dk;
ODS PATH work.templat(update) sasuser.templat(read) sashelp.tmplmst(read);
    proc template;
       define style MyStyleDefault;
       parent=Styles.Journal;

            style graphfonts from graphfonts /
            'GraphDataFont' = ("<sans-serif>, <MTsans-serif>",7pt)
            'GraphUnicodeFont' = ("<MTsans-serif-unicode>",9pt)
            'GraphValueFont' = ("<sans-serif>, <MTsans-serif>",9pt)
            'GraphLabel2Font' = ("<sans-serif>, <MTsans-serif>",14pt,bold)
            'GraphLabelFont' = ("<sans-serif>, <MTsans-serif>",14pt,bold)
            'GraphFootnoteFont' = ("<sans-serif>, <MTsans-serif>",10pt)
            'GraphTitleFont' = ("<sans-serif>, <MTsans-serif>",14pt,bold)
            'GraphTitle1Font' = ("<sans-serif>, <MTsans-serif>",14pt,bold)
            'GraphAnnoFont' = ("<sans-serif>, <MTsans-serif>",10pt);

             style GraphData1 from GraphData1 /
             markersymbol = "circle"
             linestyle = 2
             contrastcolor = CXA6A6A6
             color = CXA6A6A6;

             style GraphData2 from GraphData2 /
             markersymbol = "circlefilled"
             linestyle = 2
             contrastcolor = CX000000
             color = CX000000;

             style GraphData3 from GraphData3 /
             markersymbol = "circlefilled"
             linestyle = 2
             contrastcolor = CXF2F2F2
             color = CXF2F2F2;
        end;
    run;

    proc template;
       define style MyStyleDefault1;
       parent=MyStyleDefault;
             style GraphData1 from GraphData1 /
             markersymbol = "circlefilled"
             linestyle = 2
             contrastcolor = CXA6A6A6
             color = CXA6A6A6;
        end;
    run;

proc transreg data=&ds;
   model identity(d_last_term) = pbspline(s / aicc)/alpha=0.05;
   output out=o predicted cli clm;
run;

data o2;
	set &ds;
run;

ods latex close;
ods listing;

ods graphics on /imagename="Peak Criterion" width=10.2in;
ods listing close;
ods latex style= MyStyleDefault;

proc transreg data=&ds;
   model identity(d2_last_term) = pbspline(s / aicc)/alpha=0.05;
   output out=o predicted cli clm;
run;

data o;
    set o(where=(d2_last_term ne .));
     if ((CMUd2_last_term ge 0) and (CMLd2_last_term le 0)) then zero_flg=1;
    else zero_flg=0;
run;

%local firsts;
%local lasts;

proc sql noprint;
    select  min(s) into :firsts
    from o
    where zero_flg=1 and d2_last_term ne .;
quit;

proc sql noprint;
    select min(s) into :lasts
    from o
    where zero_flg=0 and s gt &firsts and d2_last_term ne .;
quit;

%put &firsts &lasts;
%let lasts=%sysevalf(&lasts.-&sby.);
%put lasts is &lasts;

data o;
    set o;
    length lable $10;
    label d2_last_term = "Second Derivative";
    if round(s,0.0001)=&firsts or round(s,0.0001)=&lasts then do;
        ss=s;
    label="s="||trim(left(s));
    end;
run;

    ods latex close;
    ods listing;

    ods graphics on /imagename="out" width=10.2in;
    ods listing close;
    ods latex style= MyStyleDefault  ;
    %let psmax=8;
    %let psmin=0.001;
    %local ticks;
    %let ticks = %sysevalf((&psmax-&psmin) / 14); /* 15 tick marks */

proc iml;
    vec=rowvec(do(&psmin,&psmax,&ticks));
    vec = putn(vec,"BEST6.");
    N = ncol(vec);
    string = "";
    do i = 1 to N-1;
       string = string + strip(vec[i]) + " ";
    end;
    string = string + strip(vec[N]);
    call symputx("ticklist", string, 'l');
quit;

	%put &ticklist.;

    %let ticklist=0 1 2 3 4 5 6 7 8;

  proc template;
      define statgraph bandplot;
        begingraph;
          entrytitle " ";
            layout overlay/xaxisopts= ( tickvalueattrs=(size=18 weight=bold) labelattrs=(size=18 weight=bold)
                          type=linear  linearopts=(viewmin=&psmin. viewmax=&psmax. tickvaluelist=(%str(&ticklist))))
                          yaxisopts= ( tickvalueattrs=(size=18 weight=bold) labelattrs=(size=18 weight=bold) );
/*                          type=linear linearopts=(viewmin=-0.6 viewmax=.2 tickvaluelist=(%str(-0.6 -0.4 -0.2 0 0.2))));*/
            bandplot x=s limitupper=CMUd2_last_term
                    limitlower=CMLd2_last_term /
                    name="band1" modelname="fit"
                    legendlabel="95% Confidence Limits";
            scatterplot x=s y=d2_last_term / primary=true;
            seriesplot x=s y=Pd2_last_term / name="fit"
            legendlabel="Fit Line" lineattrs=(thickness=2);;
            referenceline x=ss / curvelabel=label
            lineattrs=(color=gray pattern=dot thickness=3) curvelabelattrs=(size=18 weight=bold);
            referenceline y=0 /
            lineattrs=(color=gray pattern=solid );
            discretelegend  "fit" "band1"/valueattrs=(size=18 weight=bold) ;
          endlayout;
        endgraph;
      end;
    run;



	proc sgrender data=o(where=(d2_last_term ne .)) template=bandplot;
    run;
%mend a;

options mprint mlogic symbolgen;
data mySummary (keep= s radius);
	set summary;
run;


data mySummary2 (rename=(radius=last_term));
	set mySummary;
run;

%a(inds=mySummary2,sby=0.01);

