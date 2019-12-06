/*The following Code tracks all polymorphisms from RT region and separates them into columns*/
data HIVFull.RTPol;
	set HIVFull.C567HET1X HIVFull.C567MSM1X HIVFull.SINGLES1XALL2019 
		HIVFull.CNONB1X HIVFull.SNONB1X;
	keep Class Rt_poly lab_ID;

	if missing(Lab_ID) then
		delete;
run;

data long;
	set hivfull.rtpol;
	id=lab_ID;
	n_loop=countw(Rt_poly, ",");

	do i=1 to n_loop;
		polyMorph=scan(Rt_poly, i, ',');
		output;
	end;
run;
proc sort data=long;
	by Lab_ID;
run;

proc transpose data=long out=ordered (drop=_name_ _label_) prefix=Rt_polymorph;
	by lab_ID;
	var polymorph;
run;

proc sort data=ordered;
	by Lab_ID;
run;

proc sort data=hivfull.RiskGroups out=hivfull.sorted;
	by lab_ID;
run;

data hivfull.Rt_polymorph;
	merge hivfull.sorted work.ordered;
	by lab_id;

	if missing(Lab_ID) then
		delete;
run;
proc print dat=hivfull.rt_polymorph (obs=10);
run;

/*Using Standford DB to separate all NNRTI polymorphisms individually*/
data polytest;
	set HIVFull.rtpol;
	if (find(RT_Poly,'K 103','i') > 0 or find(RT_Poly,'103','i') > 0) then K103 = 1;
	else if find(RT_Poly,'K 101','i') > 0 or find(RT_Poly,'101','i') > 0 then K101 = 1;
	else if find(RT_Poly,'L 100','i') > 0 or find(RT_Poly,'100','i') > 0 then L100 = 1;
	else if find(RT_Poly,'A 98','i') > 0 or find(RT_Poly,'98','i') > 0 then A98 = 1;
	else if find(RT_Poly,'V 106','i') > 0 or find(RT_Poly,'106','i') > 0 then V106 = 1;
	else if find(RT_Poly,'V 108','i') > 0 or find(RT_Poly,'108','i') > 0 then V108 = 1;
	else if find(RT_Poly,'E 138','i') > 0 or find(RT_Poly,'138','i') > 0 then E138 = 1;
	else if find(RT_Poly,'V 179','i') > 0 or find(RT_Poly,'179','i') > 0 then V179 = 1;
	else if find(RT_Poly,'Y 181','i') > 0 or find(RT_Poly,'181','i') > 0 then Y181 = 1;
	else if find(RT_Poly,'Y 188','i') > 0 or find(RT_Poly,'188','i') > 0 then Y188 = 1;
	else if find(RT_Poly,'G 190','i') > 0 or find(RT_Poly,'190','i') > 0 then G190 = 1;
	else if find(RT_Poly,'H 221','i') > 0 or find(RT_Poly,'221','i') > 0 then H221 = 1;
	else if find(RT_Poly,'P 225','i') > 0 or find(RT_Poly,'225','i') > 0 then P225 = 1;
	else if find(RT_Poly,'F 227','i') > 0 or find(RT_Poly,'227','i') > 0 then F227 = 1;
	else if find(RT_Poly,'M 230','i') > 0 or find(RT_Poly,'230','i') > 0 then M230 = 1;
	else if find(RT_Poly,'K 238','i') > 0 or find(RT_Poly,'238','i') > 0 then K238 = 1;
	else if find(RT_Poly,'Y 318','i') > 0 or find(RT_Poly,'318','i') > 0 then Y318 = 1;
	else if find(RT_Poly,'N 348','i') > 0 or find(RT_Poly,'348','i') > 0 then N348 = 1;
run;


proc sort data=polytest;
	by Lab_ID;
run;
data hivfull.RT_StanfordDB (drop=RT_Poly);
	merge hivfull.sorted work.polytest;
	by lab_id;

	if missing(Lab_ID) then
		delete;

run;

data hivfull.RT_StanfordDB;
	set hivfull.RT_StanfordDB;
	where (YearInfection Between 2002 AND 2017);
	if 2002=< YearInfection <=2005 then
		Period="2002-2005";
	else if 2006=< YearInfection <=2009 then
		Period="2006-2009";
	else if 2010=< YearInfection <=2013 then
		Period="2010-2013";
	else if 2014=< YearInfection <=2017 then
		Period="2014-2017";

/*End of DATA Step creation*/


proc print data=HIVFull.RT_STANFORDDB;
	Where (YearInfection Between 2002 AND 2017);
	sum V98 K103 K101 L100 V106 V108 E138 V179 Y181 Y188 G190 H221 P225 F227 M230 K238 Y318 N348;
run;
proc sort data=HIVFULL.RT_STANFORDDB;
	by Period;
run;

proc freq data=HIVFULL.RT_STANFORDDB;
Where (YearInfection Between 2002 AND 2017) and Classification = 'Subtype B' and Clinical ^='CT';
table Class;
run;
title  "NNRTI Acquired Mutations and Natural Polymorphisms";
proc means data=HIVFull.RT_STANFORDDB N order=freq;
	Where (YearInfection Between 2002 AND 2017) and Classification = 'Subtype B' and Clinical ^='CT';
	*var A98 K103 K101 L100 V106 V108 E138 V179 Y181 Y188 G190 H221 P225 F227 M230 K238 Y318 N348;
	var K103 K101 L100 V106 V108 E138 V179 Y181 Y188 G190 H221 P225 F227 M230 K238 Y318 N348;
	class Class;
	*by Period;
	*output out=RTPol_Stats;
run;
footnote "Conditions for this: Years between 2002 and 2017, Subtype B only, no CT";
data RTPol_Stats;
	set RTPol_Stats;
	if Class = 'Large Cluster' then GroupOrder= 1;
	else if Class = 'Small Cluster' then GroupOrder=2;
	else if Class = 'Singleton' then GroupOrder = 3;
run;
proc sort data=work.RTPol_Stats out=Rt_PolStatsSorted;
	by GroupOrder;
run;
					
proc sgplot data=Rt_PolStatsSorted;
	vbar Period/ group=Class response=_Freq_ stat=pct groupdisplay=cluster grouporder=data;
	yaxis label = "Percentage of acquired NNRTI Resistance";
run;




/*Looking at K103*/
proc means data=HIVFull.RT_STANFORDDB N order=freq;
	Where Classification = 'Subtype B' and Clinical ^='CT';
	var K103;
	class Risk_group Class;
	by Period;
	output out=K103Stats (drop= _freq_) N= /Autoname;
run;
data K103Stats;
	set K103Stats;
	if Class = 'Large Cluster' then GroupOrder= 1;
	else if Class = 'Small Cluster' then GroupOrder=2;
	else if Class = 'Singleton' then GroupOrder = 3;
run;
proc sort data=work.k103stats out=K103Sorted;
	by GroupOrder;
run;
proc sgplot data=work.K103Sorted;
	vbar Period/ group=Class stat=pct response=K103_N groupdisplay=cluster grouporder=data;
	yaxis label='K103N/S/H/T/R/Q/E Mutations';
	inset "People with K103 mutations: 275"/ position=topright;
run;

proc means data=HIVFull.RT_STANFORDDB N order=freq;
	Where Classification = 'Subtype B' and Clinical ^='CT';
	var A98;
	class Risk_group Class;
	by Period;
	output out=A98Stats (drop= _freq_) N= /Autoname;
run;




/*K101*/
proc means data=HIVFull.RT_STANFORDDB N order=freq;
	Where Classification = 'Subtype B' and Clinical ^='CT';
	var K101;
	class Risk_group Class;
	by Period;
	output out=K101Stats (drop= _freq_) N= /Autoname;
run;
data K101Stats;
	set K101Stats;
	if Class = 'Large Cluster' then GroupOrder= 1;
	else if Class = 'Small Cluster' then GroupOrder=2;
	else if Class = 'Singleton' then GroupOrder = 3;
run;
proc sort data=work.K101stats out=K101Sorted;
	by GroupOrder;
run;
proc sgplot data=work.K101Sorted;
	vbar Period/ group=Class stat=pct response=K101_N groupdisplay=cluster grouporder=data;
	yaxis label='K101E/H/P/Q/R/N Mutations';
	inset "People with K101 mutations: 106"/ position=topright;
run;



proc means data=HIVFull.RT_STANFORDDB N order=freq;
	Where Classification = 'Subtype B' and Clinical ^='CT';
	var A98;
	class Risk_group Class;
	by Period;
	output out=A98Stats (drop= _freq_) N= /Autoname;
run;



/*Looking at V179*/
proc means data=HIVFull.RT_STANFORDDB N order=freq;
	Where Classification = 'Subtype B' and Clinical ^='CT';
	var V179;
	class Risk_group Class;
	by Period;
	output out=V179Stats (drop= _freq_) N= /Autoname;
run;
data V179Stats;
	set V179Stats;
	if Class = 'Large Cluster' then GroupOrder= 1;
	else if Class = 'Small Cluster' then GroupOrder=2;
	else if Class = 'Singleton' then GroupOrder = 3;
run;
proc sort data=work.V179stats out=V179Sorted;
	by GroupOrder;
run;
proc sgplot data=work.V179Sorted;
	vbar Period/ group=Class stat=pct response=V179_N groupdisplay=cluster grouporder=data;
	yaxis label='V179D/E/F/I/L/T Polymorphism';
	inset "People with V179 Polymorphism: 257"/ position=topright;
run;

/*Looking at V106*/
proc means data=HIVFull.RT_STANFORDDB N order=freq;
	Where Classification = 'Subtype B' and Clinical ^='CT';
	var V106;
	class Risk_group Class;
	by Period;
	output out=V106Stats (drop= _freq_) N= /Autoname;
run;
data V106Stats;
	set V106Stats;
	if Class = 'Large Cluster' then GroupOrder= 1;
	else if Class = 'Small Cluster' then GroupOrder=2;
	else if Class = 'Singleton' then GroupOrder = 3;
run;
proc sort data=work.V106stats out=V106Sorted;
	by GroupOrder;
run;
proc sgplot data=work.V106Sorted;
	vbar Period/ group=Class stat=pct response=V106_N groupdisplay=cluster grouporder=data;
	yaxis label='V106 Mutation';
	inset "People with V106 Mutation: 112"/ position=topright;
run;



/*Looking at G190*/
proc means data=HIVFull.RT_STANFORDDB N order=freq;
	Where Classification = 'Subtype B' and Clinical ^='CT';
	var G190;
	class Risk_group Class;
	by Period;
	output out=G190Stats (drop= _freq_) N= /Autoname;
run;
data G190Stats;
	set G190Stats;
	if Class = 'Large Cluster' then GroupOrder= 1;
	else if Class = 'Small Cluster' then GroupOrder=2;
	else if Class = 'Singleton' then GroupOrder = 3;
run;
proc sort data=work.G190stats out=G190Sorted;
	by GroupOrder;
run;
proc sgplot data=work.G190Sorted;
	vbar Period/ group=Class stat=pct response=G190_N groupdisplay=cluster grouporder=data;
	yaxis label='G190 Mutation';
	inset "People with G190 Mutation: 159"/ position=topright;
run;












/*Using Standford DB to separate all NRTI polymorphisms individually*/
data polytestNRTI;
	set HIVFull.rtpol;
	if (find(RT_Poly,'M 41','i') > 0 or find(RT_Poly,'41','i') > 0) then M41 = 1;
	else if find(RT_Poly,'A 62','i') > 0 or find(RT_Poly,'62','i') > 0 then A62 = 1;
	else if find(RT_Poly,'K 65','i') > 0 or find(RT_Poly,'65','i') > 0 then K65 = 1;
	else if find(RT_Poly,'D 67','i') > 0 or find(RT_Poly,'67','i') > 0 then D67 = 1;
	else if find(RT_Poly,'S 68','i') > 0 or find(RT_Poly,'68','i') > 0 then S68 = 1;
	else if find(RT_Poly,'T 69','i') > 0 or find(RT_Poly,'69','i') > 0 then T69 = 1;
	else if find(RT_Poly,'K 70','i') > 0 or find(RT_Poly,'70','i') > 0 then K70 = 1;
	else if find(RT_Poly,'L 74','i') > 0 or find(RT_Poly,'74','i') > 0 then L74 = 1;
	else if find(RT_Poly,'V 75','i') > 0 or find(RT_Poly,'75','i') > 0 then V75 = 1;
	else if find(RT_Poly,'F 77','i') > 0 or find(RT_Poly,'77','i') > 0 then F77 = 1;
	else if find(RT_Poly,'Y 115','i') > 0 or find(RT_Poly,'115','i') > 0 then Y115 = 1;
	else if find(RT_Poly,'F 116','i') > 0 or find(RT_Poly,'116','i') > 0 then F116 = 1;
	else if find(RT_Poly,'Q 151','i') > 0 or find(RT_Poly,'151','i') > 0 then Q151 = 1;
	else if find(RT_Poly,'M 184','i') > 0 or find(RT_Poly,'184','i') > 0 then M184 = 1;
	else if find(RT_Poly,'L 210','i') > 0 or find(RT_Poly,'210','i') > 0 then L210 = 1;
	else if find(RT_Poly,'T 215','i') > 0 or find(RT_Poly,'215','i') > 0 then T215 = 1;
	else if find(RT_Poly,'K 219','i') > 0 or find(RT_Poly,'219','i') > 0 then K219 = 1;
run;

proc sort data=polytestnrti;
	by Lab_ID;
run;
data hivfull.RT_StanfordDBNRTI (drop=RT_Poly);
	merge hivfull.sorted work.polytestnrti;
	by lab_id;

	if missing(Lab_ID) then
		delete;

run;
data hivfull.RT_StanfordDBNRTI;
	set hivfull.RT_StanfordDBNRTI;
	where (YearInfection Between 2002 AND 2017);
	if 2002=< YearInfection <=2005 then
		Period="2002-2005";
	else if 2006=< YearInfection <=2009 then
		Period="2006-2009";
	else if 2010=< YearInfection <=2013 then
		Period="2010-2013";
	else if 2014=< YearInfection <=2017 then
		Period="2014-2017";
/*End of DATA Step creation*/

proc freq data=HIVFULL.RT_STANFORDDBNRTI;
	Where (YearInfection Between 2002 AND 2017) and Classification = 'Subtype B' and Clinical ^='CT';
	table Class;
run;

title  "NNRTI Acquired Mutations and Natural Polymorphisms";
proc means data=HIVFull.rt_stanforddbnrti N order=freq;
	Where (YearInfection Between 2002 AND 2017) and Classification = 'Subtype B' and Clinical ^='CT';
	var M41 A62 K65 D67 S68 T69 K70 L74 V75 F77  M184 L210 T215 K219;
	class Class;
	*by Period;
	*output out=RTPol_Stats;
run;
footnote "Conditions for this: Years between 2002 and 2017, Subtype B only, no CT";


/*creates a column that counts the number of unique polymorphisms per row */
data Rtwant;
	set hivfull.Rt_polymorph;
	array Rt_polymorph{42}Rt_polymorph1-Rt_polymorph42;
	array new {42} $20  _temporary_;
	do _n_=1 to 42;
		new{_n_}=Rt_polymorph{_n_};
	end;
	call sortc(of new{*});
	count=(new{1}>'');

	do _n_=2 to 42;

		if new{_n_} ne new{_n_-1} then
			count + 1;
	end;
run;

proc print data=Rtwant (obs=5);
run;

proc sort data=RtWant out=RtWantOrdered;
by Class;
/*Polymorphism distribution*/
title 'Distribution of Polymorphisms in RT';
proc sgplot data=RtwantOrdered;
  Where Classification = 'Subtype B';
  histogram Count / group=Class transparency=0.5;  
  density Count / type=kernel group=Class;
  xaxis label= "Number of RT Polymorphisms" values=(0 to 42 by 1);
  styleattrs datacontrastcolors=(darkblue darkgreen darkred) datacolors=(lightblue lightgreen lightcoral) datalinepatterns=(1);
run;

/*RT Polymorphism ANOVA*/
Title;
ods noproctitle;
ods graphics / imagemap=on;

proc glm data=WORK.RTWANT plots(only)=(boxplot diagnostics);
	class Class;
	Where (YearInfection between 2002 and 2017) and Classification = "Subtype B";
	model count=Class;
	means Class / hovtest=levene welch plots=none;
	lsmeans Class / adjust=tukey pdiff alpha=.05 plots=(meanplot diffplot);
	run;
quit;