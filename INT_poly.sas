data HIVFull.INTPol;
	set HIVFull.C567HET1X HIVFull.C567MSM1X HIVFull.SINGLES1XALL2019 
		HIVFull.CNONB1X HIVFull.SNONB1X temp.int_tree;
	keep Class INT_poly Database_ID;

	if missing(Database_ID) then
		delete;
	if missing(int_poly) then 
		delete;
run;

data INTpolytest;
	set HIVFull.INTpol;
	if find(INT_Poly, '10','i') > 0 then D10 = 1;
	if find(INT_Poly, '11','i') > 0 then E11 = 1;
	if find(INT_Poly, '24','i') > 0 then S24 = 1;
	if find(INT_Poly, '20','i') > 0 then R20 = 1;
	if find(INT_Poly, '50 ','i') > 0 then M50 = 1;
	if find(INT_Poly, '51','i') > 0 then H51 = 1;
	if find(INT_Poly, '66','i') > 0 then T66 = 1;
	if find(INT_Poly, '92','i') > 0 then E92 = 1;
	if find(INT_Poly, '95','i') > 0 then Q95 = 1;
	if find(INT_Poly, '97','i') > 0 then T97 = 1;
	if find(INT_Poly, '118','i') > 0 then G118 = 1;
	if find(INT_Poly, '121','i') > 0 then F121 = 1;
	if find(INT_Poly, '138','i') > 0 then E138 = 1;
	if find(INT_Poly, '140','i') > 0 then G140 = 1;
	if find(INT_Poly, '143','i') > 0 then Y143 = 1;
	if find(INT_Poly, '146','i') > 0 then Q146 = 1;
	if find(INT_Poly, '147','i') > 0 then S147 = 1;
	if find(INT_Poly, '148','i') > 0 then Q148 = 1;
	if find(INT_Poly, '151','i') > 0 then V151 = 1;
	if find(INT_Poly, '153','i') > 0 then S153 = 1;
	if find(INT_Poly, '155','i') > 0 then N155 = 1;
	if find(INT_Poly, '157','i') > 0 then E157 = 1;
	if find(INT_Poly, '163','i') > 0 then G163 = 1;
	if find(INT_Poly, '230','i') > 0 then S230 = 1;
	if find(INT_Poly, '263','i') > 0 then R263 = 1;
run;


proc sort data=INTpolytest;
	by Database_ID;
run;
proc sort data=hivfull.sorted;
	by Database_ID;
run;
data hivfull.INT_StanfordDB (drop=INT_Poly);
	length class $30;
	merge  work.INTpolytest hivfull.sorted;
	by Database_id;
	
	if missing(Database_ID) then
		delete;

run;
data hivfull.INT_StanfordDB;
	set hivfull.INT_StanfordDB;
	where (YearInfection Between 2002 AND 2017);
	if 2002=< YearInfection <=2005 then
		Period="2002-2005";
	else if 2006=< YearInfection <=2009 then
		Period="2006-2009";
	else if 2010=< YearInfection <=2013 then
		Period="2010-2013";
	else if 2014=< YearInfection <=2017 then
		Period="2014-2017";
/*End of Data Step*/
proc sort data=HIVFULL.INT_STANFORDDB;
	by Period;
run;
title  "INI Acquired Mutations and Natural Polymorphisms";
proc means data=hivfull.INT_StanfordDB N order=freq;
	Where (YearInfection Between 2002 AND 2017) and Classification = 'Subtype B' and Clinical ^='CT';
	var D10 E11 S24 R20 M50 H51 T66 E92 Q95 T97 E138 G140 Q146 Q148 V151 N155 E157 G163 S230;
	*var D10 E11 S24 R20 M50 H51 T66 E92 Q95 T97 G118 F121 E138 G140 Y143 P145 Q146 S147 Q148 V151 S153 N155 E157 G163 S230 R263;
	class class;
	*by Period;
	*output out=INTPol_Stats;
run;

proc sql;
	SELECT D10, E11, S24, R20, M50, H51, T66, E92, Q95, T97, E138, G140, Q146, Q148, V151, N155, E157, G163, S230
	FROM hivfull.INT_StanfordDB;
	WHERE Cluster_ID in (C045, C118, C185);
quit;


title  "C185 Poly";
proc means data=HIVFull.RT_STANFORDDB N order=freq;
	Where Cluster_ID = 'C185' and Clinical ^='CT';
	var A98 K103 K101 L100 V106 V108 E138 V179 Y181 Y188 G190 H221 P225 F227 M230 K238 Y318 N348;
	*var K103 K101 L100 V106 V108 E138 V179 Y181 Y188 G190 H221 P225 F227 M230 K238 Y318 N348;
	class Risk_group Class;
	by Period;
	output out=RTPol_Stats;
run;
title  "C185 NRTI Poly";
proc means data=HIVFull.rt_stanforddbnrti N order=freq;
	Where Cluster_ID = 'C185' and Clinical ^='CT';
	var M41 A62 K65 D67 S68 T69 K70 L74 V75 F77  M184 L210 T215 K219;
	class Class;
	*by Period;
	*output out=RTPol_Stats;
run;
title  "C185 PRI Poly";
proc means data=hivfull.PR_StanfordDB N order=freq;
	Where Cluster_ID = 'C185'and Clinical ^='CT';
	var K20 L23 L24 D30 V32 L33 K43 M46 I47 G48 I50 F53 I54 Q58 G73 T74 V82 N83 L89 L90;
	class Class;
	*by Period;
	*output out=RTPol_Stats;
run;