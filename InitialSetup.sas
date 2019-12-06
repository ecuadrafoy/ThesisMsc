/*This code was used to import data from Excel into SAS libraries*/

Options validvarname=V7;
libname HIVFull XLSX '/home/ernestocuadrafoy0/HIV/HivFull/May2019Combined.xlsx';
data HIVFull.Fullycombined;
	set HIVFull.C567HET1X HIVFull.C567MSM1X HIVFull.SINGLES1XALL2019 HIVFull.CNONB1X HIVFull.SNONB1X;
	DateVisit=input(Visit_Date, 10.);
	DateVisit=DateVisit - 21916;
	NewAge=input(Age, 2.);
	State=propcase(State);
	*fixes lowercase state;
	YearInfection=Year(DateVisit);
	ViralLoad1=input(ViralLoad, 10.);
	logVL=log10(ViralLoad1);
	
	keep Patient_ID Pct_Div Database_ID Cluster_ID Clinical Sex Lab_ID Region State NewAge 
		ViralLoad1 Class Classification DateVisit YearInfection Subtype logVL Visit_Date;
	format DateVisit ddmmyy10.;
	region = tranwrd(region, 'Montérégie - CRÉ de Longueuil', 'Monteregie');
	region = tranwrd(region, 'Montérégie', 'Monteregie');
	region = tranwrd(region, 'Montérégie - CRÉ Montérégie Est', 'Monteregie');
	region = tranwrd(region, 'Monteregie - CRÉ Monteregie Est', 'Monteregie');
run;
/*The following Data Step adds the column age group according to the quartiles*/
data hivfull.ageclass;
	set hivfull.riskgroups;
	Where NewAge >=15;

	if 15 < NewAge <=31 then
		AgeGroup="15-31";
	else if 31 < NewAge <=47 then
		AgeGroup="32-47";
	else if NewAge > 47 then
		AgeGroup=">47";
	keep Lab_ID class newage YearInfection logVL Classification AgeGroup PCT_Div 
		Cluster_ID Risk_group Clinical State;
run;
/*This step adds the column for Historical Clusters and Active Cluters*/
data hivfull.ActiveClusters;
	set hivfull.ageclass;
	Where (YearInfection Between 2002 AND 2017);

	if 2002=< YearInfection <=2013 then
		YearGroup="Historical Clusters";

	if 2013 < YearInfection then
		YearGroup="Active Clusters";
run;
data hivfull.QuadrennialYears;
	set hivfull.ActiveClusters;
	where (YearInfection Between 2002 AND 2017);

	if 2002=< YearInfection <=2005 then
		Period="2002-2005";
	else if 2006=< YearInfection <=2009 then
		Period="2006-2009";
	else if 2010=< YearInfection <=2013 then
		Period="2010-2013";
	else if 2014=< YearInfection <=2017 then
		Period="2014-2017";
run;
/*This step cleans up the Clinical column and standardizes the terminology*/
data hivfull.clinical;
	set hivfull.QuadrennialYears;
	where Clinical in ("CUN","P","CT","PHI","6");
	if Clinical = "P" then Clinical = "PHI";
	else if Clinical = "6" then Clinical = "CUN";
run;
Options validvarname=V7;
libname Temp XLSX '/home/ernestocuadrafoy0/HIV/HIVRiskGroup.xlsx';

data hivfull.RiskG (keep=Lab_ID Risk_Group);
	set temp.c567iduhet1x temp.msmc5671x temp.siduhet1x temp.smsm1x;
run;

proc sort data=hivfull.fullycombined out=hivfull.sorted;
	by  Lab_ID;
run;

proc sort data= hivfull.RiskG;
	by  Lab_ID;
run;

data hivfull.RiskGroups;
	merge hivfull.riskg hivfull.sorted;
	by Lab_ID;
run;
	