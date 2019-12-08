proc sgplot data=hivfull.riskgroups;
	Where (YearInfection Between 2002 AND 2017) and (Risk_group = "MSM");

  	density Pct_Div / type=normal group=Class ;
  	xaxis label="Percent Diversity" values= (0 to 4 by 0.2);
  	inset ("Infected Patients" = "7,453") / border;
run;