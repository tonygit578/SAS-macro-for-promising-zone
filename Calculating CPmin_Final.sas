**--------------------------------------------------------------------------------------------------------------
** Promising Zone Designs for Sample Size Re-estimation in Clinical Trials: Graphical Communication using SAS---
** PharmSUG paper by Zhao Yang and Shivani Nanda ---------------------------------------------------------------
** It is a SAS implementation to obtain the CPmin based on promising zone design proposed by Mehta et al (2011)-
** --------------------------------- Details for macro parameters ---------------------------------------------
	n1: 	total # of events/sample size at the interim analysis 
	n2: 	total # of events/sample size at final analysis based on the original design 
	cpmax: 	expected maximal conditional power 
	nmax: 	the total affortable maximal # of events/sample size at final analysis 
	alpha: 	target two-sided significance level
	dsout: 	the final datasets name to be produced 
**-------------------------------------------------------------------------------------------------------------;

%macro FindingCPmin(n1 =, n2 =, cpmax = 0.9, nmax =, alpha = 0.05, dsout = );
	data KeyQuanti;
		/* this is to convert the macro parameter to variables */
		z_alpha1s = quantile('Normal', 1 - &alpha/2, 0, 1);
		z_cpmax	  = quantile('Normal', &cpmax, 0, 1);

		/* make a transformation to HR, hence, it follows a normal distribution */
		n1		  = &n1;
		n2		  = &n2;
		nmax	  = &nmax;
		%global alphaRef;
		call symput("alphaRef", z_alpha1s);
	run;

	/* iterating different value of Z1 which can be the potential z statistic calcualted at IA */
	data PotentialZstatatIA;
		set KeyQuanti;
		do z1 = -1 to 4 by 0.0001;
			/*  this is to calculate the conditional power based on the observed 
				statistic Z1 at IA assuming the original design is still followed */
			V4CP = z1 * sqrt((n2 - n1)/n1) - (z_alpha1s * sqrt(n2) - sqrt(n1)*z1)/sqrt(n2 - n1);
			condiPowerAtIA = probnorm(V4CP); 

			/* this is the incremental sample size formula from Section 5.1 */
			n2_inc_new  = (n1/z1**2)*( (z_alpha1s*sqrt(n2) - z1*sqrt(n1) )/sqrt(n2 - n1) + z_cpmax)**2; 
			n2_new 	    = n1 + n2_inc_new;
			n2_real     = min(n2_new, nmax);
			n2_real_inc = n2_real - n1;
			
			/* this is to calculate the b value in Section 5.2 */
			b = n2_real**(-0.5)*( (sqrt(n2_real_inc/(n2 - n1)))*(z_alpha1s*sqrt(n2) - z1*sqrt(n1)) + z1*sqrt(n1));
			
			if b <= z_alpha1s then flag = 1; else flag = 0;
			output;
		end;
	run;
	
	/* this is to find the record with the CPmin and CPmax */
	proc sort data = PotentialZstatatIA;
		by flag z1;
	data FindRecCPmin;
		set PotentialZstatatIA;
		by flag z1;
		if (first.flag or last.flag) & flag = 1 then output;
	run;
	proc sql noprint;
	   select condiPowerAtIA into: cp separated by ','
	   from FindRecCPmin; quit; 

	%global n_cpmin n_cpmax;
	%let n_cpmin = %sysfunc(strip(%sysfunc(scan(%bquote(&cp),1, ','))));
	%let n_cpmax = %sysfunc(strip(%sysfunc(scan(%bquote(&cp),2, ','))));
	
	/* this is to print the identified record */
	proc print data = FindRecCPmin noobs;
		var z_alpha1s n1 n2 nmax z1 condiPowerAtIA n2_new n2_real;
	run;

	proc sort data = PotentialZstatatIA;
		by z1;
	data &dsout;
		set PotentialZstatatIA;
		Nobs + 1;
	run;
%mend;

%FindingCPmin(n1 = 140, n2 = 280, cpmax = 0.9, nmax = 420, alpha = 0.05, dsout = PrepData4CPmin);

**-------------------------------------------------------------------------------------------------------------
					this is to create the graph with CP as the x-axis, and b as y-axis
**-------------------------------------------------------------------------------------------------------------;
/* this is the annotation dataset for axis label, note the key variable is drawspace = "datavalue" */
data label4axis;
	length function $4. label $400.;
	function = "text"; x1 = 0.5; y1 = 1.765; label = "CP(z^{sub '1'}, n^{unicode tilde}^{sub '2'}) evaluated at ^{unicode delta}^{unicode hat}^{sub '1'}"; 
						drawspace = "datavalue"; width = 1000; textweight='bold'; output;
	function = "text"; x1 = -0.1; y1 = 1.95; label = "b(z^{sub '1'}, n^{unicode tilde}^{sub '2'}^{sup '*'}(z^{sub '1'}))"; 
						drawspace = "datavalue"; rotate = 90; width = 1000; textweight = 'bold'; output;

	/* this is to label the CPmin and its value, please note that since rotate = 90 is used in the previous
		statement, its effect is still there unless we want to change */
	function = "text"; x1 = &n_cpmin - 0.02; y1 = 1.85; label = "CP^{sub 'min'} = &n_cpmin"; 
						drawspace = "datavalue"; width = 1000; textweight = 'bold'; output;  
run;

ods escapechar = '^';

ods listing gpath = 'G:\Papers for Publication\PharmSUG\Promising zone\Software\Output';
ods graphics / imagename = "FindingCPmin" imagefmt = png; 
proc sgplot data = PrepData4CPmin noautolegend pad = (top = 1% bottom = 1% left = 1% right = 1%) sganno = label4axis;
	/* the conditional power based on original design is graphed */
	series x = condiPowerAtIA y = b / lineattrs = (color = black thickness = 0.12cm pattern = solid);

	/* now define the x-axis and the grid lines */
	xaxis label = ' ' min = 0 max = 1 values = (0 to 1 by 0.1);  
	refline (0 to 1 by 0.1) / axis = x lineattrs = (pattern = dash color = gray) transparency = 0.1;

	/* now define the y-axis and the grid lines */
	yaxis min = 0 label = ' ' values = (1.8 to 2.1 by 0.05);
	refline (1.8 to 2.1 by 0.05) / axis = y lineattrs = (pattern = dash color = gray) transparency = 0.1;

	/* this is to define the reference line based on alpha */
	refline &alphaRef / axis = y lineattrs = (thickness = 0.12cm color = blue pattern = dash); 

	/* this draw the shaded region between CPmin and CPmax, the first shade the area, and the 2nd define the region*/
	band y = b lower = &n_cpmin upper = &n_cpmax / transparency = 0.7;
	refline &n_cpmin &n_cpmax / axis = x lineattrs = (thickness = 1 color = darkgray pattern = solid); 
run;
