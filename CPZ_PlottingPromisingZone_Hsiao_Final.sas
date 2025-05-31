**-------------------------------------------------------------------------------------------------------------
** Promising Zone Designs for Sample Size Re-estimation in Clinical Trials: Graphical Communication using SAS--
** PharmSUG paper by Zhao Yang and Shivani Nanda --------------------------------------------------------------
** It is a SAS implementation to the optimal promising zone design proposed by Hsiao et al (2019) -------------
** --------------------------------- Details for macro parameters ---------------------------------------------
	n1: 	total # of events/sample size at the interim analysis 
	n2: 	total # of events/sample size at final analysis based on the original design 
	cpmin: 	expected minimal conditional power 
	cpmax: 	expected maximal conditional power 
	hrmin: 	the minimal clinically meaningful treatment effect, here it is HR 
	nmax: 	the total affortable maximal # of events/sample size at final analysis 
	alpha: 	target two-sided significance level
	randR: 	randomization ratio between treatment vs control, i.e. ratio of treatment/control
	dsout: 	the final datasets name to be produced 
**-------------------------------------------------------------------------------------------------------------;

%macro Data4PZGraph(n1 =, n2 =, cpmin = 0.8, cpmax = 0.9, hrmin =, nmax =, alpha = 0.05, randR = 1, dsout = );
	data KeyQuanti;
		/* this is to convert the macro parameter to variables */
		z_alpha1s = quantile('Normal', 1 - &alpha/2, 0, 1);
		z_cpmin	  = quantile('Normal', &cpmin, 0, 1);
		z_cpmax	  = quantile('Normal', &cpmax, 0, 1);
		/* make a transformation to HR, hence, it follows a normal distribution */
		trt_min	  = -log(&hrmin);  
		n1		  = &n1;
		n2		  = &n2;
		nmax	  = &nmax;
		/* 	this is to convert the randomization ratio to a ratio which can  
		 	be used in the calculation of converting Z statistic to HR 		*/
		r		  = &randR/(&randR + 1);  

		/* this is to calculate the key quantities to be used in the PZ graph production based on Z1 scale */
		z1_L	  = round(z_alpha1s * sqrt(n2/n1) - sqrt((n2 - n1)/n1) * (trt_min * sqrt(nmax - n1)/2 - z_cpmin), 0.0001);
		z1_C	  = round(z_alpha1s * sqrt(n2/n1) - sqrt((n2 - n1)/n1) * (trt_min * sqrt(nmax - n1)/2 - z_cpmax), 0.0001);
		z1_U	  = round(z_alpha1s * sqrt(n2/n1) - sqrt((n2 - n1)/n1) * (trt_min * sqrt(n2   - n1)/2 - z_cpmax), 0.0001);

		/* this is to calculate the key quantities to be used in the PZ graph production based on HR scale */
		HR_U 	  = exp(-z1_L/(sqrt(n1*(r*(1 - r))))); ** smallest Z1 value corresponds to highest HR value;
		HR_C 	  = exp(-z1_C/(sqrt(n1*(r*(1 - r))))); 
		HR_L 	  = exp(-z1_U/(sqrt(n1*(r*(1 - r))))); ** highest Z1 value corresponds to smallest HR value;

		/* Declare global macros for those key quantities to be used later */
		%global LowerZ UpperZ CentrZ LowerHR UpperHR CentrHR;
		call symput("LowerZ", z1_L);
		call symput("UpperZ", z1_U);
		call symput("CentrZ", z1_C);

		call symput("LowerHR", HR_L);
		call symput("UpperHR", HR_U);
		call symput("CentrHR", HR_C);
	run;

	%put &LowerZ &UpperZ &CentrZ &LowerHR &UpperHR &CentrHR;

	/* iterating different value of Z1 which can be the potential z statistic calcualted at IA */
	data PotentialZstatatIA;
		set KeyQuanti;
		do z1 = -1 to 4 by 0.0001;
			/*  this is to calculate the conditional power based on the observed 
				statistic Z1 at IA assuming the original design is still followed */
			V4CP = trt_min * sqrt(n2 - n1)/2 - (z_alpha1s * sqrt(n2) - sqrt(n1)*z1)/sqrt(n2 - n1);
			condiPowerAtIA = probnorm(V4CP); 
			
			effHR = exp(-z1/(sqrt(n1*(r*(1 - r))))); ** this is to convert statistic Z1 at IA to HR;
			output;
		end;
	proc sort data = PotentialZstatatIA;
		by z1;
	data PotentialZstatatIA;
		set PotentialZstatatIA;
		Nobs + 1;
	run;

	data PrepData4PZ1;
		set PotentialZstatatIA;
		/* this is to break the conditional power from original design into 3 pieces 
		   it is the data to be used to generate the smooth red curve in the plot */
		if z1 < z1_L then condiPowerAtIA_P1 = condiPowerAtIA; * piece 1, lower part;
		if z1 > z1_U then condiPowerAtIA_P2 = condiPowerAtIA; * piece 2, upper part;
		if z1_L <= z1 & z1 <= z1_U then condiPowerAtIA_P3 = condiPowerAtIA; * piece 2, midle part;

		/** Now, prepare to get the relevant information due to the potential 
			sample size re-estimation. It is evaluated for different regions based on the cutoff
			values from z1_L z1_C and z1_U */
		
		/* First: check the region z1 <= z1_L */
		if z1 <= z1_L then do;
			/* the key characteristic is that at the end of the region, we boost the sample size to
				nmax, hence, we need to calculate the conditional power under this situation */
			V4CP_min = trt_min * sqrt(nmax - n1)/2 - (z_alpha1s * sqrt(n2) - sqrt(n1)*z1)/sqrt(n2 - n1);
			AdjcondiPower_min = probnorm(V4CP_min);

			** this variable will be used to pick the records at the end of region 
				based on the conditional power results for this region;
			if AdjcondiPower_min > 0 then nonmissMin = 1; else nonmissMin = 0; 

			n_SS = n2; ** it is noted that the sample for this region will remain the original planned n2;
		end;

		/* Second: check the region z1_L < z1 <= z1_C */
		if z1_L < z1 <= z1_C then do;
			/* The key characteristic is that we keep the sample size at nmax, 
				but increasing the value of Z1, then we calcualte the conditional power */
			V4CP_inc = trt_min * sqrt(nmax - n1)/2 - (z_alpha1s * sqrt(n2) - sqrt(n1)*z1)/sqrt(n2 - n1);
			AdjcondiPower_inc = probnorm(V4CP_inc);

			n_SS = nmax; ** the sample size is remained nmax for this region;
			n_SSOrig = n2; ** this will be used in the graph only since it is dashed in the plot;
		end;

		/* Third: check the region z1_C < z1 <= z1_U */
		if z1_C < z1 <= z1_U then do;
			/* the key characteristic is that we maintain the max conditional power */
			AdjcondiPower_dec = &cpmax;
			/* then, given the max conditional power and the z1 statistic, we calculate the sample size */
			n_SS = ( (2/trt_min) * ((z_alpha1s * sqrt(n2) - sqrt(n1)*z1)/sqrt(n2 - n1) + z_cpmax) ) **2 + n1;
			
			n_SSOrig = n2; ** this will be used in the graph only since it is dashed in the plot;
		end;

		/* Lastly: check the region z1_U < z1 */
		if z1_U < z1 then do;
			/* the only connecting point is at z1_U = z1 at which max conditiona power is maintained
				and the sample size now is the original planned sample size n2 */
			AdjcondiPower_dec = &cpmax; 
			n_SS = n2;
		end;

	proc sort data = PrepData4PZ1;
		by nonmissMin descending z1 AdjcondiPower_min;
	run;

	data PrepData4PZ2;
		set PrepData4PZ1;
		by nonmissMin descending z1 AdjcondiPower_min; 

		/* the main purpose is to pick two relevant records, the first one 
			is the record at z1 = z1_L; the 2nd one is the record that conditional power was calculated at
			the original planned sample size */
		if nonmissMin >. and first.nonmissMin then flagRecSSR = 1;
		if nonmissMin >. and lag1(first.nonmissMin) then flagRecOrig = 1;

	data &dsout;
		set PrepData4PZ2;

		/* then we obtain the updated conditional power based on sample size re-estimation 
		   activity, it is overall for the regaion of z1_L <= z1 <= z1_U */

		/* these two records is the conditional power at z1 = z1_L */
		if flagRecOrig = 1 then AdjcondiPower = condiPowerAtIA;
		if flagRecSSR = 1 then AdjcondiPower = AdjcondiPower_min;

		/* this is for the region z1_L < z1 <= z1_C */
		if z1_L < z1 <= z1_C then AdjcondiPower = AdjcondiPower_inc;

		/* this is for the region z1_C < z1 <= z1_U */
		if z1_C < z1 <= z1_U then AdjcondiPower = AdjcondiPower_dec;

		/* specify the label of two variables which will be used in the legend of the graph */
		label n_SS = "Sample Size" AdjcondiPower = "Conditional Power";

	proc sort data = &dsout;
		by z1;
	run;
%mend;
%Data4PZGraph(n1 = 140, n2 = 280, cpmin = 0.8, cpmax = 0.9, hrmin = 0.75, nmax = 420, alpha = 0.05, randR = 1, dsout = PrepData4PZ3);

/* printing the relevant information for the key quantities */
proc print data = KeyQuanti noobs;
	var n1 n2 nmax z1_L z1_C z1_U HR_L HR_C HR_U;
run;

/* this is the annotation dataset for Z statistic on x-axis, note the kay variable is drawspace = "datavalue"
	it will ensure that the macro variable of X1 can be used properly (i.e. the same scale as x-axis) */
data label4Z;
	length function $4. label $35.;
	function = "text"; x1 = &LowerZ; y1 = -0.05; label = "z^{sub '1'}^{sup 'L'}"; drawspace = "datavalue"; output;
	function = "text"; x1 = &CentrZ; y1 = -0.05; label = "z^{sub '1'}^{sup 'C'}"; drawspace = "datavalue"; output;
	function = "text"; x1 = &UpperZ; y1 = -0.05; label = "z^{sub '1'}^{sup 'U'}"; drawspace = "datavalue"; output;
run;

/* this is the annotation dataset for HR on x-axis, note the kay variable is drawspace = "datavalue"
	it will ensure that the macro variable of X1 can be used properly (i.e. the same scale as x-axis) */
data label4HR;
	length function $4. label $35.;
	function = "text"; x1 = &LowerHR; y1 = 0; label = "HR^{sup 'L'}"; drawspace = "datavalue"; output;
	function = "text"; x1 = &CentrHR; y1 = 0; label = "HR^{sup 'C'}"; drawspace = "datavalue"; output;
	function = "text"; x1 = &UpperHR; y1 = 0; label = "HR^{sup 'U'}"; drawspace = "datavalue"; output;
run;

ods escapechar = '^';

**-------------------------------------------------------------------------------------------------------------
							 this is to create the graph with Z1 as the x-axis
**-------------------------------------------------------------------------------------------------------------;
ods listing gpath = 'G:\Papers for Publication\PharmSUG\Promising zone\Software\Output';
ods graphics / imagename = "ConstrainedPZ_G4Z" imagefmt = png; 
proc sgplot data = PrepData4PZ3 noautolegend pad = (top = 1% bottom = 0% left = 1% right = 1%) sganno = label4Z;
	/* the conditional power based on original design is graphed via 3 series statements */
	series x = z1 y = condiPowerAtIA_P1 / lineattrs = (color = red thickness = 0.12cm pattern = solid);
	series x = z1 y = condiPowerAtIA_P2 / lineattrs = (color = red thickness = 0.12cm pattern = solid);
	series x = z1 y = condiPowerAtIA_P3 / lineattrs = (color = red thickness = 0.12cm pattern = dash);

	/* this is to graph the conditional power for region z1_L <= z1 <= z1_U, please note that the 
		option name = "RedS" will used in the statement of keylegend */
	series x = z1 y = AdjcondiPower / lineattrs = (color = red thickness = 0.12cm pattern = solid) name = "RedS" ;

	/* the following two statements are used to graph the sample size, also pay attention to the 
		option name = "BlackS" which will used in the statement of keylegend */
	series x = z1 y = n_SS / lineattrs = (color = black thickness = 0.12cm pattern = solid) y2axis name = "BlackS" ;
	series x = z1 y = n_SSOrig / lineattrs = (color = black thickness = 0.12cm pattern = dash) y2axis ;
	
	/* now define the x-axis and the grid lines */
	xaxis label = 'Z-statistic at Interim Analysis' min = -1 max = 4 values = (-1 to 4 by 1);  
	refline (-1 to 4 by 0.5) / axis = x lineattrs = (pattern = dash color = gray) transparency = 0.1;

	/* now define the y-axis and the grid lines */
	yaxis min = 0 label = 'Conditional Power' values = (0 to 1 by 0.2);
	refline (0 to 1 by 0.1) / axis = y lineattrs = (pattern = dash color = gray) transparency = 0.1;

	/* this is the second y-axis */
	y2axis min = 0 label = 'Sample Size' values = (0 to 1400 by 280);

	/* this define the legend, list item in the order you want them */	
	keylegend "RedS" "BlackS" / border down = 2 location = inside position = topleft opaque; 

	/* this draw the shaded region between z1_L and z1_U, the first shade the area, and the 2nd define the region*/
	band y = condiPowerAtIA_P1 lower = &LowerZ upper = &UpperZ / transparency = 0.7;
	refline &LowerZ &CentrZ &UpperZ / axis = x lineattrs = (thickness = 1 color = darkgray pattern = solid); 
run;


**-------------------------------------------------------------------------------------------------------------
							 this is to create the graph with HR as the x-axis
**-------------------------------------------------------------------------------------------------------------;
ods listing gpath = 'G:\Papers for Publication\PharmSUG\Promising zone\Software\Output';
ods graphics / imagename = "ConstrainedPZ_G4HR" imagefmt = png; 
proc sgplot data = PrepData4PZ3 noautolegend pad = (top = 1% bottom = 0% left = 1% right = 1%) sganno = label4HR;
	/* the conditional power based on original design is graphed via 3 series statements */
	series x = effHR y = condiPowerAtIA_P1 / lineattrs = (color = red thickness = 0.12cm pattern = solid);
	series x = effHR y = condiPowerAtIA_P2 / lineattrs = (color = red thickness = 0.12cm pattern = solid);
	series x = effHR y = condiPowerAtIA_P3 / lineattrs = (color = red thickness = 0.12cm pattern = dash);

	/* this is to graph the conditional power for region z1_L <= z1 <= z1_U, please note that the 
		option name = "RedS" will used in the statement of keylegend */
	series x = effHR y = AdjcondiPower /lineattrs = (color = red thickness = 0.12cm pattern = solid) name = "RedS" ;

	/* the following two statements are used to graph the sample size, also pay attention to the 
		option name = "BlackS" which will used in the statement of keylegend */
	series x = effHR y = n_SS / lineattrs = (color = black thickness = 0.12cm pattern = solid) y2axis name = "BlackS" ;
	series x = effHR y = n_SSOrig / lineattrs = (color = black thickness = 0.12cm pattern = dash) y2axis ;

	/* now define the y-axis and the grid lines */
	xaxis label = 'Hazard ratio (HR) at Interim Analysis' min = 0.5 max = 1 values = (0.5 to 1.2 by 0.1);  
	refline (0.5 to 1.2 by 0.05) / axis = x lineattrs = (pattern = dash color = gray) transparency = 0.1;

	/* now define the y-axis and the grid lines */
	yaxis min = 0 label = 'Conditional Power' values = (0 to 1 by 0.2);
	refline (0 to 1 by 0.1) / axis = y lineattrs = (pattern = dash color = gray) transparency = 0.1;

	/* this is the second y-axis */
	y2axis min = 0 label = 'Sample Size' values = (0 to 1400 by 280);
	
	/* this define the legend, list item in the order you want them */	
	keylegend "RedS" "BlackS" / border down = 2 location = inside position = topright opaque;  

	/* this draw the shaded region between z1_L and z1_U, the first shade the area, and the 2nd define the region*/
	band y = condiPowerAtIA_P1 lower = &LowerHR upper = &UpperHR / transparency = 0.7;
	refline &LowerHR &CentrHR &UpperHR / axis = x lineattrs = (thickness = 1 color = darkgray pattern = solid); 
run;

proc datasets lib = work kill memtype = data; run; quit;
