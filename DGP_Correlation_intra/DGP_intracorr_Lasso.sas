option autosignon = yes;
libname Lib "C:\Users\dutau\Downloads\SAS project\DGP_intracorr\Lasso";
rsubmit process=task1 wait=no sascmd='sas';
libname Lib "C:\Users\dutau\Downloads\SAS project\DGP_intracorr\Lasso";
						   /**********************************************************************/
						   /*********************                            *********************/
                           /*********************    DGP intracorrelation    *********************/
						   /*********************                            *********************/
						   /**********************************************************************/

						                             /** Lasso method **/


/* choose=AIC  stop=SBC */
Proc IML ;

/* The metrics */

count_perfect = 0;        
count_apparent = 0 ;
count_bad_over = 0;
count_good_over = 0;
count_under=0;

countX1=0;
countX2=0;
countX3=0;
countX4=0;
countX5=0;
countX6=0;
countX7=0;
countX8=0;
countX9=0;
countX10=0;
countX11=0;
countX12=0;
countX13=0;
countX14=0;
countX15=0;
countX16=0;
countX17=0;
countX18=0;
countX19=0;
countX20=0;
countX21=0;
countX22=0;
countX23=0;
countX24=0;
countX25=0;
countX26=0;
countX27=0;
countX28=0;
countX29=0;
countX30=0;
countX31=0;
countX32=0;
countX33=0;
countX34=0;
countX35=0;
countX36=0;
countX37=0;
countX38=0;
countX39=0;
countX40=0;
countX41=0;
countX42=0;
countX43=0;
countX44=0;
countX45=0;
countX46=0;
countX47=0;
countX48=0;
countX49=0;
countX50=0;

/* Initializing the model */

do i = 1 to 1000;

	corr = {1, 0.8, 0.65, 0.5, 0.3, 0.2};
	corr1 = toeplitz(corr);     		    /* Correlation matrix for X1-X6 */
	id  = I(94);                		    /* Correlation matrix for other variables */
	varcov = block(corr1, id);  		    /* The global variance-covariance matrix */
	mean = J(100,1,0);          		    /* The mean vector for our 100 variables */

X = RandNormal(100,mean,varcov);            /* Generating observations from Gaussian Law */
beta={1,2.1,0.6,2,1.5,2.9};                 /* Our model chosen coefficients */

Y = X[,1:6]*beta + normal(J(100,1,0))*0.1;  /* The model */

A = Y||X;
names_all = ("X1":"X100");
names_all = "Y"||names_all;                       /* Stores all varaibles names */
create Lib.A_1 from A [colname = names_all];      /* Create Database */
append from A;
close Lib.A_1;

submit;
proc glmselect data = Lib.A_1 outdesign = Lib.names_1 noprint; /* Selection model */
model y = X1-X50 / selection = Lasso(choose = AIC stop= SBC);      /*Training on the first 50 variables */
run;

proc contents data = Lib.names_1 out = Lib.reg_hat_1(keep=name) noprint; /*We collect regressors names*/
run;
endsubmit;

use Lib.reg_hat_1;
read all var {"name"} into reg_hat;
close Lib.reg_hat_1;

n = nrow(reg_hat);
reg_hat = remove(reg_hat,n);
reg ={"Intercept", "X1", "X2", "X3", "X4", "X5", "X6"};

X1 = {"X1"};
X2 = {"X2"};
X3 = {"X3"};
X4 = {"X4"};
X5 = {"X5"};
X6 = {"X6"};
X6 = {'X6'};
X7 = {'X7'};
X8 = {'X8'};
X9 = {'X9'};
X10 = {'X10'};
X11 = {'X11'};
X12 = {'X12'};
X13 = {'X13'};
X14 = {'X14'};
X15 = {'X15'};
X16 = {'X16'};
X17 = {'X17'};
X18 = {'X18'};
X19 = {'X19'};
X20 = {'X20'};
X21 = {'X21'};
X22 = {'X22'};
X23 = {'X23'};
X24 = {'X24'};
X25 = {'X25'};
X26 = {'X26'};
X27 = {'X27'};
X28 = {'X28'};
X29 = {'X29'};
X30 = {'X30'};
X31 = {'X31'};
X32 = {'X32'};
X33 = {'X33'};
X34 = {'X34'};
X35 = {'X35'};
X36 = {'X36'};
X37 = {'X37'};
X38 = {'X38'};
X39 = {'X39'};
X40 = {'X40'};
X41 = {'X41'};
X42 = {'X42'};
X43 = {'X43'};
X44 = {'X44'};
X45 = {'X45'};
X46 = {'X46'};
X47 = {'X47'};
X48 = {'X48'};
X49 = {'X49'};
X50 = {'X50'};

uni = unique(reg, reg_hat);
int = xsect(reg, reg_hat);        /* Intersection between preset variables and selected ones */

if ncol(int)= 7 & ncol(uni)=7 then count_perfect
= count_perfect+1;
if ncol(int) ^= 7 & ncol(reg_hat)=7 then
count_apparent = count_apparent+1;
if ncol(int) ^= 7 & ncol(reg_hat)>7 then
count_bad_over = count_bad_over+1;
if ncol(int) = 7 & ncol(reg_hat)>7 then
count_good_over = count_good_over+1;
if ncol(reg_hat)<7 then count_under = count_under
+1 ;

if isempty(xsect(reg_hat, X1)) then l = 1 ; else
countX1=countX1+1;
if isempty(xsect(reg_hat, X2)) then l = 1 ; else
countX2=countX2+1;
if isempty(xsect(reg_hat, X3)) then l = 1 ; else
countX3=countX3+1;
if isempty(xsect(reg_hat, X4)) then l = 1 ; else
countX4=countX4+1;
if isempty(xsect(reg_hat, X5)) then l = 1 ; else
countX5=countX5+1;
if isempty(xsect(reg_hat, X6)) then l = 1 ; else
countX6=countX6+1;
if isempty(xsect(reg_hat, X7)) then l = 1 ; else
countX7=countX7+1;
if isempty(xsect(reg_hat, X8)) then l = 1 ; else
countX8=countX8+1;
if isempty(xsect(reg_hat, X9)) then l = 1 ; else
countX9=countX9+1;
if isempty(xsect(reg_hat, X10)) then l = 1 ; else
countX10=countX10+1;
if isempty(xsect(reg_hat, X11)) then l = 1 ; else
countX11=countX11+1;
if isempty(xsect(reg_hat, X12)) then l = 1 ; else
countX12=countX12+1;
if isempty(xsect(reg_hat, X13)) then l = 1 ; else
countX13=countX13+1;
if isempty(xsect(reg_hat, X14)) then l = 1 ; else
countX14=countX14+1;
if isempty(xsect(reg_hat, X15)) then l = 1 ; else
countX15=countX15+1;
if isempty(xsect(reg_hat, X16)) then l = 1 ; else
countX16=countX16+1;
if isempty(xsect(reg_hat, X17)) then l = 1 ; else
countX17=countX17+1;
if isempty(xsect(reg_hat, X18)) then l = 1 ; else
countX18=countX18+1;
if isempty(xsect(reg_hat, X19)) then l = 1 ; else
countX19=countX19+1;
if isempty(xsect(reg_hat, X20)) then l = 1 ; else
countX20=countX20+1;
if isempty(xsect(reg_hat, X21)) then l = 1 ; else
countX21=countX21+1;
if isempty(xsect(reg_hat, X22)) then l = 1 ; else
countX22=countX22+1;
if isempty(xsect(reg_hat, X23)) then l = 1 ; else
countX23=countX23+1;
if isempty(xsect(reg_hat, X24)) then l = 1 ; else
countX24=countX24+1;
if isempty(xsect(reg_hat, X25)) then l = 1 ; else
countX25=countX25+1;
if isempty(xsect(reg_hat, X26)) then l = 1 ; else
countX26=countX26+1;
if isempty(xsect(reg_hat, X27)) then l = 1 ; else
countX27=countX27+1;
if isempty(xsect(reg_hat, X28)) then l = 1 ; else
countX28=countX28+1;
if isempty(xsect(reg_hat, X29)) then l = 1 ; else
countX29=countX29+1;
if isempty(xsect(reg_hat, X30)) then l = 1 ; else
countX30=countX30+1;
if isempty(xsect(reg_hat, X31)) then l = 1 ; else
countX31=countX31+1;
if isempty(xsect(reg_hat, X32)) then l = 1 ; else
countX32=countX32+1;
if isempty(xsect(reg_hat, X33)) then l = 1 ; else
countX33=countX33+1;
if isempty(xsect(reg_hat, X34)) then l = 1 ; else
countX34=countX34+1;
if isempty(xsect(reg_hat, X35)) then l = 1 ; else
countX35=countX35+1;
if isempty(xsect(reg_hat, X36)) then l = 1 ; else
countX36=countX36+1;
if isempty(xsect(reg_hat, X37)) then l = 1 ; else
countX37=countX37+1;
if isempty(xsect(reg_hat, X38)) then l = 1 ; else
countX38=countX38+1;
if isempty(xsect(reg_hat, X39)) then l = 1 ; else
countX39=countX39+1;
if isempty(xsect(reg_hat, X40)) then l = 1 ; else
countX40=countX40+1;
if isempty(xsect(reg_hat, X41)) then l = 1 ; else
countX41=countX41+1;
if isempty(xsect(reg_hat, X42)) then l = 1 ; else
countX42=countX42+1;
if isempty(xsect(reg_hat, X43)) then l = 1 ; else
countX43=countX43+1;
if isempty(xsect(reg_hat, X44)) then l = 1 ; else
countX44=countX44+1;
if isempty(xsect(reg_hat, X45)) then l = 1 ; else
countX45=countX45+1;
if isempty(xsect(reg_hat, X46)) then l = 1 ; else
countX46=countX46+1;
if isempty(xsect(reg_hat, X47)) then l = 1 ; else
countX47=countX47+1;
if isempty(xsect(reg_hat, X48)) then l = 1 ; else
countX48=countX48+1;
if isempty(xsect(reg_hat, X49)) then l = 1 ; else
countX49=countX49+1;
if isempty(xsect(reg_hat, X50)) then l = 1 ; else
countX50=countX50+1;
end;

/* Metrics tables */ 
pct_perfect = round((count_perfect / 1000) * 100, 0.01);
pct_apparent = round((count_apparent / 1000) * 100, 0.01);
pct_bad_over = round((count_bad_over / 1000) * 100, 0.01);
pct_good_over = round((count_good_over / 1000) * 100, 0.01);
pct_under = round((count_under / 1000) * 100, 0.01);

Final_fitting_1 = pct_perfect ||pct_bad_over 
||pct_good_over ||pct_under||pct_apparent;

Final_features_1 = countX1 ||countX2 ||countX3 ||countX4 
||countX5 ||countX6 ||countX7 ||countX8 ||countX9 
||countX10 ||countX11 ||countX12 ||countX13 ||countX14 
||countX15 ||countX16 ||countX17 ||countX18 ||countX19 
||countX20 ||countX21 ||countX22 ||countX23 ||countX24 
||countX25 ||countX26 ||countX27 ||countX28 ||countX29 
||countX30 ||countX31 ||countX32 ||countX33 ||countX34 
||countX35 ||countX36 ||countX37 ||countX38 ||countX39 
||countX40 ||countX41 ||countX42 ||countX43 ||countX44 
||countX45 ||countX46 ||countX47 ||countX48 ||countX49 
||countX50;

create lib.Final_fitting_1 from Final_fitting_1 [colname
= {"Perfect","Bad Over","Good Over", "Under", "Apparently"}];
append from Final_fitting_1;
close lib.Final_fitting_1;
Var50 = ("X1":"X50");
create lib.Final_features_1 from Final_features_1
[colname = Var50];
append from Final_features_1;
close lib.Final_features_1;
endrsubmit task1;

rsubmit process=task2 wait=no sascmd='sas';
libname Lib "C:\Users\dutau\Downloads\SAS project\DGP_intracorr\Lasso";


/* choose=SBC  stop=CP */
Proc IML ;

/* The metrics */

count_perfect = 0;        
count_apparent = 0 ;
count_bad_over = 0;
count_good_over = 0;
count_under=0;

countX1=0;
countX2=0;
countX3=0;
countX4=0;
countX5=0;
countX6=0;
countX7=0;
countX8=0;
countX9=0;
countX10=0;
countX11=0;
countX12=0;
countX13=0;
countX14=0;
countX15=0;
countX16=0;
countX17=0;
countX18=0;
countX19=0;
countX20=0;
countX21=0;
countX22=0;
countX23=0;
countX24=0;
countX25=0;
countX26=0;
countX27=0;
countX28=0;
countX29=0;
countX30=0;
countX31=0;
countX32=0;
countX33=0;
countX34=0;
countX35=0;
countX36=0;
countX37=0;
countX38=0;
countX39=0;
countX40=0;
countX41=0;
countX42=0;
countX43=0;
countX44=0;
countX45=0;
countX46=0;
countX47=0;
countX48=0;
countX49=0;
countX50=0;

/* Initializing the model */

do i = 1 to 1000;

	corr = {1, 0.8, 0.65, 0.5, 0.3, 0.2};
	corr1 = toeplitz(corr);     		    /* Correlation matrix for X1-X6 */
	id  = I(94);                		    /* Correlation matrix for other variables */
	varcov = block(corr1, id);  		    /* The global variance-covariance matrix */
	mean = J(100,1,0);          		    /* The mean vector for our 100 variables */

X = RandNormal(100,mean,varcov);            /* Generating observations from Gaussian Law */
beta={1,2.1,0.6,2,1.5,2.9};                 /* Our model chosen coefficients */

Y = X[,1:6]*beta + normal(J(100,1,0))*0.1;  /* The model */

A = Y||X;
names_all = ("X1":"X100");
names_all = "Y"||names_all;                       /* Stores all varaibles names */
create Lib.A_2 from A [colname = names_all];      /* Create Database */
append from A;
close Lib.A_2;

submit;
proc glmselect data = Lib.A_2 outdesign = Lib.names_2 noprint; /* Selection model */
model y = X1-X50 / selection = Lasso(choose = SBC stop=CP);          /*Training on the first 50 variables */
run;

proc contents data = Lib.names_2 out = Lib.reg_hat_2(keep=name) noprint; /*We collect regressors names*/
run;
endsubmit;

use Lib.reg_hat_2;
read all var {"name"} into reg_hat;
close Lib.reg_hat_2;

n = nrow(reg_hat);
reg_hat = remove(reg_hat,n);
reg ={"Intercept", "X1", "X2", "X3", "X4", "X5", "X6"};

X1 = {"X1"};
X2 = {"X2"};
X3 = {"X3"};
X4 = {"X4"};
X5 = {"X5"};
X6 = {"X6"};
X6 = {'X6'};
X7 = {'X7'};
X8 = {'X8'};
X9 = {'X9'};
X10 = {'X10'};
X11 = {'X11'};
X12 = {'X12'};
X13 = {'X13'};
X14 = {'X14'};
X15 = {'X15'};
X16 = {'X16'};
X17 = {'X17'};
X18 = {'X18'};
X19 = {'X19'};
X20 = {'X20'};
X21 = {'X21'};
X22 = {'X22'};
X23 = {'X23'};
X24 = {'X24'};
X25 = {'X25'};
X26 = {'X26'};
X27 = {'X27'};
X28 = {'X28'};
X29 = {'X29'};
X30 = {'X30'};
X31 = {'X31'};
X32 = {'X32'};
X33 = {'X33'};
X34 = {'X34'};
X35 = {'X35'};
X36 = {'X36'};
X37 = {'X37'};
X38 = {'X38'};
X39 = {'X39'};
X40 = {'X40'};
X41 = {'X41'};
X42 = {'X42'};
X43 = {'X43'};
X44 = {'X44'};
X45 = {'X45'};
X46 = {'X46'};
X47 = {'X47'};
X48 = {'X48'};
X49 = {'X49'};
X50 = {'X50'};

uni = unique(reg, reg_hat);
int = xsect(reg, reg_hat);        /* Intersection between preset variables and selected ones */

if ncol(int)= 7 & ncol(uni)=7 then count_perfect
= count_perfect+1;
if ncol(int) ^= 7 & ncol(reg_hat)=7 then
count_apparent = count_apparent+1;
if ncol(int) ^= 7 & ncol(reg_hat)>7 then
count_bad_over = count_bad_over+1;
if ncol(int) = 7 & ncol(reg_hat)>7 then
count_good_over = count_good_over+1;
if ncol(reg_hat)<7 then count_under = count_under
+1 ;

if isempty(xsect(reg_hat, X1)) then l = 1 ; else
countX1=countX1+1;
if isempty(xsect(reg_hat, X2)) then l = 1 ; else
countX2=countX2+1;
if isempty(xsect(reg_hat, X3)) then l = 1 ; else
countX3=countX3+1;
if isempty(xsect(reg_hat, X4)) then l = 1 ; else
countX4=countX4+1;
if isempty(xsect(reg_hat, X5)) then l = 1 ; else
countX5=countX5+1;
if isempty(xsect(reg_hat, X6)) then l = 1 ; else
countX6=countX6+1;
if isempty(xsect(reg_hat, X7)) then l = 1 ; else
countX7=countX7+1;
if isempty(xsect(reg_hat, X8)) then l = 1 ; else
countX8=countX8+1;
if isempty(xsect(reg_hat, X9)) then l = 1 ; else
countX9=countX9+1;
if isempty(xsect(reg_hat, X10)) then l = 1 ; else
countX10=countX10+1;
if isempty(xsect(reg_hat, X11)) then l = 1 ; else
countX11=countX11+1;
if isempty(xsect(reg_hat, X12)) then l = 1 ; else
countX12=countX12+1;
if isempty(xsect(reg_hat, X13)) then l = 1 ; else
countX13=countX13+1;
if isempty(xsect(reg_hat, X14)) then l = 1 ; else
countX14=countX14+1;
if isempty(xsect(reg_hat, X15)) then l = 1 ; else
countX15=countX15+1;
if isempty(xsect(reg_hat, X16)) then l = 1 ; else
countX16=countX16+1;
if isempty(xsect(reg_hat, X17)) then l = 1 ; else
countX17=countX17+1;
if isempty(xsect(reg_hat, X18)) then l = 1 ; else
countX18=countX18+1;
if isempty(xsect(reg_hat, X19)) then l = 1 ; else
countX19=countX19+1;
if isempty(xsect(reg_hat, X20)) then l = 1 ; else
countX20=countX20+1;
if isempty(xsect(reg_hat, X21)) then l = 1 ; else
countX21=countX21+1;
if isempty(xsect(reg_hat, X22)) then l = 1 ; else
countX22=countX22+1;
if isempty(xsect(reg_hat, X23)) then l = 1 ; else
countX23=countX23+1;
if isempty(xsect(reg_hat, X24)) then l = 1 ; else
countX24=countX24+1;
if isempty(xsect(reg_hat, X25)) then l = 1 ; else
countX25=countX25+1;
if isempty(xsect(reg_hat, X26)) then l = 1 ; else
countX26=countX26+1;
if isempty(xsect(reg_hat, X27)) then l = 1 ; else
countX27=countX27+1;
if isempty(xsect(reg_hat, X28)) then l = 1 ; else
countX28=countX28+1;
if isempty(xsect(reg_hat, X29)) then l = 1 ; else
countX29=countX29+1;
if isempty(xsect(reg_hat, X30)) then l = 1 ; else
countX30=countX30+1;
if isempty(xsect(reg_hat, X31)) then l = 1 ; else
countX31=countX31+1;
if isempty(xsect(reg_hat, X32)) then l = 1 ; else
countX32=countX32+1;
if isempty(xsect(reg_hat, X33)) then l = 1 ; else
countX33=countX33+1;
if isempty(xsect(reg_hat, X34)) then l = 1 ; else
countX34=countX34+1;
if isempty(xsect(reg_hat, X35)) then l = 1 ; else
countX35=countX35+1;
if isempty(xsect(reg_hat, X36)) then l = 1 ; else
countX36=countX36+1;
if isempty(xsect(reg_hat, X37)) then l = 1 ; else
countX37=countX37+1;
if isempty(xsect(reg_hat, X38)) then l = 1 ; else
countX38=countX38+1;
if isempty(xsect(reg_hat, X39)) then l = 1 ; else
countX39=countX39+1;
if isempty(xsect(reg_hat, X40)) then l = 1 ; else
countX40=countX40+1;
if isempty(xsect(reg_hat, X41)) then l = 1 ; else
countX41=countX41+1;
if isempty(xsect(reg_hat, X42)) then l = 1 ; else
countX42=countX42+1;
if isempty(xsect(reg_hat, X43)) then l = 1 ; else
countX43=countX43+1;
if isempty(xsect(reg_hat, X44)) then l = 1 ; else
countX44=countX44+1;
if isempty(xsect(reg_hat, X45)) then l = 1 ; else
countX45=countX45+1;
if isempty(xsect(reg_hat, X46)) then l = 1 ; else
countX46=countX46+1;
if isempty(xsect(reg_hat, X47)) then l = 1 ; else
countX47=countX47+1;
if isempty(xsect(reg_hat, X48)) then l = 1 ; else
countX48=countX48+1;
if isempty(xsect(reg_hat, X49)) then l = 1 ; else
countX49=countX49+1;
if isempty(xsect(reg_hat, X50)) then l = 1 ; else
countX50=countX50+1;
end;

/* Metrics tables */ 
pct_perfect = round((count_perfect / 1000) * 100, 0.01);
pct_apparent = round((count_apparent / 1000) * 100, 0.01);
pct_bad_over = round((count_bad_over / 1000) * 100, 0.01);
pct_good_over = round((count_good_over / 1000) * 100, 0.01);
pct_under = round((count_under / 1000) * 100, 0.01);

Final_fitting_2 = pct_perfect ||pct_bad_over 
||pct_good_over ||pct_under||pct_apparent;

Final_features_2 = countX1 ||countX2 ||countX3 ||countX4 
||countX5 ||countX6 ||countX7 ||countX8 ||countX9 
||countX10 ||countX11 ||countX12 ||countX13 ||countX14 
||countX15 ||countX16 ||countX17 ||countX18 ||countX19 
||countX20 ||countX21 ||countX22 ||countX23 ||countX24 
||countX25 ||countX26 ||countX27 ||countX28 ||countX29 
||countX30 ||countX31 ||countX32 ||countX33 ||countX34 
||countX35 ||countX36 ||countX37 ||countX38 ||countX39 
||countX40 ||countX41 ||countX42 ||countX43 ||countX44 
||countX45 ||countX46 ||countX47 ||countX48 ||countX49 
||countX50;

create lib.Final_fitting_2 from Final_fitting_2 [colname
= {"Perfect","Bad Over","Good Over", "Under",
"Apparently"}];
append from Final_fitting_2;
close lib.Final_fitting_2;
Var50 = ("X1":"X50");
create lib.Final_features_2 from Final_features_2
[colname = Var50];
append from Final_features_2;
close lib.Final_features_2;
endrsubmit task2;

rsubmit process=task3 wait=no sascmd='sas';
libname Lib "C:\Users\dutau\Downloads\SAS project\DGP_intracorr\Lasso";

/* choose=CV  stop=SBC */

Proc IML ;

/* The metrics */

count_perfect = 0;        
count_apparent = 0 ;
count_bad_over = 0;
count_good_over = 0;
count_under=0;

countX1=0;
countX2=0;
countX3=0;
countX4=0;
countX5=0;
countX6=0;
countX7=0;
countX8=0;
countX9=0;
countX10=0;
countX11=0;
countX12=0;
countX13=0;
countX14=0;
countX15=0;
countX16=0;
countX17=0;
countX18=0;
countX19=0;
countX20=0;
countX21=0;
countX22=0;
countX23=0;
countX24=0;
countX25=0;
countX26=0;
countX27=0;
countX28=0;
countX29=0;
countX30=0;
countX31=0;
countX32=0;
countX33=0;
countX34=0;
countX35=0;
countX36=0;
countX37=0;
countX38=0;
countX39=0;
countX40=0;
countX41=0;
countX42=0;
countX43=0;
countX44=0;
countX45=0;
countX46=0;
countX47=0;
countX48=0;
countX49=0;
countX50=0;

/* Initializing the model */

do i = 1 to 1000;

	corr = {1, 0.8, 0.65, 0.5, 0.3, 0.2};
	corr1 = toeplitz(corr);     		    /* Correlation matrix for X1-X6 */
	id  = I(94);                		    /* Correlation matrix for other variables */
	varcov = block(corr1, id);  		    /* The global variance-covariance matrix */
	mean = J(100,1,0);          		    /* The mean vector for our 100 variables */

X = RandNormal(100,mean,varcov);            /* Generating observations from Gaussian Law */
beta={1,2.1,0.6,2,1.5,2.9};                 /* Our model chosen coefficients */

Y = X[,1:6]*beta + normal(J(100,1,0))*0.1;  /* The model */

A = Y||X;
names_all = ("X1":"X100");
names_all = "Y"||names_all;                       /* Stores all varaibles names */
create Lib.A_3 from A [colname = names_all];      /* Create Database */
append from A;
close Lib.A_3;

submit;
proc glmselect data = Lib.A_3 outdesign = Lib.names_3 noprint; /* Selection model */
model y = X1-X50 / selection = Lasso(choose = CV stop=SBC);      /*Training on the first 50 variables */
run;

proc contents data = Lib.names_3 out = Lib.reg_hat_3(keep=name) noprint; /*We collect regressors names*/
run;
endsubmit;

use Lib.reg_hat_3;
read all var {"name"} into reg_hat;
close Lib.reg_hat_3;

n = nrow(reg_hat);
reg_hat = remove(reg_hat,n);
reg ={"Intercept", "X1", "X2", "X3", "X4", "X5", "X6"};

X1 = {"X1"};
X2 = {"X2"};
X3 = {"X3"};
X4 = {"X4"};
X5 = {"X5"};
X6 = {"X6"};
X6 = {'X6'};
X7 = {'X7'};
X8 = {'X8'};
X9 = {'X9'};
X10 = {'X10'};
X11 = {'X11'};
X12 = {'X12'};
X13 = {'X13'};
X14 = {'X14'};
X15 = {'X15'};
X16 = {'X16'};
X17 = {'X17'};
X18 = {'X18'};
X19 = {'X19'};
X20 = {'X20'};
X21 = {'X21'};
X22 = {'X22'};
X23 = {'X23'};
X24 = {'X24'};
X25 = {'X25'};
X26 = {'X26'};
X27 = {'X27'};
X28 = {'X28'};
X29 = {'X29'};
X30 = {'X30'};
X31 = {'X31'};
X32 = {'X32'};
X33 = {'X33'};
X34 = {'X34'};
X35 = {'X35'};
X36 = {'X36'};
X37 = {'X37'};
X38 = {'X38'};
X39 = {'X39'};
X40 = {'X40'};
X41 = {'X41'};
X42 = {'X42'};
X43 = {'X43'};
X44 = {'X44'};
X45 = {'X45'};
X46 = {'X46'};
X47 = {'X47'};
X48 = {'X48'};
X49 = {'X49'};
X50 = {'X50'};

uni = unique(reg, reg_hat);
int = xsect(reg, reg_hat);        /* Intersection between preset variables and selected ones */

if ncol(int)= 7 & ncol(uni)=7 then count_perfect
= count_perfect+1;
if ncol(int) ^= 7 & ncol(reg_hat)=7 then
count_apparent = count_apparent+1;
if ncol(int) ^= 7 & ncol(reg_hat)>7 then
count_bad_over = count_bad_over+1;
if ncol(int) = 7 & ncol(reg_hat)>7 then
count_good_over = count_good_over+1;
if ncol(reg_hat)<7 then count_under = count_under
+1 ;

if isempty(xsect(reg_hat, X1)) then l = 1 ; else
countX1=countX1+1;
if isempty(xsect(reg_hat, X2)) then l = 1 ; else
countX2=countX2+1;
if isempty(xsect(reg_hat, X3)) then l = 1 ; else
countX3=countX3+1;
if isempty(xsect(reg_hat, X4)) then l = 1 ; else
countX4=countX4+1;
if isempty(xsect(reg_hat, X5)) then l = 1 ; else
countX5=countX5+1;
if isempty(xsect(reg_hat, X6)) then l = 1 ; else
countX6=countX6+1;
if isempty(xsect(reg_hat, X7)) then l = 1 ; else
countX7=countX7+1;
if isempty(xsect(reg_hat, X8)) then l = 1 ; else
countX8=countX8+1;
if isempty(xsect(reg_hat, X9)) then l = 1 ; else
countX9=countX9+1;
if isempty(xsect(reg_hat, X10)) then l = 1 ; else
countX10=countX10+1;
if isempty(xsect(reg_hat, X11)) then l = 1 ; else
countX11=countX11+1;
if isempty(xsect(reg_hat, X12)) then l = 1 ; else
countX12=countX12+1;
if isempty(xsect(reg_hat, X13)) then l = 1 ; else
countX13=countX13+1;
if isempty(xsect(reg_hat, X14)) then l = 1 ; else
countX14=countX14+1;
if isempty(xsect(reg_hat, X15)) then l = 1 ; else
countX15=countX15+1;
if isempty(xsect(reg_hat, X16)) then l = 1 ; else
countX16=countX16+1;
if isempty(xsect(reg_hat, X17)) then l = 1 ; else
countX17=countX17+1;
if isempty(xsect(reg_hat, X18)) then l = 1 ; else
countX18=countX18+1;
if isempty(xsect(reg_hat, X19)) then l = 1 ; else
countX19=countX19+1;
if isempty(xsect(reg_hat, X20)) then l = 1 ; else
countX20=countX20+1;
if isempty(xsect(reg_hat, X21)) then l = 1 ; else
countX21=countX21+1;
if isempty(xsect(reg_hat, X22)) then l = 1 ; else
countX22=countX22+1;
if isempty(xsect(reg_hat, X23)) then l = 1 ; else
countX23=countX23+1;
if isempty(xsect(reg_hat, X24)) then l = 1 ; else
countX24=countX24+1;
if isempty(xsect(reg_hat, X25)) then l = 1 ; else
countX25=countX25+1;
if isempty(xsect(reg_hat, X26)) then l = 1 ; else
countX26=countX26+1;
if isempty(xsect(reg_hat, X27)) then l = 1 ; else
countX27=countX27+1;
if isempty(xsect(reg_hat, X28)) then l = 1 ; else
countX28=countX28+1;
if isempty(xsect(reg_hat, X29)) then l = 1 ; else
countX29=countX29+1;
if isempty(xsect(reg_hat, X30)) then l = 1 ; else
countX30=countX30+1;
if isempty(xsect(reg_hat, X31)) then l = 1 ; else
countX31=countX31+1;
if isempty(xsect(reg_hat, X32)) then l = 1 ; else
countX32=countX32+1;
if isempty(xsect(reg_hat, X33)) then l = 1 ; else
countX33=countX33+1;
if isempty(xsect(reg_hat, X34)) then l = 1 ; else
countX34=countX34+1;
if isempty(xsect(reg_hat, X35)) then l = 1 ; else
countX35=countX35+1;
if isempty(xsect(reg_hat, X36)) then l = 1 ; else
countX36=countX36+1;
if isempty(xsect(reg_hat, X37)) then l = 1 ; else
countX37=countX37+1;
if isempty(xsect(reg_hat, X38)) then l = 1 ; else
countX38=countX38+1;
if isempty(xsect(reg_hat, X39)) then l = 1 ; else
countX39=countX39+1;
if isempty(xsect(reg_hat, X40)) then l = 1 ; else
countX40=countX40+1;
if isempty(xsect(reg_hat, X41)) then l = 1 ; else
countX41=countX41+1;
if isempty(xsect(reg_hat, X42)) then l = 1 ; else
countX42=countX42+1;
if isempty(xsect(reg_hat, X43)) then l = 1 ; else
countX43=countX43+1;
if isempty(xsect(reg_hat, X44)) then l = 1 ; else
countX44=countX44+1;
if isempty(xsect(reg_hat, X45)) then l = 1 ; else
countX45=countX45+1;
if isempty(xsect(reg_hat, X46)) then l = 1 ; else
countX46=countX46+1;
if isempty(xsect(reg_hat, X47)) then l = 1 ; else
countX47=countX47+1;
if isempty(xsect(reg_hat, X48)) then l = 1 ; else
countX48=countX48+1;
if isempty(xsect(reg_hat, X49)) then l = 1 ; else
countX49=countX49+1;
if isempty(xsect(reg_hat, X50)) then l = 1 ; else
countX50=countX50+1;
end;

/* Metrics tables */ 
pct_perfect = round((count_perfect / 1000) * 100, 0.01);
pct_apparent = round((count_apparent / 1000) * 100, 0.01);
pct_bad_over = round((count_bad_over / 1000) * 100, 0.01);
pct_good_over = round((count_good_over / 1000) * 100, 0.01);
pct_under = round((count_under / 1000) * 100, 0.01);

Final_fitting_3 = pct_perfect ||pct_bad_over 
||pct_good_over ||pct_under||pct_apparent;

Final_features_3 = countX1 ||countX2 ||countX3 ||countX4 
||countX5 ||countX6 ||countX7 ||countX8 ||countX9 
||countX10 ||countX11 ||countX12 ||countX13 ||countX14 
||countX15 ||countX16 ||countX17 ||countX18 ||countX19 
||countX20 ||countX21 ||countX22 ||countX23 ||countX24 
||countX25 ||countX26 ||countX27 ||countX28 ||countX29 
||countX30 ||countX31 ||countX32 ||countX33 ||countX34 
||countX35 ||countX36 ||countX37 ||countX38 ||countX39 
||countX40 ||countX41 ||countX42 ||countX43 ||countX44 
||countX45 ||countX46 ||countX47 ||countX48 ||countX49 
||countX50;

create lib.Final_fitting_3 from Final_fitting_3 [colname
= {"Perfect","Bad Over","Good Over", "Under",
"Apparently"}];
append from Final_fitting_3;
close lib.Final_fitting_3;
Var50 = ("X1":"X50");
create lib.Final_features_3 from Final_features_3
[colname = Var50];
append from Final_features_3;
close lib.Final_features_3;
endrsubmit task3;

rsubmit process=task4 wait=no sascmd='sas';
libname Lib "C:\Users\dutau\Downloads\SAS project\DGP_intracorr\Lasso";


/* choose=PRESS  stop=SBC */
Proc IML ;

/* The metrics */

count_perfect = 0;        
count_apparent = 0 ;
count_bad_over = 0;
count_good_over = 0;
count_under=0;

countX1=0;
countX2=0;
countX3=0;
countX4=0;
countX5=0;
countX6=0;
countX7=0;
countX8=0;
countX9=0;
countX10=0;
countX11=0;
countX12=0;
countX13=0;
countX14=0;
countX15=0;
countX16=0;
countX17=0;
countX18=0;
countX19=0;
countX20=0;
countX21=0;
countX22=0;
countX23=0;
countX24=0;
countX25=0;
countX26=0;
countX27=0;
countX28=0;
countX29=0;
countX30=0;
countX31=0;
countX32=0;
countX33=0;
countX34=0;
countX35=0;
countX36=0;
countX37=0;
countX38=0;
countX39=0;
countX40=0;
countX41=0;
countX42=0;
countX43=0;
countX44=0;
countX45=0;
countX46=0;
countX47=0;
countX48=0;
countX49=0;
countX50=0;

/* Initializing the model */

do i = 1 to 1000;

	corr = {1, 0.8, 0.65, 0.5, 0.3, 0.2};
	corr1 = toeplitz(corr);     		    /* Correlation matrix for X1-X6 */
	id  = I(94);                		    /* Correlation matrix for other variables */
	varcov = block(corr1, id);  		    /* The global variance-covariance matrix */
	mean = J(100,1,0);          		    /* The mean vector for our 100 variables */

X = RandNormal(100,mean,varcov);            /* Generating observations from Gaussian Law */
beta={1,2.1,0.6,2,1.5,2.9};                 /* Our model chosen coefficients */

Y = X[,1:6]*beta + normal(J(100,1,0))*0.1;  /* The model */

A = Y||X;
names_all = ("X1":"X100");
names_all = "Y"||names_all;                       /* Stores all varaibles names */
create Lib.A_4 from A [colname = names_all];      /* Create Database */
append from A;
close Lib.A_4;

submit;
proc glmselect data = Lib.A_4 outdesign = Lib.names_4 noprint; /* Selection model */
model y = X1-X50 / selection = Lasso(lscoeffs choose = PRESS stop=SBC);      /*Training on the first 50 variables */
run;

proc contents data = Lib.names_4 out = Lib.reg_hat_4(keep=name) noprint; /*We collect regressors names*/
run;
endsubmit;

use Lib.reg_hat_4;
read all var {"name"} into reg_hat;
close Lib.reg_hat_4;

n = nrow(reg_hat);
reg_hat = remove(reg_hat,n);
reg ={"Intercept", "X1", "X2", "X3", "X4", "X5", "X6"};

X1 = {"X1"};
X2 = {"X2"};
X3 = {"X3"};
X4 = {"X4"};
X5 = {"X5"};
X6 = {"X6"};
X6 = {'X6'};
X7 = {'X7'};
X8 = {'X8'};
X9 = {'X9'};
X10 = {'X10'};
X11 = {'X11'};
X12 = {'X12'};
X13 = {'X13'};
X14 = {'X14'};
X15 = {'X15'};
X16 = {'X16'};
X17 = {'X17'};
X18 = {'X18'};
X19 = {'X19'};
X20 = {'X20'};
X21 = {'X21'};
X22 = {'X22'};
X23 = {'X23'};
X24 = {'X24'};
X25 = {'X25'};
X26 = {'X26'};
X27 = {'X27'};
X28 = {'X28'};
X29 = {'X29'};
X30 = {'X30'};
X31 = {'X31'};
X32 = {'X32'};
X33 = {'X33'};
X34 = {'X34'};
X35 = {'X35'};
X36 = {'X36'};
X37 = {'X37'};
X38 = {'X38'};
X39 = {'X39'};
X40 = {'X40'};
X41 = {'X41'};
X42 = {'X42'};
X43 = {'X43'};
X44 = {'X44'};
X45 = {'X45'};
X46 = {'X46'};
X47 = {'X47'};
X48 = {'X48'};
X49 = {'X49'};
X50 = {'X50'};

uni = unique(reg, reg_hat);
int = xsect(reg, reg_hat);        /* Intersection between preset variables and selected ones */

if ncol(int)= 7 & ncol(uni)=7 then count_perfect
= count_perfect+1;
if ncol(int) ^= 7 & ncol(reg_hat)=7 then
count_apparent = count_apparent+1;
if ncol(int) ^= 7 & ncol(reg_hat)>7 then
count_bad_over = count_bad_over+1;
if ncol(int) = 7 & ncol(reg_hat)>7 then
count_good_over = count_good_over+1;
if ncol(reg_hat)<7 then count_under = count_under
+1 ;

if isempty(xsect(reg_hat, X1)) then l = 1 ; else
countX1=countX1+1;
if isempty(xsect(reg_hat, X2)) then l = 1 ; else
countX2=countX2+1;
if isempty(xsect(reg_hat, X3)) then l = 1 ; else
countX3=countX3+1;
if isempty(xsect(reg_hat, X4)) then l = 1 ; else
countX4=countX4+1;
if isempty(xsect(reg_hat, X5)) then l = 1 ; else
countX5=countX5+1;
if isempty(xsect(reg_hat, X6)) then l = 1 ; else
countX6=countX6+1;
if isempty(xsect(reg_hat, X7)) then l = 1 ; else
countX7=countX7+1;
if isempty(xsect(reg_hat, X8)) then l = 1 ; else
countX8=countX8+1;
if isempty(xsect(reg_hat, X9)) then l = 1 ; else
countX9=countX9+1;
if isempty(xsect(reg_hat, X10)) then l = 1 ; else
countX10=countX10+1;
if isempty(xsect(reg_hat, X11)) then l = 1 ; else
countX11=countX11+1;
if isempty(xsect(reg_hat, X12)) then l = 1 ; else
countX12=countX12+1;
if isempty(xsect(reg_hat, X13)) then l = 1 ; else
countX13=countX13+1;
if isempty(xsect(reg_hat, X14)) then l = 1 ; else
countX14=countX14+1;
if isempty(xsect(reg_hat, X15)) then l = 1 ; else
countX15=countX15+1;
if isempty(xsect(reg_hat, X16)) then l = 1 ; else
countX16=countX16+1;
if isempty(xsect(reg_hat, X17)) then l = 1 ; else
countX17=countX17+1;
if isempty(xsect(reg_hat, X18)) then l = 1 ; else
countX18=countX18+1;
if isempty(xsect(reg_hat, X19)) then l = 1 ; else
countX19=countX19+1;
if isempty(xsect(reg_hat, X20)) then l = 1 ; else
countX20=countX20+1;
if isempty(xsect(reg_hat, X21)) then l = 1 ; else
countX21=countX21+1;
if isempty(xsect(reg_hat, X22)) then l = 1 ; else
countX22=countX22+1;
if isempty(xsect(reg_hat, X23)) then l = 1 ; else
countX23=countX23+1;
if isempty(xsect(reg_hat, X24)) then l = 1 ; else
countX24=countX24+1;
if isempty(xsect(reg_hat, X25)) then l = 1 ; else
countX25=countX25+1;
if isempty(xsect(reg_hat, X26)) then l = 1 ; else
countX26=countX26+1;
if isempty(xsect(reg_hat, X27)) then l = 1 ; else
countX27=countX27+1;
if isempty(xsect(reg_hat, X28)) then l = 1 ; else
countX28=countX28+1;
if isempty(xsect(reg_hat, X29)) then l = 1 ; else
countX29=countX29+1;
if isempty(xsect(reg_hat, X30)) then l = 1 ; else
countX30=countX30+1;
if isempty(xsect(reg_hat, X31)) then l = 1 ; else
countX31=countX31+1;
if isempty(xsect(reg_hat, X32)) then l = 1 ; else
countX32=countX32+1;
if isempty(xsect(reg_hat, X33)) then l = 1 ; else
countX33=countX33+1;
if isempty(xsect(reg_hat, X34)) then l = 1 ; else
countX34=countX34+1;
if isempty(xsect(reg_hat, X35)) then l = 1 ; else
countX35=countX35+1;
if isempty(xsect(reg_hat, X36)) then l = 1 ; else
countX36=countX36+1;
if isempty(xsect(reg_hat, X37)) then l = 1 ; else
countX37=countX37+1;
if isempty(xsect(reg_hat, X38)) then l = 1 ; else
countX38=countX38+1;
if isempty(xsect(reg_hat, X39)) then l = 1 ; else
countX39=countX39+1;
if isempty(xsect(reg_hat, X40)) then l = 1 ; else
countX40=countX40+1;
if isempty(xsect(reg_hat, X41)) then l = 1 ; else
countX41=countX41+1;
if isempty(xsect(reg_hat, X42)) then l = 1 ; else
countX42=countX42+1;
if isempty(xsect(reg_hat, X43)) then l = 1 ; else
countX43=countX43+1;
if isempty(xsect(reg_hat, X44)) then l = 1 ; else
countX44=countX44+1;
if isempty(xsect(reg_hat, X45)) then l = 1 ; else
countX45=countX45+1;
if isempty(xsect(reg_hat, X46)) then l = 1 ; else
countX46=countX46+1;
if isempty(xsect(reg_hat, X47)) then l = 1 ; else
countX47=countX47+1;
if isempty(xsect(reg_hat, X48)) then l = 1 ; else
countX48=countX48+1;
if isempty(xsect(reg_hat, X49)) then l = 1 ; else
countX49=countX49+1;
if isempty(xsect(reg_hat, X50)) then l = 1 ; else
countX50=countX50+1;
end;

/* Metrics tables */ 
pct_perfect = round((count_perfect / 1000) * 100, 0.01);
pct_apparent = round((count_apparent / 1000) * 100, 0.01);
pct_bad_over = round((count_bad_over / 1000) * 100, 0.01);
pct_good_over = round((count_good_over / 1000) * 100, 0.01);
pct_under = round((count_under / 1000) * 100, 0.01);

Final_fitting_4 = pct_perfect ||pct_bad_over 
||pct_good_over ||pct_under||pct_apparent;

Final_features_4 = countX1 ||countX2 ||countX3 ||countX4 
||countX5 ||countX6 ||countX7 ||countX8 ||countX9 
||countX10 ||countX11 ||countX12 ||countX13 ||countX14 
||countX15 ||countX16 ||countX17 ||countX18 ||countX19 
||countX20 ||countX21 ||countX22 ||countX23 ||countX24 
||countX25 ||countX26 ||countX27 ||countX28 ||countX29 
||countX30 ||countX31 ||countX32 ||countX33 ||countX34 
||countX35 ||countX36 ||countX37 ||countX38 ||countX39 
||countX40 ||countX41 ||countX42 ||countX43 ||countX44 
||countX45 ||countX46 ||countX47 ||countX48 ||countX49 
||countX50;

create lib.Final_fitting_4 from Final_fitting_4 [colname
= {"Perfect","Bad Over","Good Over", "Under", "Apparently"}];
append from Final_fitting_4;
close lib.Final_fitting_4;
Var50 = ("X1":"X50");
create lib.Final_features_4 from Final_features_4
[colname = Var50];
append from Final_features_4;
close lib.Final_features_4;
endrsubmit task4;
waitfor _all_ task1 task2 task3 task4;

data Fitting ; set lib.Final_fitting_1
lib.Final_fitting_2
lib.Final_fitting_3
lib.Final_fitting_4 ; 
run ;
data Features ; set lib.Final_features_1
lib.Final_features_2
lib.Final_features_3
lib.Final_features_4 ;
run ;
proc print data = Fitting ; run ;
proc print data = Features ; run ;
signoff task1;
signoff task2;
signoff task3;
signoff task4;


