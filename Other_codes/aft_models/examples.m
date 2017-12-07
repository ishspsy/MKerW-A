
close all
clear all
clc


data=[1      1  0.6  77     76     1
2      1  1.3  53     71     1
3      1  2.4  45     71     1
4      1  2.5  57     78     0
5      1  3.2  58     74     1
6      1  3.2  51     77     0
7      1  3.3  76     74     1
8      1  3.3  63     77     0
9      1  3.5  43     71     1
10     1  3.5  60     73     1
11     1  4.0  52     71     1
12     1  4.0  63     76     1
13     1  4.3  86     74     1
14     1  4.5  48     76     0
15     1  4.5  68     76     0
16     1  5.3  81     72     1
17     1  5.5  70     75     0
18     1  5.9  58     75     0
19     1  5.9  47     75     0
20     1  6.0  75     73     1
21     1  6.1  77     75     0
22     1  6.2  64     75     0
23     1  6.4  77     72     1
24     1  6.5  67     70     1
25     1  6.5  79     74     0
26     1  6.7  61     74     0
27     1  7.0  66     74     0
28     1  7.4  68     71     1
29     1  7.4  73     73     0
30     1  8.1  56     73     0
31     1  8.1  73     73     0
32     1  9.6  58     71     0
33     1 10.7  68     70     0
34     2  0.2  86     74     1
35     2  1.8  64     77     1
36     2  2.0  63     75     1
37     2  2.2  71     78     0
38     2  2.6  67     78     0
39     2  3.3  51     77     0
40     2  3.6  70     77     1
41     2  3.6  72     77     0
42     2  4.0  81     71     1
43     2  4.3  47     76     0
44     2  4.3  64     76     0
45     2  5.0  66     76     0
46     2  6.2  74     72     1
47     2  7.0  62     73     1
48     2  7.5  50     73     0
49     2  7.6  53     73     0
50     2  9.3  61     71     0
51     3  0.3  49     72     1
52     3  0.3  71     76     1
53     3  0.5  57     74     1
54     3  0.7  79     77     1
55     3  0.8  82     74     1
56     3  1.0  49     76     1
57     3  1.3  60     76     1
58     3  1.6  64     72     1
59     3  1.8  74     71     1
60     3  1.9  72     74     1
61     3  1.9  53     74     1
62     3  3.2  54     75     1
63     3  3.5  81     74     1
64     3  3.7  52     77     0
65     3  4.5  66     76     0
66     3  4.8  54     76     0
67     3  4.8  63     76     0
68     3  5.0  59     73     1
69     3  5.0  49     76     0
70     3  5.1  69     76     0
71     3  6.3  70     72     1
72     3  6.4  65     72     1
73     3  6.5  65     74     0
74     3  7.8  68     72     1
75     3  8.0  78     73     0
76     3  9.3  69     71     0
77     3 10.1  51     71     0
78     4  0.1  65     72     1
79     4  0.3  71     76     1
80     4  0.4  76     77     1
81     4  0.8  65     76     1
82     4  0.8  78     77     1
83     4  1.0  41     77     1
84     4  1.5  68     73     1
85     4  2.0  69     76     1
86     4  2.3  62     71     1
87     4  2.9  74     78     0
88     4  3.6  71     75     1
89     4  3.8  84     74     1
90     4  4.3  48     76     0];


%Larynx cancer data that are available in Klein and Moeschberger (2003)
%"Survival Analysis: Techniques for Censored and Truncated Data", Springer. 
%Variable "time" is the time to event
%Variable "age" is the age of each patient
%Variable "stage" is the stage of the cancer (can take the values 1 or or 3 or 4
%depending on the progression of the disease)
%Variable status is the status/censoring indicator. 1 for death, 0 for
%censoring. For a more detailed discussion of the data see Klein and
%Moeschberger (2003) and Kardaun Stat. Nederlandica 37 (1983), 103-126.

%The data can also be found here:
%http://www.mcw.edu/FileLibrary/Groups/Biostatistics/Publicfiles/DataFromSection/DataFromSectionTXT/Data_from_section_1.8.txt

stage=data(:,2);
age=data(:,4);
status=data(:,end);
time=data(:,3);
n=length(time);
%We ignore the other variables for this illustration.

%Note that the status variable in this data set has the coding:
%event=1, right censored=0.

%This routine is consistent with MATLAB coding for right censoring that is
%event=0, right censored=1

%So we need to transform status:
status=-status+1;

%If we also had left censored data we should somehow
%try to build status having the following coding:
%left censored=-1, event=0, right censored=1



%%

%Now we fit an AFT model with an intercept as well as one covariate, the
%age, assuming the weibull distribution. First build the covariate matrix Z:
Z=[ones(length(time),1) age];

%Then use the "aft" function:
%(Note that the "defineevent" argument is set to 1 since 1 is for death).
[pars covpars SE CI Zscores pvalues gval exitflag]=aft(time,Z,status,'weibull');
pars
covpars
SE
%The estimated regression parameters (with the scale parameter being 
%the last element of the vector pars) are given in output vector "pars".
%The corresponding covariance matrix and the standard errors are given in
%covpars and SE respectively. 

%%
clc
%Suppose we want to see the summary of the previous model in a table: 
[pars covpars SE CI Zscores pvalues gval exitflag]=aft(time,Z,status,'weibull',[],[],[],1);

%%

%Lets see how we can accomodate the stage covariate which is not continuous as the age.
%We have four  stages of the disease. We may consider stage 1 as the
%baseline stage. Hence the design matrix Z will now be of the form:
Z=[ones(length(time),1) age (stage==2) (stage==3) (stage==4)];
%Z=[ age (stage==1) (stage==2) (stage==3) ];

%and reperform analysis by using again the weibull distribution and also ask
%for a results table:
[pars covpars SE CI Zscores pvalues gval exitflag]=aft(time,Z,status,'weibull',[],[],[],1);

%%
%Now lets see how one can use his own initial values, say all values equal to 1. 
%When the Weibull is employed, and for Z as built in the previous cell the number of 
%parameters is  6: g0,g1,g2,g3,g4,scale, thus one can define
%clc

init=[3 0 -0.5 -0.5 -2 1];
[pars covpars SE CI Zscores pvalues gval exitflag]=aft(time,Z,status,'weibull',init,1,[],1);

%and we get the same results as before.

%%

%Now change the minimizer and use the fmincon with the sqp algorithm, (that is 
%set the minimizer equal to 2)
[pars covpars SE CI Zscores pvalues gval exitflag]=aft(time,Z,status,'weibull',[],2,[],1);

%we get the same results as before.

%%

%Set some options for example for the fminsearchbnd (that is minimizer=1)
options=optimset('MaxFunEvals', 10000, 'MaxIter',10000);
[pars covpars SE CI Zscores pvalues gval exitflag]=aft(time,Z,status,'weibull',[],1,options,1);

%%

%Use the broader model of the generalized model
[pars covpars SE CI Zscores pvalues gval exitflag]=aft(time,Z,status,'gamma',[],1,[],1);
%Observe that there is an extra shape parameter here. The generalized gama
%model is known for its computational problems. If convergence problems
%occur the first thing I recommend is to try all minimizers that is 1,2 or
%3. If this does not work you might want to rerun with different initial values.
%A nice idea would be to first fit the weibull or the gamma aft model and use their
%estimated parameters as initial values for the generalized gamma. (The gamma 
%and the weibull are special cases of the generalized gamma).


%The parameterization (lets say of the Weibull or the Generalized Gamma) is appear to be
%different from text to text but equivalent in many cases. The parameterization used here
%is the same used by the SAS in "proc lifereg". See also the
%documentation in the supported distributions section of "proc lifereg"
%to check on the parameterizations. A different but equivalent parameterization of the Generalized Gamma
%distribution can be found in "Statistical Methods for Survival Analysis"
%Elisa T. Lee, John Wenyu Wang - 2003, pages 277 and 278. In particular
%equation 11.6.7 shows the way one can move from the parameterization used
%here to the parameterization used in that book (and vice versa).
%See also Klein and Moeschberger (2003) for a general review of AFT models.



