addpath(genpath(pwd))


%% Figures 4-5

%pvalue
load('pheat.csv')
yvalue={'KIRP','SARC','PRAD','BRCA','KIRC','MESO','THCA','CESC','UVM','READ','BLCA','COAD','LIHC','UCEC','LUSC','ESCA','STAD','OV','PAAD','HNSC','LUAD','LGG'};
xvalue={'MKerW-A','MKer-A','K-R','Cons-A','Cons-R','iCluster-A','S-R','SIM-R','P-A','C-A','SIM-A','K-M','S-A','Cons-M','Ker-A','SIM-M','K-A','K-C','Cons-C','SIM-C','S-M','S-C'};
h=heatmap(xvalue,yvalue,pheat);
h.Colormap=flipud(hot);
h.ColorScaling = 'scaledrows';
h.GridVisible='off';
h.Title='P-value';
h.FontSize = 12;
h.ColorLimits=[0 0.5];

%min
load('minheat.csv')
yvalue={'KIRP','SARC','PRAD','BRCA','KIRC','MESO','THCA','CESC','UVM','READ','BLCA','COAD','LIHC','UCEC','LUSC','ESCA','STAD','OV','PAAD','HNSC','LUAD','LGG'};
xvalue={'MKerW-A','MKer-A','K-R','Cons-A','Cons-R','iCluster-A','S-R','SIM-R','P-A','C-A','SIM-A','K-M','S-A','Cons-M','Ker-A','SIM-M','K-A','K-C','Cons-C','SIM-C','S-M','S-C'};
h=heatmap(xvalue,yvalue,minheat);
h.Colormap=flipud(hot);
h.GridVisible='off';
h.Title='Area_Min';
h.FontSize = 12;
h.ColorLimits=[0.3 1];





