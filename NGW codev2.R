#load package
library(Momocs) #Momocs v1.2.9
library(geomorph)#Geomorph v3.2.1

#2. Methodology
#Digitisation test
#Load data and data frame
MethodCoords <- list.files("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainCoordinates\\MethodTest", full.names = TRUE)
MethodFrame <- read.csv("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainDataFrames\\MethodTest_Matrix.csv", header = TRUE)
#import coordinates
MethodTxt <- import_txt(MethodCoords , fileEncoding="UTF-8-BOM")
#create an Out-class object (i.e. set of outlines) from specified coordinate list (.txt files) and associated data frame
MethodOut <- Out(MethodTxt, fac=MethodFrame) 
#filter by grain-view and compute elliptic fourier transforms individually
MethodOut.l <- filter(MethodOut, View == "l")
MethodOut.l.efour <- efourier(MethodOut.l, nb.h=8, norm = FALSE, start = FALSE)# for why norm= false see note on ?efourier regarding roughly circular objects
MethodOut.d <- filter(MethodOut, View == "d")
MethodOut.d.efour <- efourier(MethodOut.d, nb.h=8, norm = FALSE, start = FALSE)# for why norm= false see note on ?efourier regarding roughly circular objects
#Split dataset by view, compute elliptical fourier analysis separately, recombine dataset and then run PCA 
MethodOut.d.l <- combine (MethodOut.d, MethodOut.l)
MethodTest<- MethodOut.d.l%>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% PCA (scale=FALSE, center= TRUE) 
MethodTest2<- MethodOut.d.l%>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine
#Plot PCA: corresponds to Figure 4
MethodTestPlot<-plot(MethodTest, cex= 1, zoom = 1, points= TRUE, labelspoints= FALSE,'Grain',palette = pal_qual_solarized)
#Analysis of variance using specimen (individual factor) and session as dependent variables
#Export and load coefficients of Elliptic Fourier analysis
export(MethodOut.l.efour)#Produces figures used in Coeffs3
export(MethodOut.d.efour)#Produces figures used in Coeffs4
Coeffs4_dorsal <- read.csv("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainDataFrames\\Coeffs4_dorsal.csv")
Coeffs3_lateral <- read.csv("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainDataFrames\\Coeffs3_lateral.csv")
#Calculate PCA scores for above
PCAObservation_d<-prcomp(Coeffs4_dorsal)$x 
PCAObservation_l<-prcomp(Coeffs3_lateral)$x 
#Creates 5 level session factor (consecutive)
sessionfactor<-as.factor(rep(1:5, 5))
#Create 5 level individual factor (repeated)
individualfactor<-gl(5, 5)
individualfactor
#Geomorph data frames for PCA computations of coefficients
gdf_d<-geomorph.data.frame(shape=PCAObservation_d) #ID: a vector containing the specimens ID 
gdf_l<-geomorph.data.frame(shape=PCAObservation_l) #ID: a vector containing the specimens ID 
#Procrustes ANOVA using session factor and then individual factor as dependent variables
summary(procD.lm(shape~sessionfactor, data=gdf_d))
summary(procD.lm(shape~sessionfactor, data=gdf_l))
#=session difference not significant for dorsal or lateral view
mod<-summary(procD.lm(shape~individualfactor, data= gdf_d))
mod2<-summary(procD.lm(shape~individualfactor, data= gdf_l))
mod
mod2
#=indiv difference significant for dorsal and lateral view
#dorsal measurement error (See Claude 2008)
s2within<-mswithin<-mod[[1]][2,3]
mod[[1]][2, 3]
MSamong<-mod[[1]][1, 3]
MSamong
s2among<-(MSamong-mswithin)/5
s2within/(s2within+s2among)*100#measurement error
#lateral measurement error
s2within<-mswithin<-mod2[[1]][2,3]
mod2[[1]][2, 3]
MSamong<-mod2[[1]][1, 3]
MSamong
s2among<-(MSamong-mswithin)/5
s2within/(s2within+s2among)*100#measurement error

#3. Results 
#Stage 1 Analysis: Modern Baseline Study 
#1.1 Uncharred modern grain outline analysis
#Load data and data frame
UncharredCoords <- list.files("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainCoordinates\\Stage1Uncharred", full.names = TRUE)
UncharredFrame <- read.csv("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainDataFrames\\Stage1Uncharred_Matrix.csv", header = TRUE)
#import coordinates
UncharredTxt <- import_txt(UncharredCoords, fileEncoding="UTF-8-BOM")
#create an Out-class object (i.e. set of outlines) from specified coordinate list (.txt files) and associated data frame
UncharredOut <- Out(UncharredTxt, fac=UncharredFrame)
#scale the Out-class by centroid size (default option for coo_scale)
UncharredOut1 <-coo_scale (UncharredOut)
#filter by grain-view and calibrate Calculate how many harmonics will be needed to gather x% of harmonic power
UncharredOut.l <- filter(UncharredOut1, View == "l")
calibrate_harmonicpower_efourier(UncharredOut.l,nb.h=12)
UncharredOut.d <- filter(UncharredOut1, View == "d")
calibrate_harmonicpower_efourier(UncharredOut.d,nb.h=12)
UncharredOut.p <- filter(UncharredOut1, View == "p")
calibrate_harmonicpower_efourier(UncharredOut.p,nb.h=12)
#Elliptic fourier transforms for individual grain views
UncharredOut.l.efour <- efourier(UncharredOut.l, nb.h=8, norm = FALSE, start = TRUE)
UncharredOut.d.efour <- efourier(UncharredOut.d, nb.h=8, norm = FALSE, start = TRUE)
UncharredOut.p.efour <- efourier(UncharredOut.p, nb.h=8, norm = FALSE, start = TRUE)
#Principal Components Analysis of a single grain view by taxon
#lateral view example, using 'taxon.code' as discriminating factor
UncharredOut.l.pca <- PCA(UncharredOut.l.efour, scale= FALSE, center= TRUE)
#generate plot of above
plot(UncharredOut.l.pca, cex= 1, zoom = 1, points= TRUE, 'taxon.code', labelspoints= FALSE)
#Linear discriminant analysis of single grain views (change filter group to calculate other views)
#lateral LDA by taxon to produce cross-validation table and plot 
UncharredOut.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
UncharredOut.l.lda <- LDA(UncharredOut.l.efour, "taxon.code", retain=0.95, norm= FALSE, start= FALSE)
plot(UncharredOut.l.lda,  labelspoints= FALSE)
#LDA of views in combination (example= dorsal + lateral + polar)
UncharredOut.d.l.p <- combine (UncharredOut.d, UncharredOut.l, UncharredOut.p)
#Split dataset by view, compute elliptical fourier analysis separately, recombine dataset and then run LDA with taxon code as discriminant factor
UncharredOut.d.l.p %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
#produce plot of above (corresponds to Figure 6)
UncharredOut.d.l.p.lda <- UncharredOut.d.l.p %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center = TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
plot(UncharredOut.d.l.p.lda, zoom= 1, rect.labelsgroups= TRUE, cex= 0.05, cex.labelsgroups = 0.7, center.origin= FALSE, points= FALSE)

#1.2 Charred modern grain outline analysis (see notes for 1.1- here dorsal and lateral views are combined in final LDA)
CharredCoords <- list.files("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainCoordinates\\Stage1Charred", full.names = TRUE)
CharredFrame <- read.csv("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainDataFrames\\Stage1Charred_Matrix.csv", header = TRUE)
CharredTxt <- import_txt(CharredCoords, fileEncoding="UTF-8-BOM")
CharredOut <- Out(CharredTxt, fac=CharredFrame)
CharredOut1 <-coo_scale (CharredOut)
CharredOut.l <- filter(CharredOut1, View == "l")
CharredOut.d <- filter(CharredOut1, View == "d")
CharredOut.p <- filter(CharredOut1, View == "p")
CharredOut.l.efour <- efourier(CharredOut.l, nb.h=8, norm = FALSE, start = TRUE)
CharredOut.d.efour <- efourier(CharredOut.d, nb.h=8, norm = FALSE, start = TRUE)
CharredOut.p.efour <- efourier(CharredOut.p, nb.h=8, norm = FALSE, start = TRUE)
CharredOut.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
CharredOut.l.lda <- LDA(CharredOut.l.efour, "taxon.code", retain=0.95, norm= FALSE, start= FALSE)
plot(CharredOut.l.lda,  labelspoints= FALSE)
CharredOut.d.l <- combine (CharredOut.d, CharredOut.l)
CharredOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
CharredOut.d.l.lda <- CharredOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center = TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
plot(CharredOut.d.l.lda, zoom= 1, rect.labelsgroups= TRUE, cex= 0.05, cex.labelsgroups = 0.7, center.origin= FALSE, points= FALSE)

#1.3 Addition of T. dicoccoides and use of two views only
CharredCoords_dccds <- list.files("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainCoordinates\\Stage1Charred_dccds", full.names = TRUE)
CharredFrame_dccds <- read.csv("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainDataFrames\\Stage1Charred_dccds_Matrix.csv", header = TRUE)
CharredTxt_dccds <- import_txt(CharredCoords_dccds, fileEncoding="UTF-8-BOM")
CharredOut_dccds <- Out(CharredTxt_dccds, fac=CharredFrame_dccds)
CharredOut1_dccds <-coo_scale (CharredOut_dccds)
CharredOut.l_dccds <- filter(CharredOut1_dccds, View == "l")
CharredOut.d_dccds <- filter(CharredOut1_dccds, View == "d")
CharredOut.l.efour_dccds <- efourier(CharredOut.l_dccds, nb.h=8, norm = FALSE, start = TRUE)
CharredOut.d.efour_dccds <- efourier(CharredOut.d_dccds, nb.h=8, norm = FALSE, start = TRUE)
CharredOut.d.l_dccds <- combine (CharredOut.d_dccds, CharredOut.l_dccds)
CharredOut.d.l_dccds %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
CharredOut.d.l.lda_dccds <- CharredOut.d.l_dccds %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center = TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
#corresponds to Figure 5
plot(CharredOut.d.l.lda_dccds, zoom= 1, rect.labelsgroups= TRUE, cex= 0.05, cex.labelsgroups = 0.7, center.origin= FALSE, points= FALSE)

#Stage 2 Analysis: Comparison with archaeobotanical material from Çatalhöyük
#2.1 Comparison of modern reference grains and known ('unmixed') archaeological grains from pure emmer or NGW samples are used
RefCoords <- list.files("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainCoordinates\\ReferenceGrains_CatalUnmixed", full.names = TRUE)
RefFrame <- read.csv("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainDataFrames\\ReferenceGrains_CatalUnmixed_Matrix.csv", header = TRUE)
RefTxt <- import_txt(RefCoords, fileEncoding="UTF-8-BOM")
RefOut <- Out(RefTxt, fac=RefFrame)
RefOut1 <-coo_scale (RefOut)
RefOut.l <- filter(RefOut1, View == "l")
RefOut.l.efour <- efourier(RefOut.l, nb.h=8, norm = FALSE, start = TRUE)
RefOut.d <- filter(RefOut1, View == "d")
RefOut.d.efour <- efourier(RefOut.d, nb.h=8, norm = FALSE, start = TRUE)
RefOut.d.l <- combine (RefOut.d, RefOut.l)
RefOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
RefOut.d.l.lda <- RefOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center = TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
plot(RefOut.d.l.lda, zoom= 1.5, rect.labelsgroups= TRUE, cex= 0.05, cex.labelsgroups = 0.7, center.origin= FALSE, points= FALSE)

#2.2 Use dataset from 2.1 to reclassify grains from 'mixed' deposits of less certain ID
UnclassCoords <- list.files("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainCoordinates\\CatalGrainsMixed", full.names = TRUE)
UnclassFrame <- read.csv("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainDataFrames\\CatalGrainsMixed_Matrix.csv", header = TRUE)
UnclassTxt <- import_txt(UnclassCoords, fileEncoding="UTF-8-BOM")
UnclassOut <- Out(UnclassTxt, UnclassFrame) #creates an Out object from a specified list- you can specify landmarks in "ldk"
UnclassOut1 <-coo_scale (UnclassOut)
UnclassOut.d <- filter(UnclassOut1, View == "d")
UnclassOut.l <- filter(UnclassOut1, View == "l")
UnclassOut.d.l <- combine (UnclassOut.d, UnclassOut.l)
UnclassOut.d.l.efour <- UnclassOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine
#reclassify the above using the LDA of the dorsal and lateral views of the reference grains from 2.1
reLDA(UnclassOut.d.l.efour, RefOut.d.l.lda)

#Stage 3: Stage 3: Analysis by phase
#3.1 Analysis of total combined dataset by phase
CombinedCoords <- list.files("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainCoordinates\\Combined", full.names = TRUE)
CombinedFrame <- read.csv("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainDataFrames\\Combined_Matrix_Phased.csv", header = TRUE)
CombinedTxt <- import_txt(CombinedCoords, fileEncoding="UTF-8-BOM")
CombinedOut <- Out(CombinedTxt, fac=CombinedFrame)
CombinedOut1 <-coo_scale (CombinedOut)
CombinedOut.l <- filter(CombinedOut1, View == "l")
CombinedOut.l.efour <- efourier(CombinedOut.l, nb.h=8, norm = FALSE, start = TRUE)
CombinedOut.d <- filter(CombinedOut1, View == "d")
CombinedOut.d.efour <- efourier(CombinedOut.d, nb.h=8, norm = FALSE, start = TRUE)
CombinedOut.d.l <- combine (CombinedOut.d, CombinedOut.l)
CombinedOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
CombinedOut.d.l.lda <- CombinedOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center = TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
#Corresponds to Figure 7
plot(CombinedOut.d.l.lda, zoom= 1.8, rect.labelsgroups= TRUE, cex= 0.05, cex.labelsgroups = 0.7, center.origin= FALSE, points= FALSE)

#3.2 Grains classified as NGW entered as 'unknowns' into LDA using wild and domesticated Timopheev's wheat as classifiers
#load 
TaraTimoCoords <- list.files("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainCoordinates\\TaraTimo", full.names = TRUE)
TaraTimoFrame <- read.csv("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainDataFrames\\TaraTimo_Matrix.csv", header = TRUE)
TaraTimoTxt <- import_txt(TaraTimoCoords, fileEncoding="UTF-8-BOM")
TaraTimoOut <- Out(TaraTimoTxt, fac=TaraTimoFrame)
TaraTimoOut1 <-coo_scale (TaraTimoOut)
TaraTimoOut.l <- filter(TaraTimoOut1, View == "l")
TaraTimoOut.l.efour <- efourier(TaraTimoOut.l, nb.h=8, norm = FALSE, start = TRUE)
TaraTimoOut.d <- filter(TaraTimoOut1, View == "d")
TaraTimoOut.d.efour <- efourier(TaraTimoOut.d, nb.h=8, norm = FALSE, start = TRUE)
TaraTimoOut.d.l <- combine (TaraTimoOut.d, TaraTimoOut.l)
TaraTimoOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
TaraTimoOut.d.l.lda <- TaraTimoOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center = TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
#Load grains classified as NGW
NewCoords <- list.files("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainCoordinates\\ClassifiedNew", full.names = TRUE)
NewFrame <- read.csv("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainDataFrames\\ClassifiedNew_Matrix.csv", header = TRUE)
NewTxt <- import_txt(NewCoords, fileEncoding="UTF-8-BOM")
NewOut <- Out(NewTxt, NewFrame) #creates an Out object from a specified list- you can specify landmarks in "ldk"
NewOut1 <-coo_scale (NewOut)
NewOut.d <- filter(NewOut1, View == "d")
NewOut.l <- filter(NewOut1, View == "l")
NewOut.d.l <- combine (NewOut.d, NewOut.l)
NewOut.d.l.efour <- NewOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine
#reclassify the above using the LDA of the dorsal and lateral views of the reference grains from 2.1
reLDA(NewOut.d.l.efour, TaraTimoOut.d.l.lda)
#stats of reclassifications
#load file tabulating number of grains classified as T. araraticum and T. timopheevii by period
PercentTara <- read.csv(file="C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\ClassificationFiles\\NewType_Classifications.csv")
#generate table
Taratab <- table(PercentTara$Phase2,PercentTara$NewClassific)
Taratab
#chi-square test (early vs middle and late combined)
chisq.test(Taratab) 
#3.3: Visualise grains by phase
#example: visualise NGW from early phases
# e.g. filter out dorsal-view of NGW grains from early phase 
CombinedOut.d.new_early <-filter (CombinedOut.d, taxon.code == "Early")
#centre, scale by centroid size and stack outlines
CombinedOut.d.new_early %>% coo_center %>% coo_scale %>% stack

#Stage 4: Analysis by processing stage
#Stage 3: Stage 3: Analysis by phase
ProcCoords <- list.files("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainCoordinates\\Combined_noeinkorn", full.names = TRUE)
ProcFrame <- read.csv("C:\\Users\\Tina Roushannafas\\Documents\\Supplementary Data\\GrainDataFrames\\Combined_noeinkorn_Matrix_Process.csv", header = TRUE)
ProcTxt <- import_txt(ProcCoords, fileEncoding="UTF-8-BOM")
ProcOut <- Out(ProcTxt, fac=ProcFrame)
ProcOut1 <-coo_scale (ProcOut)
ProcOut.l <- filter(ProcOut1, View == "l")
ProcOut.l.efour <- efourier(ProcOut.l, nb.h=8, norm = FALSE, start = TRUE)
ProcOut.d <- filter(ProcOut1, View == "d")
ProcOut.d.efour <- efourier(ProcOut.d, nb.h=8, norm = FALSE, start = TRUE)
ProcOut.d.l <- combine (ProcOut.d, ProcOut.l)
ProcOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
ProcOut.d.l.lda <- ProcOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center = TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
#Corresponds to Figure 10
plot(ProcOut.d.l.lda, zoom= 1.4, rect.labelsgroups= TRUE, cex= 0.05, cex.labelsgroups = 0.7, center.origin= FALSE, points= FALSE)


