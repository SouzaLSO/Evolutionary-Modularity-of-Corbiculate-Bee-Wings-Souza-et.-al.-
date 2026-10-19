####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################### ANALISES MORFOMETRICAS ################################################

#### library Packs
library("ape")
library("geomorph")
library("dplyr")
library("dispRity")

#Load the "Rdata"
#Note that some terms are in Portuguese
##### Procrustes superimposition (GPA)
marcos_GPA
plot(marcos_GPA)
#plotOutliers(marcos_GPA$coords)
#dim(marcos_GPA$coords)
#### Tabela de classificadores 
#class_ID <-read.csv("Class_ID.csv", sep=";", header=T)


#### Loading adjusted tree
Phy_tree$tip.label
plot(Phy_tree, cex = 0.5, no.margin = T)
#plot(Phy_tree, cex = 0.3, adj = 0.1, no.margin = T)

#### Identifying classifiers
length(unique(class_ID$TPS)) #n de imagens marcadas (individuos)
length(unique(class_ID$especie)) # especies 
length(unique(class_ID$genero)) # generos
length(unique(class_ID$tribo)) #tribos
length(unique(class_ID$class)) # tribos com divisão Melipona e trigoniformes


############### Calculating the mean shape per species #################
##### Selecting species as the primary classifier
spp_names<-class_ID$especie
head(spp_names)
new_coords = coords.subset(A = marcos_GPA$coords, group = spp_names)

forma_media = lapply(new_coords, mshape) #Rodar o comando "mshape" em cada linha
ordem_mshape<-names(forma_media) #ordem após forma média

new_dim<-dim(marcos_GPA$coords)

msh_sp<-array(unlist(forma_media), dim = c(new_dim[1], new_dim[2], length(unique(spp_names))))
dim(msh_sp)
#We now have a mean-shape dataset with 47 landmarks in 2 dimensions for the 99 species.

##### New classifiers based on the mean shape for each species
View(new_class)
gen<-new_class$genero
nomes_especies <- new_class$especie
dimnames(msh_sp)[[3]]<-nomes_especies
nome_tribos<-new_class$tribo
trib_class<-new_class$class

# Simple disparity between groups
gdf <- geomorph.data.frame(forma = msh_sp, taxa = new_class$class)
gdf
estrut = morphol.disparity(forma ~ taxa, groups= ~ taxa, data = gdf, iter = 999, print.progress = FALSE)
str(estrut)
disparidade = estrut$Procrustes.var
#write.csv(disparidade, "disparidade_Class.csv")


##################### Start of Analyses
#PCA_Morphospace
#TPS
#Disparity over time
#Modularity
#Evolutionary models of the wing
#Disparity and evolutionary models of the best-supported modularity hypothesis


##################################################################
################################  PCA  ##########################
##Colors for genders
{
  #COLORING BY TRIBE
  col.gp1 <- colorRampPalette (c("#F7F702", "#7A12AD", "#05F50F", "#03A9F4", "#000000"))
  col.gp2<-col.gp1(5)
  names(col.gp2) <- levels(as.factor(trib_class))
  col.gp<- col.gp2[match(trib_class, names(col.gp2))] # col.gp must NOT be a factor
  col.gp
 
  ##Cores para generos
  {
    col.gpGen <- colorRampPalette (c("#CA194E", "#F7F702", "#38fcf6", "#FF33A6", "#7A12AD", 
                                     "#33FFA6", "#d581eb", "#0A6348", "#05F50F", "#fcac0a", 
                                     "#b82e33", "#ffffff", "#37b302", "#FF1100", "#33A6FF",
                                     "#8E24AA", "#FFC107", "#009688", "#F44336", "#03A9F4", 
                                     "#FF9800", "#A6FF33", "#474b61", "#FF5722", "#cf927c", 
                                     "#CDDC39", "#000000"))
    col.gpGen<-col.gpGen(27)
    names(col.gpGen) <- levels(as.factor(gen))
    col.gpGen<- col.gpGen[match(gen, names(col.gpGen))] # col.gp must NOT be a factor
    col.gpGen
  }
  
  
  ## Cores do elipse para tribos
  colorT<-colorRampPalette (c("#F7F702", "#8E24AA","#05F50F", "blue"))
  colorT<-colorT(4)
  names(colorT)<-levels(as.factor(nome_tribos))
  
}

#### Phylomorpho PCA
# Note that the summary statistics obtained for this analysis are identical to those of standard PCA. This is because, in this case, the phylogeny is used merely for plotting and is NOT taken into account during the analytical procedures.
# https://cran.r-project.org/web/packages/geomorph/vignettes/geomorph.PCA.html

PCA_Mphy<-gm.prcomp(msh_sp, phy = Phy_tree)
summary(PCA_Mphy)

plot(PCA_Mphy, main= "PhylomorphospacePCA", phylo= TRUE, pch = 21, bg = col.gp, cex = 1.3, cex.lab = 1, font.lab = 2)

phy.par <- list(tip.labels = FALSE, node.labels = FALSE, edge.color = 1)

plot(PCA_Mphy, phylo= TRUE, main= "PhylomorphospacePCA",phylo.par = phy.par,
     pch = 21, bg = col.gp, cex = 1.3, cex.lab = 1.5, font.lab = 2)

plot(PCA_Mphy$x[, 1], main = "PC1", pch = 21, bg = col.gp, cex = 1.3, cex.lab = 1, font.lab = 2)

plot(PCA_Mphy$x[, 2], main = "PC2", pch = 21, bg = col.gp, cex = 1.3, cex.lab = 1, font.lab = 2)

#PCA with phylogenetic signal removed
Phylogenetic_PCA<-gm.prcomp(msh_sp, phy = Phy_tree, GLS = TRUE)
summary(Phylogenetic_PCA)

plot(Phylogenetic_PCA, flip = 2, phylo= TRUE, main= "Phylogenetic PCA",phylo.par = phy.par, pch = 21, bg = col.gp, cex = 1.3, cex.lab = 1.5, font.lab = 2)

#COLORIDO POR GÊNERO
{
  plot(PCA_Mphy, phylo= TRUE, main= "PhylomorphospacePCA",phylo.par = phy.par,
       pch = 21, bg = col.gpGen, cex = 1.5, cex.lab = 1.5, font.lab = 2)
  
  for(g in levels(as.factor(nome_tribos))){
   pontos <- PCA_Mphy$x[nome_tribos == g, 1:2]  # pega PC1 e PC2
    h <- chull(pontos)                           # convex hull
    h <- c(h, h[1])                              # fecha o polígono
    polygon(pontos[h, ], border = colorT[g], col = adjustcolor(colorT[g], alpha.f = 0.2), lwd = 2)}
}


###### 3D filo
plot(PCA_Mphy, time.plot = TRUE, pch = 21, bg = col.gp, rep= col.gp, cex = 1.5,
     phylo.par = list(edge.color = col.gp, edge.width = 1.5, tip.txt.cex = 0.75,
                      node.labels = F, anc.states = F))
#Salva
#library(htmlwidgets)
#rglwidget()
#widget <- rglwidget()
#saveWidget(widget, "PCAPhylo.html", selfcontained = TRUE)


######### TPS_ GRADES DE DEFORMAÇÃO dos PCA's
msh = mshape(msh_sp) #Forma média

# define which landmarks are connected by lines in wireframe graph
link.x <- c( 1,1,2,3,  4,5,6,6,      14:46,47,11,10,9,8)
link.y <- c(2,10,3,4,  5,14,14,17,  15:47,11,10,9,8,7)
links<- cbind(link.x, link.y)

# mark minimum blue and maximum red
GPBlue <- gridPar(pt.bg = "gray", link.col = "black", pt.size = 1, 
                  tar.pt.bg ="blue", tar.link.col ="blue")

GPRed <- gridPar(pt.bg = "gray", link.col = "black", pt.size = 1, 
                 tar.pt.bg ="red", tar.link.col ="red")


# PCA_ PC1
plotRefToTarget(M1 = msh,
                M2 = PCA_Mphy$shapes$shapes.comp1$min,
                gridPars=GPBlue, method = "point", links = links)

plotRefToTarget(M1 = msh,
                M2 = PCA_Mphy$shapes$shapes.comp1$max,
                gridPars=GPRed, method = "point", links = links)
#PC2
plotRefToTarget(M1 = msh,
                M2 = PCA_Mphy$shapes$shapes.comp2$min,
                gridPars=GPBlue, method = "point", links = links)

plotRefToTarget(M1 = msh,
                M2 = PCA_Mphy$shapes$shapes.comp2$max,
                gridPars=GPRed, method = "point", links = links)


################################################################
############  Disparidade al longo do tempo DTT ################

# Extract the first two components 
autoval <- PCA_Mphy$x[, 1:3]
View(autoval)

average.sq <- function(X) mean(pairwise.dist(X)^2) # metrica

dispRity_dtt <- dtt.dispRity(data = autoval, metric = average.sq, tree = Phy_tree, nsim = 200, model = "BM")
plot(dispRity_dtt, xlim = rev(range( dispRity_dtt$times))) #Inverter escala do eixo do tempo

dispRity_dtt$MDI

################################################################
####################### Modularity Analysis ################
################################################################


# 1. Criar a cofiguração dos módulos hipotéticos
# 2. Criar matriz de correlação entre os landmarks x espécies
# 3. Testar modelos de módulos com a força entre e dentro de cada modulo (EMMLi)
#Verificar o valor de CR do melhor modelo de phylo.modularity

#### Module hypotheses
{
  #H00 Null Hypothesis
  H00<- c("A","A","A","A","A",
          "A","A","A","A","A", #10
          "A","A","A","A","A",
          "A","A","A","A","A", #20
          "A","A","A","A","A",
          "A","A","A","A","A", #30
          "A","A","A","A","A",
          "A","A","A","A","A", #40
          "A","A","A","A","A",
          "A","A")
  
  ## Módulo Antero-posterior
  #H01 (Antero(A) + Posterior(T))
  H01a<-c("T","T","A","A","A", 
          "A","A","T","T","T", #10
          "T","A","T","A","A", 
          "A","A","A","A","A", #20
          "A","A","A","A","A", 
          "A","A","A","A","A", #30
          "A","T","T","T","T", 
          "T","T","T","T","T", #40
          "T","T","T","T","T",
          "T","T")
  
  H01b<-c("T","T","T","A","A",
          "A","T","T","T","T", #10
          "T","A","T","A","A",
          "A","A","A","A","A", #20
          "A","A","A","A","A",
          "A","A","T","T","T", #30
          "T","T","T","T","T",
          "T","T","T","T","T", #40
          "T","T","T","T","T",
          "T","T")
  
  ## Módulo Proximal-Distal
  #H02 (Proximal(P)= 1;2 Distal(D))
  H02a<-  c("P","P","D","D","D",
            "D","D","D","D","D",
            "D","D","D","D","D",
            "D","D","D","D","D",
            "D","D","D","D","D",
            "D","D","D","D","D",
            "D","D","D","D","D",
            "D","D","D","D","D",
            "D","D","D","D","D",
            "D","D")
  
  #H04 (Parte distal(D) menor)
  H02b<-c("P","P","P","P","P",
          "P","P","P","P","P", #10
          "P","D","D","P","P",
          "P","P","P","P","D", #20
          "D","D","D","D","D",
          "D","D","D","D","D", #30
          "D","D","D","D","D",
          "D","D","D","D","D", #40
          "D","D","D","D","D",
          "P","P")
  
  H02c<-c("P","P","P","P","P",
          "P","P","P","P","P", #10
          "P","P","D","P","P",
          "P","P","P","P","P", #20
          "P","P","P","P","P",
          "D","D","D","D","D", #30
          "D","D","D","D","D",
          "D","D","D","D","D", #40
          "D","D","D","P","P",
          "P","P")
  
  #Parte individuais como modulos
  # Modulos do estigma
  #H03 (Asa(A) + estigma(S))
  H03<-c("A","A","A","A","S",
         "S","A","A","A","A", #10
         "A","A","A","S","S",
         "S","S","A","A","A", #20
         "A","A","A","A","A",
         "A","A","A","A","A",
         "A","A","A","A","A",
         "A","A","A","A","A",
         "A","A","A","A","A",
         "A","A")
  
  #H4 Três módulos (Posterior (T) + estigma(S) + Anterior(A))
  H04<-c("T","T","T","T","S",
         "S","T","T","T","T", #10
         "T","A","T","S","S",
         "S","S","A","A","A", #20
         "A","A","A","A","A",
         "A","A","A","T","T", #30
         "T","T","T","T","T",
         "T","T","T","T","T", #40
         "T","T","T","T","T",
         "T","T")
  
  #H5 Três módulos (Proximal (P) + estigma(S) + Distal(D))
  H05<-c("P","P","P","P","S",
         "S","D","D","D","D", #10
         "D","D","D","S","S",
         "S","S","D","D","D", #20
         "D","D","D","D","D",
         "D","D","D","D","D", #30
         "D","D","D","D","D",
         "D","D","D","D","D", #40
         "D","D","D","D","D",
         "D","D")
  
  
  # Quanto modulos
  #H06 4 módulos (Anterior (A), Proximal (P) + Distal(D) + stigma(s))
  H06<-c("P","P","P","P","S",
         "S","P","P","P","P", #10
         "P","A","D","S","S",
         "S","S","A","A","A", #20
         "A","A","A","A","A",
         "D","D","D","D","D", #30
         "D","D","D","D","D",
         "D","D","D","D","D", #40
         "D","D","D","P","P",
         "P","P")
  
  
  
  #H07 4 módulos (Proxima(P) +  Anterior(A)+ Distal(D) + Traseiro (T))
  H07<-c("P","P","P","P","A",
         "A","P","P","P","T", #10
         "T","A","D","A","A",
         "A","A","A","A","A", #20
         "A","A","A","A","A",
         "A","A","D","D","D", #30
         "D","D","D","D","D",
         "D","D","T","T","T", #40
         "T","T","T","T","T",
         "T","T")
  
  
  #H08 (Proxima(P) + Estigma(S)+  Anterior(A)+ Distal(D) + Traseiro (T))
  H08<-c("P","P","P","P","S",
         "S","P","P","P","T", #10
         "T","A","D","S","S",
         "S","S","A","A","A", #20
         "A","A","A","A","A",
         "A","A","A","D","D", #30
         "D","D","D","D","D",
         "D","D","T","T","T", #40
         "T","T","T","T","T",
         "T","T")
  
  
  #Para desenhar os módulos
  #mod_four<-define.modules(spec=marcos_GPA$consensus, nmodules= 4)
}
# Criando Tabela (dataframe) com Módulos hipotéticos
Modulos_Hipoteticos<-data.frame(H00,H01a,H01b,H02a,H02b,H02c,H03,H04,H05,H06,H07,H08)
#View(Modulos_Hipoteticos)

############# Evaluating module hypotheses using a maximum likelihood approach

library("EMMLiv2")
library("paleomorph")
install.packages("caper")
library("caper")

class(msh_sp)
dim(msh_sp)

#philoMEL creates its own correlation matrix
test_hipotesesPhy <- EMMLiv2::phyloEmmli(landmarks = msh_sp, phylo = Phy_tree, method = "pgls", EMMLi = T, N_sample =  99, mod= Modulos_Hipoteticos)


test_hipotesesPhy # Se letras = Ordem alfabética (P = 3, S = 4, A = 1, D = 2,)

Best_Model<-(test_hipotesesPhy$EMMLi$rho[1])
Best_Model

#K= N de parametros estimados

###################### Checking and comparing the CR value from phylo.modularity
# CR < 1 Modular
# CR > 1 Não modular
# https://doi.org/10.1111/2041-210X.12511
# Rodando o teste para cada hipótese
{
  CR01a<-phylo.modularity(msh_sp,
                          partition.gp = H01a, #Indicar
                          phy = Phy_tree,
                          iter = 999,
                          seed = NULL,
                          print.progress = TRUE)
  
  CR01b<-phylo.modularity(msh_sp,
                          partition.gp = H01b, #Indicar
                          phy = Phy_tree,
                          iter = 999,
                          seed = NULL,
                          print.progress = TRUE)
  
  CR02a<-phylo.modularity(msh_sp,
                          partition.gp = H02a, #Indicar
                          phy = Phy_tree,
                          iter = 999,
                          seed = NULL,
                          print.progress = TRUE)
  
  CR02b<-phylo.modularity(msh_sp,
                          partition.gp = H02b, #Indicar
                          phy = Phy_tree,
                          iter = 999,
                          seed = NULL,
                          print.progress = TRUE)
  
  CR02c<-phylo.modularity(msh_sp,
                          partition.gp = H02c, #Indicar
                          phy = Phy_tree,
                          iter = 999,
                          seed = NULL,
                          print.progress = TRUE)
  
  CR03<-phylo.modularity(msh_sp,
                         partition.gp = H03, #Indicar
                         phy = Phy_tree,
                         iter = 999,
                         seed = NULL,
                         print.progress = TRUE)
  
  
  CR04<-phylo.modularity(msh_sp,
                         partition.gp = H04, #Indicar
                         phy = Phy_tree,
                         iter = 999,
                         seed = NULL,
                         print.progress = TRUE)
  
  CR05<-phylo.modularity(msh_sp,
                         partition.gp = H05, #Indicar
                         phy = Phy_tree,
                         iter = 999,
                         seed = NULL,
                         print.progress = TRUE)
  
  CR06<-phylo.modularity(msh_sp,
                         partition.gp = H06, #Indicar
                         phy = Phy_tree,
                         iter = 999,
                         seed = NULL,
                         print.progress = TRUE)
  
  CR07<-phylo.modularity(msh_sp,
                         partition.gp = H07, #Indicar
                         phy = Phy_tree,
                         iter = 999,
                         seed = NULL,
                         print.progress = TRUE)
  
  CR08<-phylo.modularity(msh_sp,
                         partition.gp = H08, #Indicar
                         phy = Phy_tree,
                         iter = 999,
                         seed = NULL,
                         print.progress = TRUE)
}
# Compare the modulation effects, including the null model
resultado_comparacao <- compare.CR(CR01a, CR01b, CR02a, CR02b, CR02c, CR03, CR04, CR05, CR06, CR07, CR08, CR.null = TRUE)


# Exibir o resumo dos resultados
summary(resultado_comparacao)
str(resultado_comparacao)


#####################################################################
################### Evolutionary models ###########################

################### Modelos evolutivos ###########################

#install.packages("geiger")
library(mvMORPH)
library(geiger)
#install.packages("RPANDA")



# Defina os dados de tamanho
size <- new_class$max_size

# Divide into quintiles using the cut function
#size_categories <- cut(size, quantile(size, probs = seq(0, 1, 0.2)), labels = c("Very Small", "Small", "Medium", "Large", "Very Large"), include.lowest = TRUE)
#size_categories

#new_class$size<-size_categories

#Creating SIZE configuration
tamanho<-size
tamanho<-new_class$size
tamanho
names(tamanho)<-(new_class$especie)
tamanho

#Link the tree
tree5 = make.simmap(Phy_tree, tamanho)
col5=c("blue","gray", "black" ,"orange","red")
names(col5)=c("Very Small", "Small", "Medium", "Large", "Very Large")
plotSimmap(tree5, col5, fsize=0.3, node.numbers = FALSE, lwd=2)

# Eusocial (s) non-eusocial (n)
Eussocial<-new_class$eussocial
Eussocial
names(Eussocial)<-(new_class$especie) #Vincula grupo ao nome da sp
Eussocial

#Link the tree 
tree2 = make.simmap(Phy_tree, Eussocial)
col2=c("red","blue")
names(col2)=c("n","s")
plotSimmap(tree2, col2, fsize=0.3, node.numbers = FALSE, lwd=2)


## Second Test
#Behavior (Solitary, subsocial, eusocial)
social<-new_class$social
social
names(social)<-(new_class$especie)
social

#Link the tree
tree3 = make.simmap(Phy_tree, social)
col3=c("blue","green","red")
names(col3)=c("eussocial","solitaria","s.primitiva")
plotSimmap(tree3, col3, fsize=0.3, node.numbers = FALSE, lwd=2)

## Third (Solitary, isolated)
solitario<-new_class$solitario
solitario
names(solitario)<-(new_class$especie) #Vincula grupo ao nome da sp
solitario

#Link the tree 
tree4 = make.simmap(Phy_tree, solitario)
col4=c("green","blue")
names(col4)=c("solitaria","social")
plotSimmap(tree4, col4, fsize=0.3, node.numbers = FALSE, lwd=2)

### Create a table with the key components for analysis
summary(PCA_Mphy)
autoval <- PCA_Mphy$x[,1:3] #PC1 PC2 PC3

#### Model testing

### Individual models
#BM - one taxa
BM = mvBM(tree2, autoval, model = "BM1")
#EB
EB = mvEB(tree2, autoval)
#OU
OU1 = mvOU(tree2, autoval, model = "OU1")

###### Isolated eusocial behavior
#BMM - múltiplas taxas
Eusoc_BMM = mvBM(tree2, autoval, model = "BMM")
#OUM com múltiplos ótimos
Eusoc_OUM = mvOU(tree2, autoval, model = "OUM")

##### Individual Social Behavior (A+M)x(B)x(E)
#BMM - múltiplas taxas
Soc_BMM = mvBM(tree3, autoval, model = "BMM")
#OU com múltiplos ótimos
Soc_OUM = mvOU(tree3, autoval, model = "OUM")
str(Soc_OUM)
##### Comportamento solitario isolado
#BMM - múltiplas taxas
Sol_BMM = mvBM(tree4, autoval, model = "BMM")
#OUM com múltiplos ótimos
Sol_OUM = mvOU(tree4, autoval, model = "OUM")

##### SIZE
size_BMM = mvBM(tree5, autoval, model = "BMM")
#OUM com múltiplos ótimos
size_OUM = mvOU(tree5, autoval, model = "OUM")



#Comparing models
nomes_modelos<-c("BM","EB","OU1","Eusoc_BMM","Eusoc_OUM","Soc_BMM","Soc_OUM", "Sol_BMM","Sol_OUM","size_BMM","size_OUM") ## CUIDADO COM A ORDEM

#Table with the results
Model_tes<-as.data.frame(c("AIC","AICc"))
rownames(Model_tes)<-Model_tes[, 1]
Model_tes<-Model_tes[,-1]

## BEWARE OF THE ORDER
Model_tes$BM<-c(BM$AIC,BM$AICc)
Model_tes$EB<-c(EB$AIC,EB$AICc)
Model_tes$OU1<-c(OU1$AIC,OU1$AICc)


## Eussocial isolado
Model_tes$Eusoc_BMM<-c(Eusoc_BMM$AIC,Eusoc_BMM$AICc)
Model_tes$Eusoc_OUM<-c(Eusoc_OUM$AIC,Eusoc_OUM$AICc)
# Solcial individual
Model_tes$Soc_BMM<-c(Soc_BMM$AIC,Soc_BMM$AICc)
Model_tes$Soc_OUM<-c(Soc_OUM$AIC,Soc_OUM$AICc)
## Solitario isolado
Model_tes$Sol_BMM<-c(Sol_BMM$AIC,Sol_BMM$AICc)
Model_tes$Sol_OUM<-c(Sol_OUM$AIC,Sol_OUM$AICc)
## Tamanho
Model_tes$size_BMM<-c(size_BMM$AIC,size_BMM$AICc)
Model_tes$size_OUM<-c(size_OUM$AIC,size_OUM$AICc)

Model_tabela<-as.data.frame(t(Model_tes)) #Inverter linhas e colunas
View(Model_tabela)
# Calcular ΔAICc, w e AICw 

Result_AICc<-mvMORPH::aicw(Model_tabela$AICc, aicc= TRUE, models=nomes_modelos) #calculo
Result_AICc
Result_AICc$models<-nomes_modelos ### CUIDADO COM A ORDEM
Result_AICc

#Adding the table

Model_tabela$AIC<-round(Model_tabela$AIC, 0)
Model_tabela$AICc<-round(Result_AICc$AIC, 0)
Model_tabela$ΔAICc<-round(Result_AICc$diff,2)
Model_tabela$w<-round(Result_AICc$wi,3)
Model_tabela$AICcw<-round(Result_AICc$aicweights,3)
# Visualizar o resultado
print(Model_tabela)
#write.csv(Model_tabela,"Modelos_Evolutivos com 6 eixos.csv") #Salvando tabela
#Medir o peso de Akaike para ver qual o melhor
#Resultado OU+Multiplos ótimos para as duas análises

# Melhor modelo foi o Soc_OUM 
#Análises de Ornstein-Uhlenbeck
meia_vida<-halflife(Soc_OUM)
meia_vida
### A meia-vida filogenética descreve o tempo para passar da metade do estado ancestral até o ótimo primário (Hansen, 1997).
estacionario<-stationary(Soc_OUM)
estacionario
#The expected variance when the process is stationary. The multivariate normal stationary distribution of the Ornstein-Uhlenbeck process is calculated following Bartoszek et al. (2012).

plot(Soc_OUM$param)
str(Soc_OUM)

Soc_OUM

Soc_OUM$theta
Soc_OUM$alpha
Soc_OUM$sigma

Soc_OUM$param$sigma
Soc_OUM$param$alpha

str(Soc_OUM$param)
Soc_OUM$param$names_regimes
#Theta<- Soc_OUM$theta

Theta #valores de theta

########## Now I will plot the PCA using the Theta values.
plot(PCA_Mphy, phylo = TRUE, main = "PCA_phy", phylo.par = phy.par,
     pch = 21, bg = col.gp, cex = 1.3, cex.lab = 1.5, font.lab = 2)

points(Theta,  pch = 23, bg = "red", cex = 1.8)
text(Theta, labels = rownames(Theta), pos = 3, cex = 1.1, font = 2)



##############################################################################
##############################################################################
##############################################################################
################## INDIVIDUAL MODULES Script ###############################

{
# Separar Marcos dos módulos
Mod_A<- msh_sp[c(12,18:25),,]
Mod_D<- msh_sp[c(13,26:38),,]
Mod_P<- msh_sp[c(1:4,7:11,39:47),,]
Mod_S<- msh_sp[c(5:6,14:17),,]

dim(Mod_A)
dim(Mod_D)
dim(Mod_P)
dim(Mod_S)

# Modelos já estabelecidos
social
plotSimmap(tree3, col3, fsize=0.3, node.numbers = FALSE, lwd=2)

Eussocial
plotSimmap(tree2, col2, fsize=0.3, node.numbers = FALSE, lwd=2)

solitario
plotSimmap(tree4, col4, fsize=0.3, node.numbers = FALSE, lwd=2)

tamanho
plotSimmap(tree5, col5, fsize=0.3, node.numbers = FALSE, lwd=2)


### Crie uma tabela com os componentes principais para análises
PCA_Mod_A<-gm.prcomp(Mod_A, phy = Phy_tree)
PCA_Mod_D<-gm.prcomp(Mod_D, phy = Phy_tree)
PCA_Mod_P<-gm.prcomp(Mod_P, phy = Phy_tree)
PCA_Mod_S<-gm.prcomp(Mod_S, phy = Phy_tree)


#Plot PCA Modulos individuais
par(mfrow=c(2,2))
phy.par <- list(tip.labels = FALSE, node.labels = FALSE, edge.color = 1)

plot(PCA_Mod_A,flip = 1, phylo= TRUE, main= "Module A",phylo.par = phy.par,
     pch = 21, bg = col.gp, cex = 1.3, cex.lab = 1.5, font.lab = 2)

plot(PCA_Mod_D,flip = 1, phylo= TRUE, main= "Module D",phylo.par = phy.par,
     pch = 21, bg = col.gp, cex = 1.3, cex.lab = 1.5, font.lab = 2)

plot(PCA_Mod_P,flip = 1, phylo= TRUE, main= "Module P",phylo.par = phy.par,
     pch = 21, bg = col.gp, cex = 1.3, cex.lab = 1.5, font.lab = 2)

plot(PCA_Mod_S,flip = 1, phylo= TRUE, main= "Module S",phylo.par = phy.par,
     pch = 21, bg = col.gp, cex = 1.3, cex.lab = 1.5, font.lab = 2)

par(mfrow=c(1,1))

# Autovalores para modelos evolutivos dos modulos
auto_Mod_A<-PCA_Mod_A$x[,1:3]
auto_Mod_D<-PCA_Mod_D$x[,1:3] 
auto_Mod_P<-PCA_Mod_P$x[,1:4]
auto_Mod_S<-PCA_Mod_S$x[,1:3] 

# Rodando modelos para modulos

#Modulo A #Result Soc_OUM 
{
  #BM - uma taxa
  BM_ModA = mvBM(tree2, auto_Mod_A, model = "BM1")
  #EB
  EB_ModA = mvEB(tree2, auto_Mod_A)
  #OU
  OU1_ModA = mvOU(tree2, auto_Mod_A, model = "OU1")
  
  ###### Eussocial ou Não Eussocial (A+M)x(B+E)
  #BMM - múltiplas taxas
  Eusoc_BMM_ModA = mvBM(tree2, auto_Mod_A, model = "BMM")
  #OUM com múltiplos ótimos
  Eusoc_OUM_ModA = mvOU(tree2, auto_Mod_A, model = "OUM")
  
  
  ##### Social (A+M)x(B)x(E)
  #BMM - múltiplas taxas
  Soc_BMM_ModA = mvBM(tree3, auto_Mod_A, model = "BMM")
  #OU com múltiplos ótimos
  Soc_OUM_ModA = mvOU(tree3, auto_Mod_A, model = "OUM")
  
  Sol_BMM_ModA = mvBM(tree4, auto_Mod_A, model = "BMM")
  #OUM com múltiplos ótimos
  Sol_OUM_ModA = mvOU(tree4, auto_Mod_A, model = "OUM")
  
  ##### Tamanho
  size_BMM_ModA = mvBM(tree5, auto_Mod_A, model = "BMM")
  #OUM com múltiplos ótimos
  size_OUM_ModA = mvOU(tree5, auto_Mod_A, model = "OUM")
  
  
  
  #Comparando modelos
  nomes_modelos<-c("BM","EB","OU1","Eusoc_BMM","Eusoc_OUM","Soc_BMM","Soc_OUM", "Sol_BMM","Sol_OUM","size_BMM","size_OUM") ## CUIDADO COM A ORDEM
  
  #Tabela com os resultados
  Modu_Ates<-as.data.frame(c("Lik","AIC","AICc"))
  rownames(Modu_Ates)<-Modu_Ates[, 1]
  Modu_Ates<-Modu_Ates[,-1]
  
  Modu_Ates$BM_ModA<-c(BM_ModA$LogLik,BM_ModA$AIC,BM_ModA$AICc)
  Modu_Ates$EB_ModA<-c(EB_ModA$LogLik,EB_ModA$AIC,EB_ModA$AICc)
  Modu_Ates$OU1_ModA<-c(OU1_ModA$LogLik,OU1_ModA$AIC,OU1_ModA$AICc)
  
  Modu_Ates$Eusoc_BMM_ModA<-c(Eusoc_BMM_ModA$LogLik,Eusoc_BMM_ModA$AIC,Eusoc_BMM_ModA$AICc)
  Modu_Ates$Eusoc_OUM_ModA<-c(Eusoc_OUM_ModA$LogLik,Eusoc_OUM_ModA$AIC,Eusoc_OUM_ModA$AICc)
  
  Modu_Ates$Soc_BMM_ModA<-c(Soc_BMM_ModA$LogLik,Soc_BMM_ModA$AIC,Soc_BMM_ModA$AICc)
  Modu_Ates$Soc_OUM_ModA<-c(Soc_OUM_ModA$LogLik,Soc_OUM_ModA$AIC,Soc_OUM_ModA$AICc)
  
  Modu_Ates$Sol_BMM_ModA<-c(Sol_BMM_ModA$LogLik,Sol_BMM_ModA$AIC,Sol_BMM_ModA$AICc)
  Modu_Ates$Sol_OUM_ModA<-c(Sol_OUM_ModA$LogLik,Sol_OUM_ModA$AIC,Sol_OUM_ModA$AICc)
  
  Modu_Ates$size_BMM_ModA<-c(size_BMM_ModA$LogLik,size_BMM_ModA$AIC,size_BMM_ModA$AICc)
  Modu_Ates$size_OUM_ModA<-c(size_OUM_ModA$LogLik,size_OUM_ModA$AIC,size_OUM_ModA$AICc)
  
  Modu_Atabela<-as.data.frame(t(Modu_Ates)) #Inverter linhas e colunas
  
  # Calcular ΔAICc, w e AICw 
  
  Mod_AAICc<-mvMORPH::aicw(Modu_Atabela$AICc) #calculo
  Mod_AAICc$models<-nomes_modelos
  Mod_AAICc
  
  #Adicionando a tabela
  
  Modu_Atabela$AIC<-round(Modu_Atabela$AIC, 0)
  Modu_Atabela$AICc<-round(Mod_AAICc$AIC, 0)
  Modu_Atabela$ΔAICc<-round(Mod_AAICc$diff,2)
  Modu_Atabela$w<-round(Mod_AAICc$wi,3)
  Modu_Atabela$AICcw<-round(Mod_AAICc$aicweights,3)
  
  print(Modu_Atabela)
  
  write.csv(Modu_Atabela,"ModEv_Mod_A.csv") #Salvando tabela
  Mod_AAICc
  
  Parametros_ModA_SocOUM<- Soc_OUM_ModA$theta
  
}

#Modulo D #Result Soc OUM (diff < 2)
{
  #BM - uma taxa
  BM_ModD = mvBM(tree2, auto_Mod_D, model = "BM1")
  #EB
  EB_ModD = mvEB(tree2, auto_Mod_D)
  #OU
  OU1_ModD = mvOU(tree2, auto_Mod_D, model = "OU1")
  
  ###### Eussocial ou Não Eussocial (A+M)x(B+E)
  #BMM - múltiplas taxas
  Eusoc_BMM_ModD = mvBM(tree2, auto_Mod_D, model = "BMM")
  #OUM com múltiplos ótimos
  Eusoc_OUM_ModD = mvOU(tree2, auto_Mod_D, model = "OUM")
  
  
  ##### Social (A+M)x(B)x(E)
  #BMM - múltiplas taxas
  Soc_BMM_ModD = mvBM(tree3, auto_Mod_D, model = "BMM")
  #OU com múltiplos ótimos
  Soc_OUM_ModD = mvOU(tree3, auto_Mod_D, model = "OUM")
  
  Sol_BMM_ModD = mvBM(tree4, auto_Mod_D, model = "BMM")
  #OUM com múltiplos ótimos
  Sol_OUM_ModD = mvOU(tree4, auto_Mod_D, model = "OUM")
  
  ##### Tamanho
  size_BMM_ModD = mvBM(tree5, auto_Mod_D, model = "BMM")
  #OUM com múltiplos ótimos
  size_OUM_ModD = mvOU(tree5, auto_Mod_D, model = "OUM")
  
  
  
  #Comparando modelos
  nomes_modelos<-c("BM","EB","OU1","Eusoc_BMM","Eusoc_OUM","Soc_BMM","Soc_OUM", "Sol_BMM","Sol_OUM","size_BMM","size_OUM") ## CUIDADO COM A ORDEM
  
  #Tabela com os resultados
  Modu_Dtes<-as.data.frame(c("Lik","AIC","AICc"))
  rownames(Modu_Dtes)<-Modu_Dtes[, 1]
  Modu_Dtes<-Modu_Dtes[,-1]
  
  Modu_Dtes$BM_ModD<-c(BM_ModD$LogLik,BM_ModD$AIC,BM_ModD$AICc)
  Modu_Dtes$EB_ModD<-c(EB_ModD$LogLik,EB_ModD$AIC,EB_ModD$AICc)
  Modu_Dtes$OU1_ModD<-c(OU1_ModD$LogLik,OU1_ModD$AIC,OU1_ModD$AICc)
  
  Modu_Dtes$Eusoc_BMM_ModD<-c(Eusoc_BMM_ModD$LogLik,Eusoc_BMM_ModD$AIC,Eusoc_BMM_ModD$AICc)
  Modu_Dtes$Eusoc_OUM_ModD<-c(Eusoc_OUM_ModD$LogLik,Eusoc_OUM_ModD$AIC,Eusoc_OUM_ModD$AICc)
  
  Modu_Dtes$Soc_BMM_ModD<-c(Soc_BMM_ModD$LogLik,Soc_BMM_ModD$AIC,Soc_BMM_ModD$AICc)
  Modu_Dtes$Soc_OUM_ModD<-c(Soc_OUM_ModD$LogLik,Soc_OUM_ModD$AIC,Soc_OUM_ModD$AICc)
  
  Modu_Dtes$Sol_BMM_ModD<-c(Sol_BMM_ModD$LogLik,Sol_BMM_ModD$AIC,Sol_BMM_ModD$AICc)
  Modu_Dtes$Sol_OUM_ModD<-c(Sol_OUM_ModD$LogLik,Sol_OUM_ModD$AIC,Sol_OUM_ModD$AICc)
  
  Modu_Dtes$size_BMM_ModD<-c(size_BMM_ModD$LogLik,size_BMM_ModD$AIC,size_BMM_ModD$AICc)
  Modu_Dtes$size_OUM_ModD<-c(size_OUM_ModD$LogLik,size_OUM_ModD$AIC,size_OUM_ModD$AICc)
  
  Modu_Dtabela<-as.data.frame(t(Modu_Dtes)) #Inverter linhas e colunas
  
  # Calcular ΔAICc, w e AICw 
  
  Mod_DAICc<-mvMORPH::aicw(Modu_Dtabela$AICc) #calculo
  Mod_DAICc$models<-nomes_modelos
  Mod_DAICc
  
  #Adicionando a tabela
  
  Modu_Dtabela$AIC<-round(Modu_Dtabela$AIC, 0)
  Modu_Dtabela$AICc<-round(Mod_DAICc$AIC, 0)
  Modu_Dtabela$ΔAICc<-round(Mod_DAICc$diff,2)
  Modu_Dtabela$w<-round(Mod_DAICc$wi,3)
  Modu_Dtabela$AICcw<-round(Mod_DAICc$aicweights,3)
  
  print(Modu_Dtabela)
  
  write.csv(Modu_Dtabela,"ModEv_Mod_D.csv") #Salvando tabela
  Mod_DAICc
  
  
  Parametros_ModD_SocOUM<- Soc_OUM_ModD$theta
  
  
}

#Modulo P #Result Soc_BMM
{ 
  #BM - uma taxa
  BM_ModP = mvBM(tree2, auto_Mod_P, model = "BM1")
  #EB
  EB_ModP = mvEB(tree2, auto_Mod_P)
  #OU
  OU1_ModP = mvOU(tree2, auto_Mod_P, model = "OU1")
  
  ###### Eussocial ou Não Eussocial (A+M)x(B+E)
  #BMM - múltiplas taxas
  Eusoc_BMM_ModP = mvBM(tree2, auto_Mod_P, model = "BMM")
  #OUM com múltiplos ótimos
  Eusoc_OUM_ModP = mvOU(tree2, auto_Mod_P, model = "OUM")
  
  
  ##### Social (A+M)x(B)x(E)
  #BMM - múltiplas taxas
  Soc_BMM_ModP = mvBM(tree3, auto_Mod_P, model = "BMM")
  #OU com múltiplos ótimos
  Soc_OUM_ModP = mvOU(tree3, auto_Mod_P, model = "OUM")
  
  #Solitário
  Sol_BMM_ModP = mvBM(tree4, auto_Mod_P, model = "BMM")
  Sol_OUM_ModP = mvOU(tree4, auto_Mod_P, model = "OUM")
  
  ##### Tamanho
  size_BMM_ModP = mvBM(tree5, auto_Mod_P, model = "BMM")
  #OUM com múltiplos ótimos
  size_OUM_ModP = mvOU(tree5, auto_Mod_P, model = "OUM")
  
  
  
  #Comparando modelos
  nomes_modelos<-c("BM","EB","OU1","Eusoc_BMM","Eusoc_OUM","Soc_BMM","Soc_OUM", "Sol_BMM","Sol_OUM","size_BMM","size_OUM") ## CUIDADO COM A ORDEM
  
  #Tabela com os resultados
  Modu_P_tes<-as.data.frame(c("Lik","AIC","AICc"))
  rownames(Modu_P_tes)<-Modu_P_tes[, 1]
  Modu_P_tes<-Modu_P_tes[,-1]
  
  Modu_P_tes$BM_ModP<-c(BM_ModP$LogLik,BM_ModP$AIC,BM_ModP$AICc)
  Modu_P_tes$EB_ModP<-c(EB_ModP$LogLik,EB_ModP$AIC,EB_ModP$AICc)
  Modu_P_tes$OU1_ModP<-c(OU1_ModP$LogLik,OU1_ModP$AIC,OU1_ModP$AICc)
  
  Modu_P_tes$Eusoc_BMM_ModP<-c(Eusoc_BMM_ModP$LogLik,Eusoc_BMM_ModP$AIC,Eusoc_BMM_ModP$AICc)
  Modu_P_tes$Eusoc_OUM_ModP<-c(Eusoc_OUM_ModP$LogLik,Eusoc_OUM_ModP$AIC,Eusoc_OUM_ModP$AICc)
  
  Modu_P_tes$Soc_BMM_ModP<-c(Soc_BMM_ModP$LogLik,Soc_BMM_ModP$AIC,Soc_BMM_ModP$AICc)
  Modu_P_tes$Soc_OUM_ModP<-c(Soc_OUM_ModP$LogLik,Soc_OUM_ModP$AIC,Soc_OUM_ModP$AICc)
  
  Modu_P_tes$Sol_BMM_ModP<-c(Sol_BMM_ModP$LogLik,Sol_BMM_ModP$AIC,Sol_BMM_ModP$AICc)
  Modu_P_tes$Sol_OUM_ModP<-c(Sol_OUM_ModP$LogLik,Sol_OUM_ModP$AIC,Sol_OUM_ModP$AICc)
  
  Modu_P_tes$size_BMM_ModP<-c(size_BMM_ModP$LogLik,size_BMM_ModP$AIC,size_BMM_ModP$AICc)
  Modu_P_tes$size_OUM_ModP<-c(size_OUM_ModP$LogLik,size_OUM_ModP$AIC,size_OUM_ModP$AICc)
  
  Modu_P_tabela<-as.data.frame(t(Modu_P_tes)) #Inverter linhas e colunas
  
  # Calcular ΔAICc, w e AICw 
  
  Mod_PAICc<-mvMORPH::aicw(Modu_P_tabela$AICc) #calculo
  Mod_PAICc$models<-nomes_modelos
  Mod_PAICc
  
  #Adicionando a tabela
  
  Modu_P_tabela$AIC<-round(Modu_P_tabela$AIC, 0)
  Modu_P_tabela$AICc<-round(Mod_PAICc$AIC, 0)
  Modu_P_tabela$ΔAICc<-round(Mod_PAICc$diff,2)
  Modu_P_tabela$w<-round(Mod_PAICc$wi,3)
  Modu_P_tabela$AICcw<-round(Mod_PAICc$aicweights,3)
  
  print(Modu_P_tabela)
  
  write.csv(Modu_P_tabela,"ModEv_Mod_P.csv") #Salvando tabela
  Mod_PAICc
  
  Parametros_ModP_SocBMM<- Soc_BMM_ModP$theta
  
}

#Modulo S #Result Size_BMM
{
  #BM - uma taxa
  BM_ModS = mvBM(tree2, auto_Mod_S, model = "BM1")
  #EB
  EB_ModS = mvEB(tree2, auto_Mod_S)
  #OU
  OU1_ModS = mvOU(tree2, auto_Mod_S, model = "OU1")
  
  ###### Eussocial ou Não Eussocial (A+M)x(B+E)
  #BMM - múltiplas taxas
  Eusoc_BMM_ModS = mvBM(tree2, auto_Mod_S, model = "BMM")
  #OUM com múltiplos ótimos
  Eusoc_OUM_ModS = mvOU(tree2, auto_Mod_S, model = "OUM")
  
  
  ##### Social (A+M)x(B)x(E)
  #BMM - múltiplas taxas
  Soc_BMM_ModS = mvBM(tree3, auto_Mod_S, model = "BMM")
  #OU com múltiplos ótimos
  Soc_OUM_ModS = mvOU(tree3, auto_Mod_S, model = "OUM")
  
  Sol_BMM_ModS = mvBM(tree4, auto_Mod_S, model = "BMM")
  #OUM com múltiplos ótimos
  Sol_OUM_ModS = mvOU(tree4, auto_Mod_S, model = "OUM")
  
  ##### Tamanho
  size_BMM_ModS = mvBM(tree5, auto_Mod_S, model = "BMM")
  #OUM com múltiplos ótimos
  size_OUM_ModS = mvOU(tree5, auto_Mod_S, model = "OUM")
  
  
  
  #Comparando modelos
  nomes_modelos<-c("BM","EB","OU1","Eusoc_BMM","Eusoc_OUM","Soc_BMM","Soc_OUM", "Sol_BMM","Sol_OUM","size_BMM","size_OUM") ## CUIDADO COM A ORDEM
  
  #Tabela com os resultados
  Modulo_Stes<-as.data.frame(c("Lik","AIC","AICc"))
  rownames(Modulo_Stes)<-Modulo_Stes[, 1]
  Modulo_Stes<-Modulo_Stes[,-1]
  
  Modulo_Stes$BM_ModS<-c(BM_ModS$LogLik,BM_ModS$AIC,BM_ModS$AICc)
  Modulo_Stes$EB_ModS<-c(EB_ModS$LogLik,EB_ModS$AIC,EB_ModS$AICc)
  Modulo_Stes$OU1_ModS<-c(OU1_ModS$LogLik,OU1_ModS$AIC,OU1_ModS$AICc)
  
  Modulo_Stes$Eusoc_BMM_ModS<-c(Eusoc_BMM_ModS$LogLik,Eusoc_BMM_ModS$AIC,Eusoc_BMM_ModS$AICc)
  Modulo_Stes$Eusoc_OUM_ModS<-c(Eusoc_OUM_ModS$LogLik,Eusoc_OUM_ModS$AIC,Eusoc_OUM_ModS$AICc)
  
  Modulo_Stes$Soc_BMM_ModS<-c(Soc_BMM_ModS$LogLik,Soc_BMM_ModS$AIC,Soc_BMM_ModS$AICc)
  Modulo_Stes$Soc_OUM_ModS<-c(Soc_OUM_ModS$LogLik,Soc_OUM_ModS$AIC,Soc_OUM_ModS$AICc)
  
  Modulo_Stes$Sol_BMM_ModS<-c(Sol_BMM_ModS$LogLik,Sol_BMM_ModS$AIC,Sol_BMM_ModS$AICc)
  Modulo_Stes$Sol_OUM_ModS<-c(Sol_OUM_ModS$LogLik,Sol_OUM_ModS$AIC,Sol_OUM_ModS$AICc)
  
  Modulo_Stes$size_BMM_ModS<-c(size_BMM_ModS$LogLik,size_BMM_ModS$AIC,size_BMM_ModS$AICc)
  Modulo_Stes$size_OUM_ModS<-c(size_OUM_ModS$LogLik,size_OUM_ModS$AIC,size_OUM_ModS$AICc)
  
  Modulo_Stabela<-as.data.frame(t(Modulo_Stes)) #Inverter linhas e colunas
  
  # Calcular ΔAICc, w e AICw 
  
  Mod_SAICc<-mvMORPH::aicw(Modulo_Stabela$AICc) #calculo
  Mod_SAICc$models<-nomes_modelos
  Mod_SAICc
  
  #Adicionando a tabela
  
  Modulo_Stabela$AIC<-round(Modulo_Stabela$AIC, 0)
  Modulo_Stabela$AICc<-round(Mod_SAICc$AIC, 0)
  Modulo_Stabela$ΔAICc<-round(Mod_SAICc$diff,2)
  Modulo_Stabela$w<-round(Mod_SAICc$wi,3)
  Modulo_Stabela$AICcw<-round(Mod_SAICc$aicweights,3)
  
  print(Modulo_Stabela)
  
  write.csv(Modulo_Stabela,"ModEv_Mod_S.csv") #Salvando tabela
  Mod_SAICc
  
  Parametros_ModS_SizeBMM<- size_BMM_ModS$theta
}



}

# Parametros dos modelos dos modulos
Parametros_ModA_SocOUM #Theta Modulo A OUM
Parametros_ModD_SocOUM #Theta Modulo D OUM
Parametros_ModP_SocBMM #Theta Modulo P BMM
Parametros_ModS_SizeBMM ##Theta Modulo S BMM

#Plots dos modulos
#Paramatros modulos

par(mfrow=c(2,2))

#Modulo A
plot(PCA_Mod_A, flip = 1, phylo = TRUE, main = "PCA_AM", phylo.par = phy.par,
     pch = 21, bg = col.gp, cex = 1.3, cex.lab = 1.5, font.lab = 2,
     xlim = c(-0.11, 0.1),   # limites para PC1
     ylim = c(-0.04, 0.06))  # limites para PC2
points(cbind(-Parametros_ModA_SocOUM[,1],  # inverte eixo 1 (X)
             Parametros_ModA_SocOUM[,2]),
       pch = 23, bg = "red", cex = 1.8)
text(x = -Parametros_ModA_SocOUM[,1],  # inverte X
     y = Parametros_ModA_SocOUM[,2],
     labels = rownames(Parametros_ModA_SocOUM),
     pos = 3, cex = 1.1, font = 2)
#3D
#plot3d(PCA_Mod_A$x[,1:3], main = "PCA_AM", col = col.gp, size = 10)
#points3d(Parametros_ModA_SocOUM[,1:3], col = "red", size = 25, shape = "tetra",)
#text3d(Parametros_ModA_SocOUM[,1:3],texts = rownames(Parametros_ModA_SocOUM), adj = 1.2, cex = 1.1)

#Modulo D
plot(PCA_Mod_D, flip = 1, phylo = TRUE, main = "PCA_DM", phylo.par = phy.par,
     pch = 21, bg = col.gp, cex = 1.3, cex.lab = 1.5, font.lab = 2)
points(cbind(-Parametros_ModD_SocOUM[,1],  # inverte eixo 1 (X)
             Parametros_ModD_SocOUM[,2]),
       pch = 23, bg = "red", cex = 1.8)
text(x = -Parametros_ModD_SocOUM[,1],  # inverte X
     y = Parametros_ModD_SocOUM[,2],
     labels = rownames(Parametros_ModD_SocOUM),
     pos = 3, cex = 1.1, font = 2)
#3D
#plot3d(PCA_Mod_D$x[,1:3], main = "PCA_DM", col = col.gp, size = 10)
#points3d(Parametros_ModD_SocOUM[,1:3], col = "red", size = 25, shape = "tetra",)
#text3d(Parametros_ModD_SocOUM[,1:3],texts = rownames(Parametros_ModD_SocOUM), adj = 1.2, cex = 1.1)

#Modulo P NÃO APARECE
plot(PCA_Mod_P,flip = 1, phylo = TRUE, main = "PCA_PM", phylo.par = phy.par,
     pch = 21, bg = col.gp, cex = 1.3, cex.lab = 1.5, font.lab = 2)
points(Parametros_ModP_SocBMM[,1:3],  pch = 23, bg = "red", cex = 1.8)
text(Parametros_ModP_SocBMM[,1:3], labels = rownames(Parametros_ModP_SocBMM[,1:3]), pos = 3, cex = 1.1, font = 2)
#3D
#plot3d(PCA_Mod_P$x[,1:3], main = "PCA_PM", col = col.gp, size = 10)
#points3d(Parametros_ModP_SocBMM[,1:3], col = "red", size = 25, shape = "tetra",)
#text3d(Parametros_ModP_SocBMM[,1:3],texts = "BMM_PM", adj = 1.2, cex = 1.1)

#Modulo S NÃO APARECE
plot(PCA_Mod_S, phylo = TRUE, main = "PCA_SM", phylo.par = phy.par,
     pch = 21, bg = col.gp, cex = 1.3, cex.lab = 1.5, font.lab = 2)
points(Parametros_ModS_SizeBMM[,1:3],  pch = 23, bg = "red", cex = 1.8)
text(Parametros_ModS_SizeBMM[,1:3], labels = rownames(Parametros_ModS_SizeBMM[,1:3]), pos = 3, cex = 1.1, font = 2)
#3D
#plot3d(PCA_Mod_S$x[,1:3], main = "PCA_SM", col = col.gp, size = 10)
#points3d(Parametros_ModS_SizeBMM[,1:3], col = "pink", size = 25, shape = "tetra",)
#text3d(Parametros_ModS_SizeBMM[,1:3],texts = "BMM_SM", adj = 1.2, cex = 1.1)

par(mfrow=c(1,1))

################################################################
 # DIsparidade por módulo

{

# Extrair os dois primeiros componentes
autoval <- PCA_Mphy$x[, 1:3]
View(autoval)

################################ Disparidade ##########################
average.sq <- function(X) mean(pairwise.dist(X)^2)

dispRity_dtt <- dtt.dispRity(data = autoval, metric = average.sq, tree = Phy_tree, nsim = 200)

plot(dispRity_dtt, xlim = rev(range( dispRity_dtt$times)))

#DTT dos modulos
autoval_Mod_A<-PCA_Mod_A$x[, 1:3]
autoval_Mod_D<-PCA_Mod_D$x[, 1:3]
autoval_Mod_P<-PCA_Mod_P$x[, 1:3]
autoval_Mod_S<-PCA_Mod_S$x[, 1:3]

#Modulo A
dispRity_dtt_Mod_A <- dtt.dispRity(data = autoval_Mod_A, metric = average.sq, tree = Phy_tree, nsim = 200)
#Modulo D
dispRity_dtt_Mod_D <- dtt.dispRity(data = autoval_Mod_D, metric = average.sq, tree = Phy_tree, nsim = 200)
#Modulo P
dispRity_dtt_Mod_P <- dtt.dispRity(data = autoval_Mod_P, metric = average.sq, tree = Phy_tree, nsim = 200)
#Modulo S
dispRity_dtt_Mod_S <- dtt.dispRity(data = autoval_Mod_S, metric = average.sq, tree = Phy_tree, nsim = 200)

#PLOTS
par(mfrow=c(2,2))
plot(dispRity_dtt_Mod_A, main = "Module AM", xlim = rev(range( dispRity_dtt_Mod_A$times)))
plot(dispRity_dtt_Mod_D, main = "Module DM",xlim = rev(range( dispRity_dtt_Mod_D$times)))
plot(dispRity_dtt_Mod_P, main = "Module PM",xlim = rev(range( dispRity_dtt_Mod_P$times)))
plot(dispRity_dtt_Mod_S, main = "Module SM",xlim = rev(range( dispRity_dtt_Mod_S$times)))
par(mfrow=c(1,1))

dispRity_dtt_Mod_A$MDI
dispRity_dtt_Mod_D$MDI
dispRity_dtt_Mod_P$MDI
dispRity_dtt_Mod_S$MDI
}
