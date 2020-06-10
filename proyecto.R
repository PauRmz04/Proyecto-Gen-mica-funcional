###MURCIELAGOS
#Los datos fueron obtenidos del artículo "Network analysis of host-virus communities in bats and rodents reveals determinants of cross-species transmission"
#Se cargan los datos
batvirusdata = read.csv("ele12491-sup-0003-batvirus data.csv")
battraits = read.csv("ele12491-sup-0005-battraitdata.csv")


battraits$cave=factor(battraits$cave)
battraits$migration=factor(battraits$migration)
battraits$torpor.category=factor(battraits$torpor.category)

batnames=sort(unique(batvirusdata$binomial))



bat.host.virus.list=list(NA) 
for(i in 1:length(batnames)){
  bat.host.virus.list[[i]]= batvirusdata $virus.species[which(batvirusdata $binomial==batnames[i])]
}
names(bat.host.virus.list)=batnames	


##Se cre una matriz de adjacencia
batsharedvirusesmat=matrix(NA,ncol=length(batnames),nrow=length(batnames))
colnames(batsharedvirusesmat)= batnames
rownames(batsharedvirusesmat)= batnames

for (i in 1:length(batnames)){
  for (j in 1:length(batnames)){	 
    batsharedvirusesmat[i,j]=length(which(is.finite(match(unique(bat.host.virus.list[[i]]),unique(bat.host.virus.list[[j]])))))
  }
}


##Librería utilizada
library(igraph)

batsharedvirusesgraph=graph.adjacency(batsharedvirusesmat,mode="undirected",weighted=TRUE,diag=FALSE)
batsharedvirusesgraph
summary(batsharedvirusesgraph) 
plot(batsharedvirusesgraph, vertex.color = "chartreuse2", vertex.size = 20, vertex.label.color= "chartreuse4")

##Propiedades de la red
class(batsharedvirusesgraph)


adj <- get.adjacency(batsharedvirusesgraph, attr="weight")
adj 

edge_density(batsharedvirusesgraph)
reciprocity(batsharedvirusesgraph)
transitivity(batsharedvirusesgraph, type="global")

diameter(batsharedvirusesgraph)
degree(batsharedvirusesgraph)
hist(degree(batsharedvirusesgraph), col= "orange")
degree.distribution(batsharedvirusesgraph)
plot(degree.distribution(batsharedvirusesgraph), col="orangered")


##Centralidad
degree(batsharedvirusesgraph)
centr_degree(batsharedvirusesgraph)
closeness(batsharedvirusesgraph)
centr_clo(batsharedvirusesgraph)
eigen_centrality(batsharedvirusesgraph)
centr_eigen(batsharedvirusesgraph)
betweenness(batsharedvirusesgraph)
edge_betweenness(batsharedvirusesgraph)
centr_betw(batsharedvirusesgraph)


hub.score(batsharedvirusesgraph)

##Clusterización
 
ceb <- cluster_edge_betweenness(batsharedvirusesgraph)
dendPlot(ceb, mode="hclust") 

## ##ahora para hacer lo grafico  
plot(ceb, batsharedvirusesgraph) 
##
clp <- cluster_label_prop(batsharedvirusesgraph)  
plot(clp, batsharedvirusesgraph) 
##
cfg <- cluster_fast_greedy(as.undirected(batsharedvirusesgraph))  
plot(cfg, as.undirected(batsharedvirusesgraph))

modularity(batsharedvirusesgraph, membership(ceb))
modularity_matrix(batsharedvirusesgraph, membership(ceb))




###ROEDORES
rodentvirusdata = read.csv("ele12491-sup-0004-rodentvirusdata.csv")
rodtraits = read.csv("ele12491-sup-0006-rodenttraitdata.csv")


rodtraits$cave=factor(rodtraits$cave)
rodtraits$migration=factor(rodtraits$migration)
rodtraits$torpor.category=factor(rodtraits$torpor.category)

rodnames=sort(unique(rodentvirusdata$binomial))



rod.host.virus.list=list(NA) 
for(i in 1:length(rodnames)){
  rod.host.virus.list[[i]]= rodentvirusdata $virus.species[which(rodentvirusdata $binomial==rodnames[i])]
}
names(rod.host.virus.list)=rodnames	



rodsharedvirusesmat=matrix(NA,ncol=length(rodnames),nrow=length(rodnames))
colnames(rodsharedvirusesmat)= rodnames
rownames(rodsharedvirusesmat)= rodnames

for (i in 1:length(rodnames)){
  for (j in 1:length(rodnames)){	 
    rodsharedvirusesmat[i,j]=length(which(is.finite(match(unique(rod.host.virus.list[[i]]),unique(rod.host.virus.list[[j]])))))
  }
}


library(igraph)

rodsharedvirusesgraph=graph.adjacency(rodsharedvirusesmat,mode="undirected",weighted=TRUE,diag=FALSE)
rodsharedvirusesgraph
summary(rodsharedvirusesgraph) 
plot(rodsharedvirusesgraph, vertex.color= "darkseagreen", vertex.size=20, vertex.label.color="darksalmon")


class(rodsharedvirusesgraph)

adjr<-get.adjacency(rodsharedvirusesgraph, attr="weight")
adjr

##
##
deg <- degree(rodsharedvirusesgraph, mode="all")
deg

plot(rodsharedvirusesgraph, vertex.size=deg)
##Histograma con el degree
hist(deg, main="Histogram of node degree the rodents", col = "turquoise3")

##degree distribution plot 
deg.dist <- degree_distribution(rodsharedvirusesgraph, cumulative=T, mode="all")

plot(deg.dist, col = "palevioletred4")
##
hs <- hub_score(rodsharedvirusesgraph, weights=NA)$vector


diameter(rodsharedvirusesgraph)
edge_density(rodsharedvirusesgraph)
reciprocity(rodsharedvirusesgraph)
transitivity(rodsharedvirusesgraph)
degree.distribution(rodsharedvirusesgraph)
betweenness(rodsharedvirusesgraph,directed=FALSE)


#para hacer cluster con la de roedores 
ceb <- cluster_edge_betweenness(rodsharedvirusesgraph)
dendPlot(ceb, mode="hclust") 
## 
plot(ceb, rodsharedvirusesgraph) 
##
clp <- cluster_label_prop(rodsharedvirusesgraph)  
plot(clp, rodsharedvirusesgraph) 
##
cfg <- cluster_fast_greedy(as.undirected(rodsharedvirusesgraph))  
plot(cfg, as.undirected(rodsharedvirusesgraph))

modularity(rodsharedvirusesgraph, membership(ceb))
modularity_matrix(rodsharedvirusesgraph, membership(ceb))


