library(lifecycle)

sim.bdh.age <-function(age,numbsim,
                      lambda,mu,
                      nu, hybprops,
                      hyb.inher.fxn,
                      frac=1,twolineages=FALSE,complete=TRUE,stochsampling=FALSE,
                      hyb.rate.fxn=NULL,
                      trait.model=NULL,
                      mrca=deprecated(),Ngene){ #Ngene how many genes to simulation, if Ngene=0, all possible topologies will be listed
  if (is_present(mrca)) {
    lifecycle::deprecate_warn("1.1.0", "sim.bdh.age(mrca)", "sim.bdh.age(twolineages)")
    twolineages <- mrca
  }
    out<-lapply(1:numbsim,sim.bdh.age.help2,
                age=age,
                lambda=lambda,mu=mu,
                nu=nu, hybprops=hybprops,
                hyb.rate.fxn=hyb.rate.fxn,
                hyb.inher.fxn=hyb.inher.fxn,
                frac=frac,mrca=twolineages,
                complete=complete, stochsampling=stochsampling,
                trait.model=trait.model,Ngene=Ngene)
    class(out)<-c('list','multiPhylo')
    out
}
sim.bdh.age.help2 <-function(dummy,age,lambda,mu,
                           nu, hybprops,hyb.rate.fxn,
                           hyb.inher.fxn,
                           frac=1,mrca=FALSE,complete=TRUE,stochsampling=FALSE,
                           trait.model,Ngene){
	out<-sim.bdh.age.loop2(age=age,
	                      lambda=lambda,mu=mu,
	                      nu=nu, hybprops=hybprops,
	                      hyb.inher.fxn=hyb.inher.fxn,
	                      hyb.rate.fxn=hyb.rate.fxn,
	                      frac=frac,mrca=mrca,
	                      complete=complete,stochsampling=stochsampling,
	                      trait.model=trait.model,Ngene=Ngene)
	out
}

sim.bdh.age.loop2 <- function(age,
                             lambda,mu,
                             nu,hybprops,
                             hyb.inher.fxn,
                             hyb.rate.fxn,
                             frac=1,mrca,complete,stochsampling,
                             trait.model,Ngene) {
    phy <- sim2.bdh.origin2(m=0,n=0,
                           age=age,
                           lambda=lambda,mu=mu,
                           nu=nu,
                           hyb.inher.fxn=hyb.inher.fxn,
                           hybprops=hybprops,hyb.rate.fxn=hyb.rate.fxn,
                           mrca=mrca,
                           trait.model=trait.model,Ngene=Ngene)

    if( "phylo" %in% class(phy)){
      nleaves<-phy$nleaves
      phy$nleaves<-NULL

      ###Determine how many tips need to be deleted###
      if(stochsampling){
        ntips_remaining<-rbinom(1,nleaves,frac)
      }else{
        ntips_remaining<-round(nleaves*frac)
      }
      phy<-SiPhyNetwork:::handleTipsTaxa(phy=phy,complete=complete,target_ntaxa = ntips_remaining, current_n = nleaves)


    }
    phy
  }

sim2.bdh.origin2 <- function(m=0,n=0,age,lambda,mu,nu,hyb.inher.fxn,hybprops,hyb.rate.fxn,mrca,trait.model,Ngene){

    ##TODO make edge a list and convert to a matrix at the end for quicker building
    hyb_edge<-list()    #list of hybrid edges
    num_hybs<-0         #number of hybridizations
    inheritance<- c()   #Inheritance probabilities of hyb edges
    time <-0				#time after origin
    extinct <- vector()		#list of extinct leaves
    extincttree = 0
    stop = 0
    time_in_n<-c()##records the times we have n number of lineages

    is_null_trait<- is.null(trait.model)
    if(!is_null_trait){

      if( !all(c("initial","hyb.event.fxn","time.fxn","spec.fxn","hyb.compatibility.fxn") %in% names(trait.model))){
        stop('the trait model does not have all the needed named components')
      }



      trait_state <- trait.model[['initial']]
      ext_trait_state<-c()
      if(length(trait_state)!= (mrca+1)){
        stop(paste('there are ',(mrca+1),' starting species. There were ',length(trait_state),' initial polyploid states',sep = ''))
      }
    }

    if(!mrca){#Start with one lineage
      edge <- matrix(c(-1,-2),nrow=1,ncol=2)		#matrix of edges
      leaves <- c(-2)			#list of extant leaves
      nleaves<-1
      timecreation <-c(0,0)		#time when species -2, -3, ... -n was created after origin
      maxspecies <- -2		#smallest species
      edge.length <- c(0)		#edge length. if 0: leaf which didn't speciate /extinct yet
      genetic_dists<-list('-2'=list('-2'=0.0)) ##list of genetic distances between species
      genetic_dists2<-list('-1'=list('-1'=0.0,'-2'=0.0),'-2'=list('-1'=0.0,'-2'=0.0))
      tip <- "-2"
    }else{#start with two lineages
      edge<-rbind(c(-1,-2),c(-1,-3))
      leaves<-c(-2,-3)
      nleaves<-2
      timecreation<-c(0,0,0)
      maxspecies<- -3
      edge.length<- c(0,0)
      genetic_dists<-list('-2'=list('-2'=0.0,'-3'=0.0),'-3'=list('-2'=0.0,'-3'=0.0))
      genetic_dists<-list('-1'=list('-1'=0.0,'-2'=0.0,'-3'=0.0),'-2'=list('-1'=0.0,'-2'=0.0,'-3'=0.0),'-3'=list('-1'=0.0,'-2'=0.0,'-3'=0.0))
      nleaves<-2
      tip <- c("-2","-3")
    }
    edge.length2 <- edge.length #edge.length2 is the tree that over estimate pairwise distance, while minimize the difference
    node <- "-1"
    hyb_edge.length <- NULL
    #make a list of trees
    trees <- list()
    trees[[1]] <- list(edge=edge)
    trees.prob <- 1
    if (Ngene>0) {
        trees <- rep(trees,Ngene)
    }
    while (stop == 0 ){
      if (nleaves == 0){
        #phy2 = 0
        if (age >0) {
          phy=0
        }
        extincttree=1
        stop = 1
      } else {
        spec_rate<-nleaves*lambda
        ext_rate<-nleaves*mu
        hyb_rate<-choose(nleaves,2)*nu
        total_rate<-spec_rate+ext_rate+hyb_rate
        timestep <- rexp(1, (total_rate) )    #time since last event
        if ((age > 0 && (time+timestep) < age) && !(nleaves==m) ) {
          time = time+timestep			#time after origin

          if(nleaves==n){
            time_in_n<-c(time_in_n,time-timestep,time)
          }

          genetic_dists<-lapply(genetic_dists, function(x) lapply(x, function(x) x+(2*timestep))) ##all species increase genetic distance by twice the timestep
          for(i in names(genetic_dists)){ ##lineages shouldn't have a distance from themselves
            genetic_dists[[i]][[i]]<- genetic_dists[[i]][[i]]-(2*timestep)
          }
          if (!is.null(node)) {
          for(i in node){ 
          	for (j in tip) {
            	genetic_dists2[[i]][[j]]<- genetic_dists2[[i]][[j]]+timestep
            	genetic_dists2[[j]][[i]]<- genetic_dists2[[j]][[i]]+timestep
          	}
          }
          }
          for(i in tip){ 
          	for (j in tip) {
          		if (i!=j) {
            	genetic_dists2[[i]][[j]]<- genetic_dists2[[i]][[j]]+2*timestep
            	}
          	}
          }
          
          if(!is_null_trait){
            trait_state<-trait.model[['time.fxn']](trait_state,timestep)
          }

          species <- sample(leaves,1+(nleaves>1))
          del <- match(species,leaves)
          edgespecevent <- match(species,edge) - length(edge.length)

          randevent <- runif(1,0,1)		#event speciation or extinction or hybridization
          if ( randevent <= (spec_rate/total_rate) ) {##Speciation Event
            ##only use the first species
            species<-species[1]
            del<-del[1]
            edgespecevent<-edgespecevent[1]

            edge.length[edgespecevent] <- time-timecreation[- species]
            edge <- rbind(edge,c(species,maxspecies-1))
            edge <- rbind(edge,c(species,maxspecies-2))
            edge.length <- c(edge.length,0,0)
            edge.length2 <- c(edge.length2,0,0)
            leaves <- c(leaves,maxspecies-1,maxspecies-2)
            
            leaves <- leaves[- del]
            timecreation <- c(timecreation,time,time)
            timecreation[-species]<-time
            
            #update tree list
            for (i in 1:length(trees)) {
                trees[[i]]$edge <- rbind(trees[[i]]$edge,c(species,maxspecies-1))
                trees[[i]]$edge <- rbind(trees[[i]]$edge,c(species,maxspecies-2))
            }
            
            #if(!is.null(hyb.rate.fxn)){
              ##update genetic distances
              new_species_rw<-genetic_dists[[as.character(species)]]
              new_species_rw[[as.character(species)]]<-NULL
              new_species_rw[[as.character(maxspecies-1)]]<-0.0
              new_species_rw[[as.character(maxspecies-2)]]<-0.0
              node <- c(node,as.character(species)) #add it to node
              tip <- tip[-which(tip==as.character(species))] #remove it from tip
              tip <- c(tip,as.character(maxspecies-1),as.character(maxspecies-2)) #add new species to tip
              genetic_dists[[as.character(species)]]<-NULL ##remove the row with species
              for(i in names(genetic_dists)){ ##add column for maxspecies-1 and maxspecies-2. Also delete column with species
                genetic_dists[[i]][[as.character(maxspecies-1)]]<-new_species_rw[[i]]
                genetic_dists[[i]][[as.character(maxspecies-2)]]<-new_species_rw[[i]]
                genetic_dists[[i]][[as.character(species)]]<-NULL
              }
              genetic_dists[[as.character(maxspecies-1)]]<-new_species_rw ##add row for maxspecies-1
              genetic_dists[[as.character(maxspecies-2)]]<-new_species_rw ##add row for maxspecies-2
              new_species_rw[[as.character(species)]]<-0.0
              
              new_species_rw<-genetic_dists2[[as.character(species)]]
              new_species_rw[[as.character(maxspecies-1)]]<-0.0
              new_species_rw[[as.character(maxspecies-2)]]<-0.0
              for(i in names(genetic_dists2)){ ##add column for maxspecies-1 and maxspecies-2.
                genetic_dists2[[i]][[as.character(maxspecies-1)]]<-new_species_rw[[i]]
                genetic_dists2[[i]][[as.character(maxspecies-2)]]<-new_species_rw[[i]]
              }
              genetic_dists2[[as.character(maxspecies-1)]]<-new_species_rw ##add row for maxspecies-1
              genetic_dists2[[as.character(maxspecies-2)]]<-new_species_rw ##add row for maxspecies-2
            #}
            if(!is_null_trait){
              trait_state<-c( trait_state, trait.model[['spec.fxn']](trait_state[del]) )
              trait_state<-trait_state[-del]
            }

            maxspecies <- maxspecies-2
            nleaves<-nleaves+1
          } else if( randevent <= ((spec_rate+ext_rate)/total_rate) ){##Extinction Event
            ##only use the first species
            species<-species[1]
            del<-del[1]
            edgespecevent<-edgespecevent[1]

            ##update leaves and edges
            extinct <- c(extinct,leaves[del])
            leaves <- leaves[- del]

            ##update edge lengths
            edge.length[edgespecevent] <- time-timecreation[- species]
            edge.length2[edgespecevent] <- time-timecreation[- species]
            
            timecreation[-species]<-time

			tip <- tip[-which(tip==species)]
            #if(!is.null(hyb.rate.fxn)){
              ##update genetic distances
              genetic_dists[[as.character(species)]]<-NULL ##remove the row with species
              genetic_dists2[[as.character(species)]]<-NULL
              for(i in names(genetic_dists)){ ##delete column with species
                genetic_dists[[i]][[as.character(species)]]<-NULL
              }
              for(i in names(genetic_dists2)){ ##delete column with species
                genetic_dists2[[i]][[as.character(species)]]<-NULL
              }
            #}
            if(!is_null_trait){
              ext_trait_state<-c(ext_trait_state,trait_state[del])
              trait_state<-trait_state[-del]
            }
            nleaves<-nleaves-1
          } else{ ##Hybridization Event
            hyb_occurs<-T
            inher<-hyb.inher.fxn()

            if (!is_null_trait){ ##use this for hybridization of similar ploidy
              hyb_trait <-trait.model[['hyb.event.fxn']](trait_state[del],inher)
              hyb_occurs<-trait.model[['hyb.compatibility.fxn']](trait_state[del],hyb_trait)
            }

            if (hyb_occurs && !is.null(hyb.rate.fxn)){ ##Use this to restrict hybridizations based on genetic distance
              randevent<-runif(1,0,1)
              if(randevent > hyb.rate.fxn(genetic_dists[[as.character(species[1])]][[as.character(species[2])]])){
                hyb_occurs<-F
              }
            }
            if(hyb_occurs){
              ##Update things that are constant between all hybridization types
              num_hybs<-num_hybs+1
              inheritance[num_hybs]<- inher

              randevent <- runif(1,0,1)
              if( randevent<=(hybprops[1]/sum(hybprops)) ){ #Lineage Generating
                ##update leaves
                leaves <- c(leaves,(maxspecies- 2:4))
                leaves <- leaves[- del]
                ##Update the edge matrices
                parent <- ifelse((1-inher)>0.5,1,2) ##which parent has higher inheritance
                edge <- rbind(edge,c(species[1],maxspecies-2))  ##Create new leaf for a parent lineage
                edge <- rbind(edge,c(species[2],maxspecies-3)) ##Create a new leaf for the other parent lineage
                #edge <- rbind(edge,c(species[1],maxspecies-1)) ##link the parent lineage to the hybrid node
                edge <- rbind(edge,c(species[parent],maxspecies-1)) ##link the parent lineage with higher inheritance to the hybrid node
                edge <- rbind(edge,c(maxspecies-1,maxspecies-4)) ##create a leaf for the hybrid lineage
                #hyb_edge[[num_hybs]]<-c(species[2],maxspecies-1) ##link the other parent lineage to the hybrid node
                hyb_edge[[num_hybs]]<-c(species[-parent],maxspecies-1) ##link the parent lineage with lower inheritance  to the hybrid node as reticulate edge
                ##Update edge lengths
                edge.length[edgespecevent] <- time-timecreation[- species] ##update the edge length of the lineages with the event
                edge.length <- c(edge.length,0,0,0,0)
                edge.length2[edgespecevent] <- time-timecreation[- species] ##update the edge length of the lineages with the event
                edge.length2 <- c(edge.length2,0,0,0,0)
                
                #update tree list
                if (Ngene==0) {
                trees2 <- trees
                for (i in 1:length(trees)) {
                    trees[[i]]$edge <- rbind(trees[[i]]$edge,c(species[1],maxspecies-2))
                    trees[[i]]$edge <- rbind(trees[[i]]$edge,c(species[2],maxspecies-3))
                    trees[[i]]$edge <- rbind(trees[[i]]$edge,c(species[parent],maxspecies-1))
                    trees[[i]]$edge <- rbind(trees[[i]]$edge,c(maxspecies-1,maxspecies-4))
                }
                trees.prob <- trees.prob * inher #proportion of genes from species[parent]
                for (i in 1:length(trees2)) {
                    trees2[[i]]$edge <- rbind(trees2[[i]]$edge,c(species[1],maxspecies-2))
                    trees2[[i]]$edge <- rbind(trees2[[i]]$edge,c(species[2],maxspecies-3))
                    trees2[[i]]$edge <- rbind(trees2[[i]]$edge,c(species[-parent],maxspecies-1))
                    trees2[[i]]$edge <- rbind(trees2[[i]]$edge,c(maxspecies-1,maxspecies-4))
                }
                trees <- c(trees,trees2)
                trees.prob2 <- trees.prob * (1-inher) #proportion of genes from species[-parent]
                trees.prob <- c(trees.prob,trees.prob2)
                }
                if (Ngene>0) {
                    it <- sample.int(2,Ngene,replace=T,prob=c(inher,1-inher))
                    for (i in 1:length(trees)) {
                        if (it[i]==1) {
                            trees[[i]]$edge <- rbind(trees[[i]]$edge,c(species[1],maxspecies-2))
                            trees[[i]]$edge <- rbind(trees[[i]]$edge,c(species[2],maxspecies-3))
                            trees[[i]]$edge <- rbind(trees[[i]]$edge,c(species[parent],maxspecies-1))
                            trees[[i]]$edge <- rbind(trees[[i]]$edge,c(maxspecies-1,maxspecies-4))
                        } else {
                            trees[[i]]$edge <- rbind(trees[[i]]$edge,c(species[1],maxspecies-2))
                            trees[[i]]$edge <- rbind(trees[[i]]$edge,c(species[2],maxspecies-3))
                            trees[[i]]$edge <- rbind(trees[[i]]$edge,c(species[-parent],maxspecies-1))
                            trees[[i]]$edge <- rbind(trees[[i]]$edge,c(maxspecies-1,maxspecies-4))
                        }
                    }
                }
                #if(!is.null(hyb.rate.fxn)){
                  ##update genetic distances
                  species1<-as.character(species[1])
                  species2<-as.character(species[2])
                  species3<-as.character(maxspecies-1) #hybrid node
                  maxspecies2<-as.character(maxspecies-2) #lineage of one parent
                  maxspecies3<-as.character(maxspecies-3) #lineage of the other parent
                  maxspecies4<-as.character(maxspecies-4) #hybrid lineage
                  node <- c(node, species1, species2, species3) #add the hybridization lienages to create three internal nodes
                  tip <- tip[-which(tip==species1)] #remove the hybridization lineages from tip
                  tip <- tip[-which(tip==species2)]
                  tip <- c(tip,maxspecies2,maxspecies3,maxspecies4) #add renamed hybridization lienages to tip
                  ##First we need to create the hybrid lineage
                  hybrid_species_row<-list()
                  for(i in names(genetic_dists)){
                    hybrid_species_row[[i]]<- ( (1-inher)*genetic_dists[[species1]][[i]] ) + ( inher*genetic_dists[[species2]][[i]] ) ##add row values for the recipient
                    genetic_dists[[i]][[maxspecies4]]<-hybrid_species_row[[i]] ##add column values for the recipient
                  }
                  genetic_dists[[maxspecies4]]<-hybrid_species_row
                  genetic_dists[[maxspecies4]][[maxspecies4]]<-0.0
                  ##update the species1 lineage with a new name essentially
                  genetic_dists[[maxspecies2]]<-genetic_dists[[species1]] ##add the new lineage as a row
                  genetic_dists[[species1]]<-NULL
                  for(i in names(genetic_dists)){ ##add the new lineage as a column
                    genetic_dists[[i]][[maxspecies2]]<-genetic_dists[[i]][[species1]]
                    genetic_dists[[i]][[species1]]<-NULL
                  }
                  ##update the species2 lineage with a new name essentially
                  genetic_dists[[maxspecies3]]<-genetic_dists[[species2]] ##add the new lineage as a row
                  genetic_dists[[species2]]<-NULL
                  for(i in names(genetic_dists)){ ##add the new lineage as a column
                    genetic_dists[[i]][[maxspecies3]]<-genetic_dists[[i]][[species2]]
                    genetic_dists[[i]][[species2]]<-NULL
                  }
                  #update edge length between parent and hybrid node
                  edge.length2[match(maxspecies-1,edge[,2])] <- genetic_dists[[maxspecies4]][[c(maxspecies2,maxspecies3)[parent]]]
                  #update length of reticulate edge
                  hyb_edge.length <- c(hyb_edge.length,genetic_dists[[maxspecies4]][[c(maxspecies2,maxspecies3)[-parent]]])
                  ##update genetic_dists2 in the same way
				  hybrid_species_row<-list()
                  for(i in names(genetic_dists2)){
                    hybrid_species_row[[i]]<- ( (1-inher)*genetic_dists2[[species1]][[i]] ) + ( inher*genetic_dists2[[species2]][[i]] ) ##add row values for the recipient
                    genetic_dists2[[i]][[maxspecies4]]<-hybrid_species_row[[i]] ##add column values for the recipient
                    genetic_dists2[[i]][[species3]]<-hybrid_species_row[[i]]
                  }
                  genetic_dists2[[maxspecies4]]<-hybrid_species_row
                  genetic_dists2[[maxspecies4]][[maxspecies4]]<-0.0
                  genetic_dists2[[maxspecies4]][[species3]]<-0.0
                  genetic_dists2[[species3]]<-hybrid_species_row
                  genetic_dists2[[species3]][[maxspecies4]]<-0.0
                  genetic_dists2[[species3]][[species3]]<-0.0
                  ##update the species1 lineage with a new name essentially
                  genetic_dists2[[maxspecies2]]<-genetic_dists2[[species1]] ##add the new lineage as a row
                  for(i in names(genetic_dists2)){ ##add the new lineage as a column
                    genetic_dists2[[i]][[maxspecies2]]<-genetic_dists2[[i]][[species1]]
                  }
                  ##update the species2 lineage with a new name essentially
                  genetic_dists2[[maxspecies3]]<-genetic_dists2[[species2]] ##add the new lineage as a row
                  for(i in names(genetic_dists2)){ ##add the new lineage as a column
                    genetic_dists2[[i]][[maxspecies3]]<-genetic_dists2[[i]][[species2]]
                  }
                #}
                if(!is_null_trait){
                  trait_state<-c( trait_state,trait_state[del],hyb_trait)
                  trait_state<-trait_state[-del]
                }

                timecreation <- c(timecreation,time,time,time,time)
                maxspecies <- maxspecies-4
                nleaves<-nleaves+1

              }else if( randevent<=(sum(hybprops[1:2])/sum(hybprops)) ){ #Lineage Degenerative
                ##update leaves and extinct species
                parent <- ifelse((1-inher)>0.5,1,2) ##which parent has higher inheritance
                #extinct <- c(extinct,species[2])
                extinct <- c(extinct,species[-parent]) ## the parent with lower inheritance go extinct
                leaves<-c(leaves,maxspecies-1)
                leaves <- leaves[- del]
                ##update edge matrices
                edge <- rbind(edge,c(species[parent],maxspecies-1))
                hyb_edge[[num_hybs]]<-c(species[-parent],species[parent])
                #edge <- rbind(edge,c(species[1],maxspecies-1))
                #hyb_edge[[num_hybs]]<-c(species[2],species[1])
                ##update edge lengths
                edge.length[edgespecevent] <- time-timecreation[- species]
                edge.length<- c(edge.length,0) ##new edge for hybrid lineage
                edge.length2[edgespecevent] <- time-timecreation[- species]
                edge.length2<- c(edge.length2,0) ##new edge for hybrid lineage
                
                #update tree list
                if (Ngene==0) {
                trees2 <- trees
                for (i in 1:length(trees)) {
                    trees[[i]]$edge <- rbind(trees[[i]]$edge,c(species[parent],maxspecies-1))
                }
                trees.prob <- trees.prob * inher #proportion of genes from species[parent]
                for (i in 1:length(trees2)) {
                    trees2[[i]]$edge <- rbind(trees2[[i]]$edge,c(species[-parent],maxspecies-1))
                }
                trees <- c(trees,trees2)
                trees.prob2 <- trees.prob * (1-inher) #proportion of genes from species[-parent]
                trees.prob <- c(trees.prob,trees.prob2)
                }
                if (Ngene>0) {
                    it <- sample.int(2,Ngene,replace=T,prob=c(inher,1-inher))
                    for (i in 1:length(trees)) {
                        if (it[i]==1) {
                            trees[[i]]$edge <- rbind(trees[[i]]$edge,c(species[parent],maxspecies-1))
                        } else {
                            trees[[i]]$edge <- rbind(trees[[i]]$edge,c(species[-parent],maxspecies-1))
                        }
                    }
                }
                
                #if(!is.null(hyb.rate.fxn)){
                  ##update genetic distances
                  species1<-as.character(species[1])
                  species2<-as.character(species[2])
                  maxspecies1<-as.character(maxspecies-1)
                  node <- c(node, species[parent]) #add the hybridization lienages to create two internal nodes
                  tip <- tip[-which(tip==species1)] #remove the hybridization lineages from tip
                  tip <- tip[-which(tip==species2)]
                  tip <- c(tip,maxspecies1) #add the new hybridization lienage to tip
                  if (parent==1) {
                  ##we need to update the recipient lineage
                  hybrid_species_row<-list()
                  genetic_dists[[species1]]<-NULL ##remove species 1 row
                  for(i in names(genetic_dists)){
                    hybrid_species_row[[i]]<- ( (1-inher)*genetic_dists[[i]][[species1]] ) + ( inher*genetic_dists[[species2]][[i]] ) ##add row values for the recipient
                    genetic_dists[[i]][[maxspecies1]]<-hybrid_species_row[[i]] ##add column values for the recipient
                    genetic_dists[[i]][[species1]]<-NULL ##remove species 1 column
                  }
                  genetic_dists[[maxspecies1]]<-hybrid_species_row
                  genetic_dists[[maxspecies1]][[maxspecies1]]<-0.0
                  ##Now we need to get rid of the donor lineage
                  genetic_dists[[species2]]<-NULL
                  for(i in names(genetic_dists)){ ##add the new lineage as a column
                    genetic_dists[[i]][[species2]]<-NULL
                  }
                  edge.length2[match(maxspecies1,edge[,2])] <- hybrid_species_row[[species1]]
                  hyb_edge.length <- c(hyb_edge.length,hybrid_species_row[[species2]])
                  }
                  if (parent==2) {
                  ##we need to update the recipient lineage
                  hybrid_species_row<-list()
                  genetic_dists[[species2]]<-NULL ##remove species 1 row
                  for(i in names(genetic_dists)){
                    hybrid_species_row[[i]]<- ( (1-inher)*genetic_dists[[i]][[species1]] ) + ( inher*genetic_dists[[species2]][[i]] ) ##add row values for the recipient
                    genetic_dists[[i]][[maxspecies1]]<-hybrid_species_row[[i]] ##add column values for the recipient
                    genetic_dists[[i]][[species2]]<-NULL ##remove species 1 column
                  }
                  genetic_dists[[maxspecies1]]<-hybrid_species_row
                  genetic_dists[[maxspecies1]][[maxspecies1]]<-0.0
                  ##Now we need to get rid of the donor lineage
                  genetic_dists[[species1]]<-NULL
                  for(i in names(genetic_dists)){ ##add the new lineage as a column
                    genetic_dists[[i]][[species1]]<-NULL
                  }
                  edge.length2[match(maxspecies1,edge[,2])] <- hybrid_species_row[[species2]]
                  hyb_edge.length <- c(hyb_edge.length, hybrid_species_row[[species1]])
                  }
                  #update genetic_dist2 in the same way
                  hybrid_species_row<-list()
                  for(i in names(genetic_dists2)){
                    hybrid_species_row[[i]]<- ( (1-inher)*genetic_dists2[[i]][[species1]] ) + ( inher*genetic_dists2[[species2]][[i]] ) ##add row values for the recipient
                    genetic_dists2[[i]][[maxspecies1]]<-hybrid_species_row[[i]] ##add column values for the recipient
                    genetic_dists2[[i]][[c(species1,species2)[parent]]]<-hybrid_species_row[[i]]
                  }
                  genetic_dists2[[maxspecies1]]<-hybrid_species_row
                  genetic_dists2[[maxspecies1]][[maxspecies1]]<-0.0
                  genetic_dists2[[maxspecies1]][[c(species1,species2)[parent]]]<-0.0
                  genetic_dists2[[c(species1,species2)[parent]]]<-genetic_dists2[[maxspecies1]]
                #}
                if(!is_null_trait){
                  trait_state<-c( trait_state, hyb_trait)
                  ext_trait_state<-c(ext_trait_state,trait_state[del])
                  trait_state<-trait_state[-del]
                }

                timecreation <- c(timecreation,time)
                maxspecies <- maxspecies-1
                nleaves<-nleaves-1

              }else{ ##Lineage Neutral
                parent <- ifelse((1-inher)<0.5,1,2) ##which parent has lower inheritance will be overflown by the other lineage
                ##update leaves
                leaves<-c(leaves,maxspecies-1:2) ##maxspecies-1 will be the hybrid and 'recipient' lineage while maxspecies-2 will be the 'donor' lineage
                leaves <- leaves[- del]
                ##update edge matrices
                edge <- rbind(edge,c(species[1],maxspecies-1))
                edge <- rbind(edge,c(species[2],maxspecies-2))
                hyb_edge[[num_hybs]]<-c(species[parent],species[-parent])
                ##update edge lengths
                edge.length[edgespecevent] <- time-timecreation[- species]
                edge.length<- c(edge.length,0,0)
                edge.length2[edgespecevent] <- time-timecreation[- species]
                edge.length2<- c(edge.length2,0,0)
                
                #update tree list
                #trees2 <- trees
                #for (i in 1:length(trees)) {
                #    if (parent==1) { #maxspecies-1 is the hybrid lineage
                #        trees2[[i]]$edge <- rbind(trees2[[i]]$edge,c(species[1],maxspecies-1))
                #        trees2[[i]]$edge <- rbind(trees2[[i]]$edge,c(species[2],maxspecies-2))
                #    }
                #    if (parent==2) { #maxspecies-2 is the hybrid lineage
                #        trees2[[i]]$edge <- rbind(trees2[[i]]$edge,c(species[1],maxspecies-1))
                #        trees2[[i]]$edge <- rbind(trees2[[i]]$edge,c(species[1],maxspecies-2))
                #    }
                #}
                #trees.prob <- trees.prob * (1-inher) #proportion of genes from species[1]
                #for (i in 1:length(trees2)) {
                #    if (parent==1) { #maxspecies-1 is the hybrid lineage
                #        trees2[[i]]$edge <- rbind(trees2[[i]]$edge,c(species[2],maxspecies-1))
                #        trees2[[i]]$edge <- rbind(trees2[[i]]$edge,c(species[2],maxspecies-2))
                #    }
                #    if (parent==2) { #maxspecies-2 is the hybrid lineage
                #        trees2[[i]]$edge <- rbind(trees2[[i]]$edge,c(species[1],maxspecies-1))
                #        trees2[[i]]$edge <- rbind(trees2[[i]]$edge,c(species[2],maxspecies-2)) #combine this scenario with the trees_parent==1
                #    }
                #}
                #trees <- c(trees,trees2)
                #trees.prob2 <- trees.prob * inher  #proportion of genes from species[2]
                #trees.prob <- c(trees.prob,trees.prob2)
                for (i in 1:length(trees)) {
                    if (parent==1) { #maxspecies-1 is the hybrid lineage
                        trees2[[i]]$edge <- rbind(trees2[[i]]$edge,c(species[1],maxspecies-1))
                        trees2[[i]]$edge <- rbind(trees2[[i]]$edge,c(species[2],maxspecies-2))
                    }
                    if (parent==2) { #maxspecies-2 is the hybrid lineage
                        trees2[[i]]$edge <- rbind(trees2[[i]]$edge,c(species[1],maxspecies-1))
                        trees2[[i]]$edge <- rbind(trees2[[i]]$edge,c(species[1],maxspecies-2))
                        trees.prob[i] <- trees.prob[i] * (1-inher) #proportion of genes from species[1]
                    }
                }
                if (parent==1) {
                trees2 <- trees
                trees.prob2 <- trees.prob * inher #proportion of genes from species[2]
                    for (i in 1:length(trees2)) {
                        trees2[[i]]$edge <- rbind(trees2[[i]]$edge,c(species[2],maxspecies-1))
                        trees2[[i]]$edge <- rbind(trees2[[i]]$edge,c(species[2],maxspecies-2))
                    }
                trees <- c(trees,trees2)
                trees.prob <- c(trees.prob,trees.prob2)
                }
                
                if (Ngene>0) {
                    it <- sample.int(2,Ngene,replace=T,prob=c(inher,1-inher))
                    for (i in 1:length(trees)) {
                        if (it[i]==1) {
                            if (parent==1) {
                                trees[[i]]$edge <- rbind(trees[[i]]$edge,c(species[2],maxspecies-1))
                                trees[[i]]$edge <- rbind(trees[[i]]$edge,c(species[2],maxspecies-2))
                            }
                            if (parent==2) {
                                trees2[[i]]$edge <- rbind(trees[[i]]$edge,c(species[1],maxspecies-1))
                                trees2[[i]]$edge <- rbind(trees[[i]]$edge,c(species[1],maxspecies-2))
                            }
                        } else {
                            if (parent==1) {
                                    trees[[i]]$edge <- rbind(trees[[i]]$edge,c(species[1],maxspecies-1))
                                    trees[[i]]$edge <- rbind(trees[[i]]$edge,c(species[2],maxspecies-2))
                            }
                            if (parent==2) {
                                    trees[[i]]$edge <- rbind(trees[[i]]$edge,c(species[1],maxspecies-1))
                                    trees[[i]]$edge <- rbind(trees[[i]]$edge,c(species[1],maxspecies-2))
                            }
                        }
                    }
                }
                            
                
                #if(!is.null(hyb.rate.fxn)){
                  ##update genetic distances
                  species1<-as.character(species[1])
                  species2<-as.character(species[2])
                  maxspecies1<-as.character(maxspecies-1)
                  maxspecies2<-as.character(maxspecies-2)
                  node <- c(node,species1,species2) #add the two hybridisation lineages to internal nodes
                  tip <- tip[-which(tip==species1)] #remove them from tip
                  tip <- tip[-which(tip==species2)]
                  tip <- c(tip,maxspecies1,maxspecies2) #add the two resulting new lienages to tip
                  ##First we need to update the donor lineage with a new name essentially
                  if (parent==1) {
                  genetic_dists[[maxspecies2]]<-genetic_dists[[species2]] ##add the new lineage as a row
                  genetic_dists[[species2]]<-NULL
                  for(i in names(genetic_dists)){ ##add the new lineage as a column
                    genetic_dists[[i]][[maxspecies2]]<-genetic_dists[[i]][[species2]]
                    genetic_dists[[i]][[species2]]<-NULL
                  }
                  ##Now we need to update the recipient lineage
                  hybrid_species_row<-list()
                  genetic_dists[[species1]]<-NULL ##remove species 1 row
                  for(i in names(genetic_dists)){
                    hybrid_species_row[[i]]<- ( (1-inher)*genetic_dists[[i]][[species1]] ) + ( inher*genetic_dists[[maxspecies2]][[i]] ) ##add row values for the recipient
                    genetic_dists[[i]][[maxspecies1]]<-hybrid_species_row[[i]] ##add column values for the recipient
                    genetic_dists[[i]][[species1]]<-NULL ##remove species 1 column
                  }
                  genetic_dists[[maxspecies1]]<-hybrid_species_row
                  genetic_dists[[maxspecies1]][[maxspecies1]]<-0.0
                  }
                  if (parent==2) {
                  genetic_dists[[maxspecies1]]<-genetic_dists[[species1]] ##add the new lineage as a row
                  genetic_dists[[species1]]<-NULL
                  for(i in names(genetic_dists)){ ##add the new lineage as a column
                    genetic_dists[[i]][[maxspecies1]]<-genetic_dists[[i]][[species1]]
                    genetic_dists[[i]][[species1]]<-NULL
                  }
                  ##Now we need to update the recipient lineage
                  hybrid_species_row<-list()
                  genetic_dists[[species2]]<-NULL ##remove species 1 row
                  for(i in names(genetic_dists)){
                    hybrid_species_row[[i]]<- ( (1-inher)*genetic_dists[[i]][[species2]] ) + ( inher*genetic_dists[[maxspecies1]][[i]] ) ##add row values for the recipient
                    genetic_dists[[i]][[maxspecies2]]<-hybrid_species_row[[i]] ##add column values for the recipient
                    genetic_dists[[i]][[species2]]<-NULL ##remove species 1 column
                  }
                  genetic_dists[[maxspecies2]]<-hybrid_species_row
                  genetic_dists[[maxspecies2]][[maxspecies2]]<-0.0
                  }
                  #update genetic_dists2 in the same way
                  if (parent==1) {
                  genetic_dists2[[maxspecies2]]<-genetic_dists2[[species2]] ##add the new lineage as a row
                  for(i in names(genetic_dists2)){ ##add the new lineage as a column
                    genetic_dists2[[i]][[maxspecies2]]<-genetic_dists2[[i]][[species2]]
                  }
                  ##Now we need to update the recipient lineage
                  hybrid_species_row<-list()
                  for(i in names(genetic_dists2)){
                    hybrid_species_row[[i]]<- ( (1-inher)*genetic_dists2[[i]][[species1]] ) + ( inher*genetic_dists2[[maxspecies2]][[i]] ) ##add row values for the recipient
                    genetic_dists2[[i]][[maxspecies1]]<-hybrid_species_row[[i]] ##add column values for the recipient
                    genetic_dists2[[i]][[species1]]<-hybrid_species_row[[i]]
                  }
                  genetic_dists2[[maxspecies1]]<-hybrid_species_row
                  genetic_dists2[[maxspecies1]][[maxspecies1]]<-0.0
                  genetic_dists2[[maxspecies1]][[species1]]<-0.0
                  genetic_dists2[[species1]]<-genetic_dists2[[maxspecies1]]
                  edge.length2[match(maxspecies1,edge[,2])] <- hybrid_species_row[[species1]]
                  hyb_edge.length <- c(hyb_edge.length,hybrid_species_row[[species2]])
                  }
                  if (parent==2) {
                      genetic_dists2[[maxspecies1]]<-genetic_dists2[[species1]] ##add the new lineage as a row
                      for(i in names(genetic_dists2)){ ##add the new lineage as a column
                        genetic_dists2[[i]][[maxspecies1]]<-genetic_dists2[[i]][[species1]]
                      }
                      ##Now we need to update the recipient lineage
                      hybrid_species_row<-list()
                      for(i in names(genetic_dists2)){
                        hybrid_species_row[[i]]<- ( (1-inher)*genetic_dists2[[i]][[species2]] ) + ( inher*genetic_dists2[[maxspecies1]][[i]] ) ##add row values for the recipient
                        genetic_dists2[[i]][[maxspecies2]]<-hybrid_species_row[[i]] ##add column values for the recipient
                        genetic_dists2[[i]][[species2]]<-hybrid_species_row[[i]]
                      }
                      genetic_dists2[[maxspecies2]]<-hybrid_species_row
                      genetic_dists2[[maxspecies2]][[maxspecies2]]<-0.0
                      genetic_dists2[[maxspecies2]][[species2]]<-0.0
                      genetic_dists2[[species2]]<-genetic_dists2[[maxspecies2]]
                      edge.length2[match(maxspecies2,edge[,2])] <- hybrid_species_row[[species2]]
                      hyb_edge.length <- c(hyb_edge.length,hybrid_species_row[[species1]])
                  }
                #}
                if(!is_null_trait){
                  trait_state<-c( trait_state, hyb_trait, trait_state[del[2]] )
                  trait_state<-trait_state[-del]
                }

                timecreation <- c(timecreation,time,time)
                maxspecies <- maxspecies-2
              }
              timecreation[-species]<-time
            }
          }
        } else {
          stop = 1
        }
      }
    } ##end while loop

    if (extincttree==0 || age ==0){

      if (extincttree==0){	#length pendant edge
        if(m!=0){
          time<- time+timestep
          age<-time
        }else{
          timestep<-age-time ##we need the length of time between age and the previous event for updating genetic distances
          time = age
        }

        #assign pendant edge length
        for (j in (1:length(leaves))){
          k = which( edge == leaves[j]  ) - length(edge.length)
          edge.length[ k ] <- time - timecreation[- leaves[j]]
          edge.length2[ k ] <- time - timecreation[- leaves[j]]
        }
        timecreation[-leaves]<-age
      }

      if (!is.null(node)) {
      for(i in node){
          for (j in tip) {
            genetic_dists2[[i]][[j]]<- genetic_dists2[[i]][[j]]+timestep
            genetic_dists2[[j]][[i]]<- genetic_dists2[[j]][[i]]+timestep
          }
      }
      }
      for(i in tip){
          for (j in tip) {
              if (i!=j) {
            genetic_dists2[[i]][[j]]<- genetic_dists2[[i]][[j]]+2*timestep
            }
          }
      }
      
      ##used for debugging and testing of the distance matrix creation
      # if(!is.null(hyb.rate.fxn)){
      #   ##update genetic distances
      #   genetic_dists<-lapply(genetic_dists, function(x) lapply(x, function(x) x+(2*timestep))) ##add twice the timestep to all indices
      #   for(i in names(genetic_dists)){ ##lineages shouldn't have a distance from themselves
      #     genetic_dists[[i]][[i]]<- genetic_dists[[i]][[i]]-(2*timestep)
      #   }
      #   print(unlist(genetic_dists))
      #   genetic_dists<-matrix(unlist(genetic_dists),nrow=length(genetic_dists),ncol = length(genetic_dists),byrow=TRUE,dimnames=list(names(genetic_dists),names(genetic_dists)))
      #   colnames(genetic_dists)<- (1:length(colnames(genetic_dists))) [order(as.integer(colnames(genetic_dists)),decreasing=T)]
      #   rownames(genetic_dists)<- colnames(genetic_dists)
      # }

      ##convert hyb_edge to a data frame
      if(num_hybs==0){
        hyb_edge<-matrix(nrow=0,ncol=2)
      }else{
        hyb_edge<-matrix(unlist(hyb_edge),nrow = num_hybs,ncol = 2,byrow = T)
      }
      colnames(hyb_edge)<-c('from','to')
      num_nodes <- abs(maxspecies)
      ## we want to create things into a 'phylo' to report what happened
      ##Convert 'maxspecies' node numberings to ape friendly numberings
      leaf=1 ##Start numbering for leaves
      interior=nleaves+length(extinct)+1 ##Start numberings. start at n+1 because 1:n is reserved for tips
      extinct_tips<-rep(NA,length(extinct)) ##record all extinct tips
      timecreation_order<-rep(NA,length(timecreation)) ##reorder timecreation so node numbers and creation match
      Nextinct<-1
      
      if(!is_null_trait){
        tip_states<-c()
      }
      
      name_table <- matrix(NA,nrow=1,ncol=2)
      for (j in (1:num_nodes)){
        inleaves <- (leaves  %in% (-j))
        inextinct<- (extinct %in% (-j))
        if (!(any(inleaves) | any(inextinct))){ ##we are at an internal node
          replaced_value <- interior
          interior <- interior +1
        } else { ## We are at a tip of some sort (extant,extinct, or hybrid)
          
          if(!is_null_trait){##if a trait model
            
            #record the trait value of the leaf we are at
            
            ##determine if we are at an extant or extinct tip
            if(any(inleaves)){ #we are in an extant leaf
              trait_value <- trait_state[which(-j==leaves)]
            }else{ ##we are extinct leaf
              trait_value <- ext_trait_state[which(-j==extinct)]
            }
            tip_states<-c(tip_states,trait_value)
          }
          
          replaced_value <- leaf
          leaf <- leaf +1
          if(sum(match(extinct,-j,0))==1){
            extinct_tips[Nextinct]<-replaced_value
            Nextinct<-Nextinct+1
          }
        }
        timecreation_order[replaced_value] <- timecreation[j]
        hyb_edge[hyb_edge== - j]<-replaced_value
        edge[edge == -j]<-replaced_value
        name_table <- rbind(name_table,c(-j,replaced_value))
      }

      ##create a list with all the 'phylo' bits
      phy <- list(edge = edge)
      phy$tip.label <- paste("t", sample(nleaves+length(extinct)), sep = "")
      phy$edge.length <- edge.length2
      phy$edge.length.time <- edge.length
      phy$Nnode <- num_nodes-(nleaves+length(extinct)) ##Nnode reports the number of internal nodes
      class(phy) <- c("evonet","phylo")
      phy$reticulation<-hyb_edge
      phy$reticulation.length <- hyb_edge.length
      #phy$dists<-genetic_dists
      phy$inheritance<-inheritance
      phy$nleaves<-nleaves

      hyb_tip_indices<-match(extinct_tips,phy$reticulation[,1],0)>0
      phy$hyb_tips<-extinct_tips[hyb_tip_indices] ## Keep track of tips that have an edge going into a hybrid one. These will be fused later
      phy$extinct<-phy$tip.label[extinct_tips[!hyb_tip_indices]] ##Keep track of extinct tips. Will be used if using the reconstructed tree

      if(!is_null_trait){
        phy$tip.states<-tip_states
      }

      if(n!=0){ ##we use these for GSA
        phy$timecreation <- timecreation_order
        if(is.null(time_in_n)){
          time_in_n<-NA ##We were never in the correct n
        } else{
          time_in_n<- time_in_n[!(time_in_n %in% unique(time_in_n[duplicated(time_in_n)]))]
          time_in_n<-matrix(time_in_n,ncol=2,byrow=T)
        }
        phy$time_in_n<-time_in_n
      }
    }
    if (!is.null(phy$edge)) {
    phy$node.label <- name_table[match(phy$edge[1,1]-1+c(1:phy$Nnode),name_table[,2]),1]
    phy$tip.label <- name_table[match(1:(phy$edge[1,1]-1),name_table[,2]),1]
    }
    #translate trees into newick format
    trees.new <- trees
    for (i in 1:length(trees.new)) {
        for (j in 1:dim(name_table)[1]) {
            trees.new[[i]]$edge[trees.new[[i]]$edge==name_table[j,1]] <- name_table[j,2]
        }
        trees.new[[i]]$Nnode <- phy$Nnode
        trees.new[[i]]$tip.label <- phy$tip.label
        class(trees.new[[i]]) <- c("evonet","phylo")
        trees.new[[i]] <- ape::write.tree(trees.new[[i]])
    }
    
    list(phy=phy,name_table=name_table,distance=genetic_dists2,trees=trees,trees.prob=trees.prob,trees.new=trees.new)
  }

#example
#library(ape)
#library(SiPhyNetwork)
#inheritance.fxn <- make.uniform.draw()
#sim <- sim.bdh.age(age=2,numbsim=1,lambda=1,mu=0.2,nu=0.25,hybprops=c(1,0,0),hyb.inher.fxn=inheritance.fxn,complete=F,Ngene=100)
#phy <- sim[[1]]$phy
#plot(phy,show.node.label=T)
#edgelabels(phy$edge.length,frame = "none", adj = c(0.5, -0.5))

#get unique trees and count how many genes have the same tree
#unique_trees <- unique(trees.new)
#n_trees <- length(unique_trees)
#count <- numeric(n_trees)
#for (i in 1:length(trees.new)) {
#    for (j in 1:n_trees)
#        count[j] <- count[j] + all(trees.new[[i]]==unique_trees[[j]])
#}

