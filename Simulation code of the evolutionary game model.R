library(sna)
library(network)
DataF=NULL
DataD=NULL

# Setting parameters
N = 100 # network size
d = 1 # original network density
round_numbers = 200 
num_replications = 100
b=2
Beta = seq(55,95,10)/1000
Delta = seq(0,40,2)

for (beta in Beta){
  for (delta in Delta){
    data=NULL
    for (cases in 1:num_replications){
      
      ## Initial conditions: round zero
      # Action sets
      actions = rep(0,N)
      actions[sample(seq(1,N),N/2, replace = FALSE)] =1 # 50% are cooperators (denoted by 1)
      reputations = actions
      
      # Individual profiles
      accounts = NULL
      strategy = NULL
      
      # Networks: fully linked in the original state
      original_Net=matrix(rep(0,N^2),ncol=N,nrow=N)
      upper_tri_index = which(upper.tri(original_Net, diag = FALSE) == TRUE)
      original_Net[sample(upper_tri_index, 
                          d*length(upper_tri_index), 
                          replace = FALSE)] = 1
      original_Net=original_Net+t(original_Net)
      game_network = original_Net 
      
      # Game dynamics
      for (r in 1:round_numbers){
        
        # Play game and determine payoff
        payoff = rep(0,N)
        for (i in 1:N){
          neighbors = which(game_network[i,]==1)
          if (length(neighbors)>0){ # if linked to any neighbor
            for (j in neighbors){
              payoff[j] = payoff[j]+ actions[i]*b # recall action=1 means cooperation
              payoff[i] = payoff[i]+ actions[j]*b + abs(actions[i]-1)
            }
          }
        }
        accounts = rbind(accounts,payoff) # keep track of payoffs of each round
        strategy = rbind(strategy, actions) # keep track of behavioral strategies of each round
        
        # Strategy update
        for (i in 1:N){
          # Randomly select a neighbor for social learning of behavioral strategy
          neighbors = which(game_network[i,]==1)
          if (length(neighbors)>0){
            if (length(neighbors)>1){
              target = sample(neighbors,1,replace = FALSE)
            }
            else{
              target = neighbors
            }
            # Decide whether to adopt the strategy of the (randomly) selected neighbor
            tao = 1/(1+exp(beta*(payoff[i]-payoff[target])))
            if (runif(1) < tao){
              actions[i] = actions[target]
            }
          }
        }
        reputations = reputations + actions # update reputation records
        
        # Network updating: formation
        formed_network=matrix(rep(0,N^2),ncol=N,nrow=N)
        if (r==2){# Beginning from round 2, (game) network starts to evolve from a null network  
          game_network = formed_network
        }
        
        for (i in 1:N){
          non_neighbors = which(game_network[i,]==0)# identify players who are not linked to the focal player
          non_neighbor_cooperativeness = reputations[non_neighbors]/(r+1) # assessing their cooperativeness
          formation_prob = 1/(1+exp(-1*delta*(non_neighbor_cooperativeness-0.5))) # calculate the probability of proposing a tie
          unif_set = runif(length(non_neighbors)) # grab a set of random numbers between 0 and 1
          to_be_formed = non_neighbors[which(unif_set < formation_prob)]
          formed_network[i,to_be_formed]=1 # unilateral proposals of links
        }
        game_network = game_network + formed_network*t(formed_network) # bilateral decisions of newly formed links
        diag(game_network) = 0 # self-link not allowed
      }# loop over rounds
      
      
      #-----The following code is for demonstration only!-----
      #-----It shows how we process network evolution under the dissolution rule
      #-----Please do not execute it along with the previous lines of code on 'network formation'
      #-----Other than this part, the 'network formation' and 'network dissolution' games share the same code
    
      for (i in 1:N){
        neighbors = which(game_network[i,]==1) # keep track of extant network neighbors
        neighbor_cooperativeness = reputations[neighbors]/(r+1) # calculate their cooperativeness
        deletion_prob = 1/(1+exp(-1*delta*((1-neighbor_cooperativeness)-0.5))) # calculate the probabilities of deleting these neighbors
        unif_set = runif(length(neighbors)) # grab a set of random numbers between 0 and 1
        to_be_deleted = neighbors[which(unif_set < deletion_prob)]
        game_network[i,to_be_deleted]=0
        game_network[to_be_deleted,i]=0 # network dissolution is a unilateral decision
      }
      #-----end of this portion-----
      
      # Record network attributes in the final round
      # The following functions draw on the packages of 'network' and 'sna' 
      net=network(game_network,matrix.type="adjacency") # processing network into an object
      density_mean=mean(degree(net,gmode="graph"))/(N-1) # calculating average degree
      compnt_max=max(component.dist(net, connected="strong")$csize)/N # calculating the (relative) maximum size of component 
      clustering_mean=mean(gtrans(net)) # calculating network clustering (the probability of transitivity)
      
    }#cases
    
  } # Delta
} # Beta

