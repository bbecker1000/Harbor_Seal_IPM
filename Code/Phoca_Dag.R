##DAGs for Phoca MARSS models

library(dagitty)

## Dependent Var
# counts of seals at age class by site

## Independent covariates
# Year
# Site  (only a few early dates so probably skip)
# MEI
# PDO
# ENSO
# Coyote


## Random Effects?

# build the DAG

set.seed(123)

DAG_PHOCA <- dagitty("dag{ 
  MOCI -> Prey -> Adult -> Pup;
  Coyote -> Adult -> Pup ;
  Adult -> Pup ;
  Adult -> Molt ;
  Disturbance -> Adult -> Pup;
  Disturbance -> Molt ;
  
  
  

  MOCI [exposure] ;
  Coyote [exposure] ;
  Adult [outcome] ;
  Pup [outcome] ;
  Molt [outcome] ;
  Prey [unobserved] ;
  Disturbance [exposure] ;
  
  

}")



## 2024-06-11
## add sharks and eSeals (and Orca)?-

## add in zone?
par(mfrow = c(1,1))

#plot(DAG_PHOCA)
impliedConditionalIndependencies(DAG_PHOCA)


# to pretty up the plot 
coordinates(DAG_PHOCA) <- list(x=c(Coyote=2.5,
                                   MOCI=2.25,
                                   Disturbance=2.75,
                                   Pup=2.5,
                                   Adult=2.5, 
                                   Molt=2.75, 
                                   Prey=2.25),
                              y=c(Coyote=-4,
                                  MOCI=-4,
                                  Disturbance=-4,
                                  Pup=-2,
                                  Adult=-3, 
                                  Molt=-2, 
                                  Prey=-3.5))
plot(DAG_PHOCA)





#####----------------------------------
#try with ggdag
library(ggdag)



Phoca_dag <- dagify(
  POPULATION ~ Year,
  Prey ~ Oceanography,
  Pup ~ Prey,
  Adult ~ Prey,
  POPULATION ~ Pup,
  POPULATION ~ Adult,
  POPULATION ~ Molt,
  Pup ~ Coyote,
  Adult ~ Coyote,
  Coyote ~ Site,
  
  exposure = "Coyote",#"Oceanography",
  outcome = "POPULATION"

)
ggdag(Phoca_dag, text_col = "white", 
      stylized = TRUE) + 
  theme_dag_blank()
ggdag_paths(Phoca_dag, 
            #adjust_for = c("Breach"),
            text_col = "black")

ggdag_adjustment_set(Phoca_dag, text_col = "black")

library(dagitty)

adjustmentSets(goby_dag)