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
  Year -> POPULATION ;
  Oceanography -> Prey -> Adult -> Pup -> POPULATION ;
  Oceanography -> Prey -> Pup -> POPULATION ;
  eSeal -> Prey -> Pup ;
  eSeal -> Prey -> Adult ;
  eSeal -> Prey -> Molt ;
  Site -> Coyote -> Pup -> POPULATION ;
  Site -> Pup -> POPULATION;
  Site -> Adult -> POPULATION;
  Site -> Coyote -> Adult -> POPULATION ;
  Site -> Coyote -> Molt -> POPULATION ;
  Adult -> Pup -> POPULATION;
  Adult -> Molt -> POPULATION;
  
  
  Year [exposure] ;
  Site [exposure] ;
  Oceanography [exposure] ;
  Coyote [exposure] ;
  Adult [exposure] ;
  Pup [exposure] ;
  Molt [exposure] ;
  POPULATION [outcome] ;
  Prey [unobserved] ;
  eSeal [exposure] ;
  
  

}")

## 2024-06-11
## add sharks and eSeals (and Orca)?-

## add in zone?
par(mfrow = c(1,1))

#plot(DAG_PHOCA)
impliedConditionalIndependencies(DAG_PHOCA)


# to pretty up the plot 
coordinates(DAG_PHOCA) <- list(x=c(Year=2,
                                   Site=2.5,
                                   Coyote=3,
                                   Oceanography=2,
                                   Pup=2,
                                   Adult=3, 
                                   Molt=4, 
                                   POPULATION=3,
                                   Prey=2),
                              y=c(Year=-1,
                                  Site=-4,
                                  Coyote=-4,
                                  Oceanography=-4,
                                  Pup=-2,
                                  Adult=-2, 
                                  Molt=-2, 
                                  POPULATION=-1,
                                  Prey=-3))
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