# eNODEconstr
eNODEconstr: An ecology-guided and knowledge-constrained framework for inferring microbial community dynamics

## Framework of eNODEconstr
![](https://github.com/tanxiaoxiu/eNODEconstr/blob/master/Framework.png)


## a. Learning 
eNODEconstr models the dynamics of microbial and metabolite abundances using the Microbial Consumer Resource Model (MiCRM), which is represented by a set of ordinary differential equations. These ecological equations are incorporated into a deep learning framework known as Neural Ordinary Differential Equations (NeuralODEs). The ecology-guided NeuralODEs are trained on time-series abundance data for microbes and metabolites, with prior knowledge of inter-species cross-feeding interactions encoded as constraints and integrated into the learning process via a regularization penalty.

## b. Inference
Using the trained model, eNODEconstr predicts unknown cross-feeding interactions that are missing from the current knowledge base. Additionally, by solving the initial value problem with an ODE solver, it forecasts the trajectories of microbial and metabolite abundances. 

## c. Manipulation
Through a thought experiment in which a microbe is hypothetically removed from the community, eNODEconstr calculates the composition score and function score to quantify the contribution of each species.

