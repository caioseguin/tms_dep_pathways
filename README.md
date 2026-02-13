%%% Neuroanatomical pathways of TMS therapy for depression
%%% Seguin et al. 2025, medRxiv (https://www.medrxiv.org/content/10.1101/2025.02.10.25322034v2)

% This script runs the main analyses of the manuscript. 
% It computes analyses and generates the figures of the manuscript.
% The script uses publicly available data from Weigand et al 2018 (Biol Psychiatry), named Cohort I in our study. 
% Data from Cohort II is not publicly available and therefore analyses using this patient group are not reproduced here. 
% We perform analyses on two normative connectomes constructed from publicly available datasets from 
% Griffa et al 2019 (Lausanne connectome) and Mansour et al 2021 (Schaefer connectome). 
% We use functions from the BCT (Rubinov & Sporns 2010) to compute network communication measures.

%%% Caio Seguin, 22 Oct 2025