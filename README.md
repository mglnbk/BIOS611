# BIOS611 Project: Pocker Hand Pattern Recognition


## Problem Settings

Problem Statement: "While Poker hands are governed by deterministic rules, the goal of this project is to evaluate the capacity of ML method (XGboost here) to approximate complex non-linear functions based on our re-sampled training data.

Specifically, this project aims to do:

- Discover combinatory logic (e.g., pairs, sets) in unsupervised manner.

- Overcome extreme class imbalance through carefull feature engineering

- See how big the influence of prior knowledge on model performance.


## Dataset description
Each record is an example of a hand consisting of five playing cards drawn from a standard deck of 52. Each card is described using two attributes (suit and rank), for a total of 10 predictive attributes. There is one Class attribute that describes the "Poker Hand". The order of cards is important, which is why there are 480 possible Royal Flush hands as compared to 4 (one for each suit - explained in ftp://ftp.ics.uci.edu/pub/machine-learning-databases/poker/poker-hand.names). The class is purely based on the rules of Texas Holdem.

## Model setup

I trained a XGboost model combined with manual feature engineering since the model performs not very well without feature engineering.



## Instructions for replicating the project

1. git clone this repo, `git clone ...`
2. checkout proj branch: `git checkout proj`
3. In the root directory of git repo, do `make build` to build a Docker env called my name
4. This step is optional, you can call `make clean` before step 5 to ensure you can re-run all the analysis
5. Do `make run` to generate all the figures and report

Note: 

1. Some steps may take long time
2. I use `FROM rocker/tidyverse:latest` this repo solves **TONS** of dependency issues rather than the one provided in project announcement


## Result

See report/03\_report.pdf for a conclusional summary.

See saved_model/lastest_model.model for trained XGboost model.

See plot/*.png for some exploratory visualization.
