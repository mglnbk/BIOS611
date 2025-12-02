# BIOS611 Project: Pocker Hand Pattern Recognition


## Problem Settings

Problem Statement: "While Poker Hands are governed by deterministic rules (logic), the goal of this project is to evaluate the capacity of Machine Learning classifiers to approximate complex non-linear functions solely from observed data.

Specifically, this project challenges the model to:

- Discover combinatory logic (e.g., pairs, sets) without explicit rule specification.

- Handle invariant properties (e.g., card order does not matter).

- Overcome extreme class imbalance (Class 9 represents < 0.001% of instances)."


## Dataset description
Each record is an example of a hand consisting of five playing cards drawn from a standard deck of 52. Each card is described using two attributes (suit and rank), for a total of 10 predictive attributes. There is one Class attribute that describes the "Poker Hand". The order of cards is important, which is why there are 480 possible Royal Flush hands as compared to 4 (one for each suit - explained in ftp://ftp.ics.uci.edu/pub/machine-learning-databases/poker/poker-hand.names).

## Model setup

I trained a XGboost model combined with manual feature engineering since the model performs not very well without feature engineering.

Feature engineering process:


Model parameter settings:


## Instructions for replicating the project



## Result


