# simplet_counting

These are source codes for the paper "Estimating Simplet Counts via Sampling".

# **SC3: Simplet Counting via Color Coding**
* **SC3** is the algorithm for counting simplets of size 4 and 5 based on the Color Coding method where simplets are isomorphic classes of connected induced subcomplexes.
* **SC3**, which extends CC (Color Coding) method from graphs to simplicial complexes, has the following properties:
  1. accurate compared to baselines
  2. fast and scalable
* Especially, the counts of simplets by SC3 have strong power to characterize SCs domain by domain.

# **SCRW: Simplet Counting via Random Walk**
* **SCRW** is the restricted access model on simplicial complexes for estimating simplet concentration.
* **SCRW**, which extends RW (Random Walk) method from graphs to simplicial complexes, has the following properties:
  1. comparable speed and accuracy compared to SC3
  2. memory efficient

## Datasets 
1. Real-world datasets are available [here](https://www.cs.cornell.edu/~arb/data/). The format of the data is exactly as provided, and for SC3, only "{data}-nverts.txt" and "{data}-simplices.txt" are necessary.
2. Each dataset we used in our SCRW experiments is the LCC of the original one. This was pre-processed to satisfy maximality and is available [here](https://postechackr-my.sharepoint.com/personal/hhyy0401_postech_ac_kr/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fhhyy0401%5Fpostech%5Fac%5Fkr%2FDocuments%2F%EC%97%B0%EA%B5%AC%2D%EC%98%A8%EB%9D%BC%EC%9D%B8%2Fmaximal%5FLCC&ga=1).
3. [Statistics of datasets we used]

| dataset  | # of vertices | # of maximal simplices |
|----|---:|---:|
| coauth-DBLP | 1,924,991 | 1,730,664 |
| coauth-MAG-Geology | 1,256,385 | 925,027 |
| coauth-MAG-History | 1,014,734 | 774,495 |
| congress-bills  | 1,718 | 48,898 |
| contact-high-school | 327 | 4,862 |
| contact-primary-school | 242 | 8,010 |
| DAWN  | 2,558 | 72,421 |
| email-Eu | 998 | 8,102 |
| email-Enron | 143 |433 |
| NDC-classes | 1,161 | 563 |
| NDC-substances | 5,311 | 6,555 |
| tags-ask-ubuntu | 3,029 | 95,639 |
| tags-stak-overflow  | 49,998 |3,781,574 |
| threads-ask-ubuntu |125,602 | 149,025 |
| threads-math-sx  |  176,445 | 519,573  |
| threads-stack-overflow  |  2,675,955 |  8,694,667 |
  
## Outputs
1. (csv) timestamp for each step (directory: *result\timestamp*)
2. (csv) count of every simplet of size k with a format ({id of a simplet}, {counts}\n) (directory: *result\CC* for SC3, *result\RW* for SCRW)


## Running Codes
1. Select algorithms between SC3 (directory: *src_CC*) and SCRW (directory: *src*). 
2. command "bash run.sh". Feel free to change the variable **k, ss, datas, trials, threads** in the file "run.sh"
```
k="4 5 6" # a list of the size of simplets
ss="100 1000 10000 100000" #a list of the number of samples
datas="toy" #a list of the names of datasets
trials=5 #the number of trials for running SC3/SCRW
threads=6 #the number of multi-threads
```
3. For SC3, choose a version between SC3 and SC3-E (variable: **ver** in the file "run.sh").
```
ver=1 #origial SC3
ver=2 #memory efficient SC3: SC3-E
```
