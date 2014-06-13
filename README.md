
# Details of _H. impatiens_ analyses #

Initially this document was linear and merged together details about the
analyses for (1) all sequences, (2) just COI for doing the GMYC and (3) the
*BEAST analyses on the ESU1 complex. I'm now reformatting it so it's divided in
three sections (and an additional one for early analyses that I didn't keep). 

# Starting over the analyses

* For PartitionFinder, even if I specify the models that are implemented by
  BEAST as an option, some of the models selected are not directly implemented
  with BEAUTI, so I re-run PartitionFinder saying that I wanted a model that was
  JC, HKY, GTR with or without I and/or Gamma.
* For *BEAST, the trouble of having to deal with the phased dataset doesn't seem
  totally worth it given that there are so few ambiguities and that the extra
  data generated is going to slow down the analyses,  therefore I decided to go
  with the raw sequences.
* Given that the reduced dataset may reduce the ability to estimate some of the
  parameters for the models of molecular evolution, re-submitted partitionfinder
  analysis to UF HPC (submit queue) with only a subset of the sequences.
* Given that I don't want to wait for the results fo the former (still hasn't
  been submitted after 3+ hours), I started the *BEAST analyses using the same
  kind of models as with the full dataset, just specified "empirical" for the
  base frequencies.
* Tested briefly the code for path and stepping stone sampling on my laptop with
  a really short run, and it seems to be working. Submitted 1 XML file with all
  impatiens ESU specified (`20130929.starBeast_all.xml`) to CIPRES with maximum
  allowed time (not sure it's going to be enough); and one file with RedSea ESU
  replaced with ESU1 (`20130929.starBeast_noRedSea.xml`) to CIPRES (under
  sinqueso account). For both I'm using strict clocks with COI clock fixed
  at 0.037 and the other estimated with Uniform prior 0-10 (or 20).
* NOTE: when editing manually the *BEAST file to change species attribution,
  make sure to also change speciesTree.splitPopSize so that there is the same
  number of 1 as of species, and remove 2 '0' for each species removed.

# *BEAST analyses

## 20131002
* on CIPRES (under `fmichonneau@iplant` account) running:
  - noWpacHawaiiRedSea: ESU1, ESU3, gracilis
  - noWpacRedSea: ESU1, ESU3, gracilis, Hawaii
  - noWpacHawaii: ESU1, ESU3, gracilis, RedSea
  - noWpac:       ESU1, ESU3, gracilis, RedSea, Hawaii
  - allESUs:      ESU1, ESU3, gracilis, RedSea, Hawaii, WPac
* on CIPRES (under sinqueso account) running:
  - noRedSea:     ESU1, ESU3, gracilis, Hawaii, WPac
* on UF's HPC (submit queue):
  - noHawaii:     ESU1, ESU3, gracilis, RedSea, WPac

## 20131007 (not kept)
* The first analyses finished today (20131006), definitely problems with IO of
  the log files, portions of the sampled posterior are missing completely. This
  problem was due to too much logging (I said every 5000 instead of 50,000)
  [CIPRES sysadmin says there was a problem on their end]. The calculations for
  the marginal likelihood seem fine. Haven't checked the trees yet. Quick
  overview of the output in tracer seems to indicate that mixing is correct.
* Edited the XML files:
  - change prefix of files to 20131007
  - change the logging to every 50,000
* The noRedSea analyses also didn't run as expected, and I haven't been able to
  retrieve output, therefore, re-running it. Exact same thing as before but
  logging every 50000 instead (didn't keep copy of old file.)
  - noRedSea          --  seed 19820327 -- label: `20131002_impatiens_noRedSea` (on CIPRES)
* Analyses submitted (on CIPRES under fmichonneau@iplant account):
  - noWpacHawaii       -- seed 19830304 -- label: `20131007_impatiens_noWpacHawaii`
  - noWpacHawaiiRedSea -- seed 19830304 -- label: `20131007_impatiens_noWpacHawaiiRedSea`
  - allESU1impatiens   -- seed 27031982 -- label: `20131007_impatiens_allESU1`
  - assigned individuals to random species -- seed 27031982 -- label: `20131007_impatiens_random`
  - noWpacRedSea       -- seed 27031982 -- label: `20131007_impatiens_noWpacRedSea`
  - noWpac             -- seed 19820304 -- label: `20131007_impatiens_noWpac`
  - noRedSea           -- seed 27031982 -- label: `20131007_impatiens_noRedSea`
* Analyses submitted on submit queue (UF's HPC)
  - noHawaii           -- seed 27031982

## 20131009 (in `impatiens_analyses/001.noMt`)
* redoing analyses without using mt data
  - kept the same partition names as for the full dataset for the subst rates,
    based on partitionfinder output
  - using emprical base frequencies
  - Yule process -- piecewise linear & constant root

## 20131226 -- starBEAST
* To test if the differences in marginal likelihoods is due to convergence
  problems, I am re-running the estimation of the likelihoods for all the
  species in ESU1 with the exact same parameters as before, except that, I set
  the number of iterations to 3e6 with 2 independent runs, 1 with the same seed
  as one of the previous runs, and one with a random seed.




- - - - - - -

# COI tree for GMYC analyses

## 20131008 -- COI tree
* the tree produced using only COI sequences for GMYC wasn't too good
  (20131008), redoing it with:
  - remove duplicate sequences (using RAxML, 75 sequences are exactly
    identical), file 20131009.impatiens\_COI\_uniq.phy
  - use HKY+I+G model
  - use relaxed clock (uniform prior [0,10000] on ucld.mean)
  - file 20131009\_impatiens\_COI.xml

## 20131118
* trying again to generate COI tree suitable for GMYC analysis.
* Kept unique haplotypes
* use HKY+I+G model
* use random local clock (following email on BEAST mailing list from A. Drummond
  on March 5th 2012), with clock.rate prior uniform [0,1e100]
* file `20131118_impatiens_COI_uniq.xml`

## 20131120
* the random clock analysis poses problem with BEAGLE and crashes. Running
  20131118 on HiPerGator with single thread to avoid problem
* decided to also re-run strict clock to see how it compares
* analysis on:
  - unique haplotypes
  - HKY+I+G model
  - strict clock, fixed at 0.037
  - file `20131120_impatiens_COIuniq_strict.xml`
* running on CIPRES

## 20131125 COI
* redoing COI tree with unique sequences and no outgroup
* running it for both strict and lognormal clock
* HKY+I+G model
* coalescent: constant size model
* for strict:
  - fixed at 0.037
* for relaxed:
  - lognormal estimated
  - put prior on age of nodes for tmrca(WA+Gala+EP) with lognormal [1.5, 1]

## 20131126 COI
* same as before but changing tree prior to coalescent growth.

## 20140421 -- Partitionfinder and BEAST on COI sequences
* Given that the number of estimated species seem to vary depending on the
  parameters used for the BEAST analysis, I decided to use a more methodological
  approach to understand the relative importance of each.
* First, I re-ran partitionfinder on the COI alignment to figure out best
  partitioning scheme.
* To create alignment, got latest from files `20130923-185350-COI.fas` wrote
  script  to make:
  - PHYLIP file for partitionfinder
  - NEXUS file for BEAST
  

```r
  library(seqManagement)
  setwd("~/Documents/Impatiens/20140421.COI_partitionFinder")
  system("muscle -in 20130923-185350-COI.fas -out 20130923-185350-COI.afa")
  fas2phy(file="20130923-185350-COI.afa", toDrop="S0213", overwrite=TRUE)
  file.rename("20130923-185350-COI.phy", "20140421.COI.phy")
  alg2nex(file="20140421.COI.phy", format="seq", interleaved=FALSE,
            partition.file="20140421.COI_partitions.part")
  ## get extraction numbers for EP+WA+Gala
  alg <- ape::read.dna(file="20140421.COI.phy", format="seq")
  nm <- dimnames(alg)[[1]]
  m <- match(nm, impDB$Extract)
  tmpDB <- impDB[m, ]
  sort(subset(tmpDB, consensusESU %in% c("EP", "WA", "Gala"))$Extract)
```

* According to partitionfinder best model is:
  - 1st codon position: HKY+I+G
  - 2nd codon position: GTR+G
  - 3rd codon position: HKY
* BEAST parameters to vary:
  - include all sequences vs. include only different sequences
  - use strict or relaxed clock model
  - use constant or exponential growth with coalescent prior; or speciation Yule
* BEAST parameters that won't vary
  - calibrate divergence of isthmus of Panama
  - allSeq\_strict/relaxed\_yule/coalcst/coalexp
  - each in their own folder
  - allSeq\_strict\_yule 1e7 generations with sampling every 1e3, all others with
    5e7 generations and sampling every 5e3.

## 20140422 -- BEAST on COI sequences
* It turns out that the the HiPer Gator cluster doesn't have BEAST 2.1.2 yet,
  and the files I generated don't work with BEAST 2.0.2. I'm thus going to use
  BEAST 1.8.0 so I can use Cipres and submitted all my analyses at once.
* Define EP+WA+Gala as group with MRCA set to lognormal distribution mean (1.5),
  sd (0.75) and offset (2.5) 
* Strict clock: estimated with lognormal prior mean (-4), SD (1)
  Relaxed clock: estimated with lognormal prior mean (-4), SD (1)
  
### Results for these analyses
* The trace files look good. However the values of the marginal likelihoods are
  all over the place and indicate that something is going on. I emailed on
  20140430 Guy Baele to see what could be the problem here. 

## 20140428 -- BEAST on COI sequences
* To make code easier
* I now need to do the same analyses as for 20140422 but with removing the
  duplicated sequences. It also seems that using HKY+I+G is causing some
  convergence problems, so changing to HKY+G
* Furthermore, the path-sampling stepping-stone is for some reason very slow and
  completely inaccurate, and not that important for this analysis. So I'm
  dropping it. If it turns out that it's something I want to do in the future,
  it won't take too long. At this stage, I want to know if the other parameters
  are influencing the number of species estimated by the various GMYC
  algorithms.
* To generate the NEXUS file for the analyses:
  - get `20140421.COI.phy` from the partitionfinder folder
  - run RAxML to remove duplicated sequences and obtain file
    `20140428.COI_noDup.phy`
  - create partition file for each codon position
  - converts to NEXUS using `alg2nex()`


```r
  library(seqManagement)  
  setwd("~/Documents/Impatiens/20140428.noDup_datapreparation")
  alg2nex(file="20140428.COI_noDup.phy", format="seq", interleaved=FALSE,
            partition.file="20140428.COI_noDup.part")
```

- then in this folder, generate all the BEAST input files (using
    BEAST 1.8.0)
  - to get the tips in the WA+EP+Gala clade:
  

```r
  library(ape)
  setwd("~/Documents/Impatiens/20140428.noDup_datapreparation")
  impDB <- read.csv(file="~/Documents/Impatiens/impatiens_phylogeography/data/impatiensDB.csv",
                      stringsAsFactors=FALSE)
  alg <- ape::read.dna(file="20140428.COI_noDup.phy", format="seq")
  nm <- dimnames(alg)[[1]]
  m <- match(nm, impDB$Extract)
  tmpDB <- impDB[m, ]
  sort(subset(tmpDB, consensusESU %in% c("EP", "WA", "Gala"))$Extract)
```

- Used the same parameters as in the previous analyses (`allSeq`) EXCEPT for
    WA+EP+Gala constrained as monophyletic in this analysis and not in the
    other!
	
### Results for these analyses
* *noDup_relaxed_coalcst* no convergence, needs to be run longer
* *noDup_relaxed_coalexp* no convergence, needs to be run longer
* *noDup_relaxed_yule* no convergence, needs to be run longer
* *noDup_strict_coalcst* no problem, almost all ESS >> 300
* *noDup_strict_coalexp* no convergence, needs to be run longer
* *noDup_strict_yule* barely converged, all ESS just about 100


## 20140429 -- BEAST on COI sequences

### Rationale

* Even after re-running all the analyses on the dataset without the duplicated
  sequences for another 50e7 generations, I can't get convergence and combining
  the results from the two runs are not making things better.
* It's a little unclear whether this lack of convergence/mixing changes anything
  (the bad mixing mostly affects "prior") but it looks like there is
  overparametrization.
* An easy way to reduce the overparametrization is by making the partitioning
  scheme simpler. So I'm re-doing a partitionfinder on the alignemnt to see if I
  can do a simpler model.

### Results

* The exact same model as before is being selected... so not the best way to
  reduce number of parameters.


## 20140430 -- BEAST on COI sequences

### Rationale

Given the results of 20140422, 20140428 and 20140429, I'm going to re-run the
analyses using a simple model of molecular evolution for COI (just 1 partition
with GTR + G) to see if that makes things easier (better convergence and
possible to estimate marginal likelihoods). I'm going to run it on (1) all
sequences and (2) no duplicated sequences; for relaxed clock and coalcst models.

### Methods

* For _noDup_ use file in 20140428.noDup_datapreparation
* For _allSeq_ use file in 20140421.COI_partitionfinder
* Same priors as before:
  - GTR+G model of molecular evolution
  - lognormal relaxed clock with prior mean -4, SD 1
  - tmrca for EP+WA+Gala set as lognormal with mean 1.5, SD 0.75, offset 2.5
  - 5e7 generations with sampling every 5000
  - estimation of log marginal
  - CIPRES is down, submitted to diagriid (not sure they sent notifications when
    job is done...)

### Update (20140501)

* job on diagriid still running
* submitted jobs on CIPRES nonetheless
* Guy Baele suggested that the underflow computations were coming from the
  beagle scaling, and suggested to use beagle_scaling always; HiPerGator is
  still running the all impatiens analysis, so running on my computer. Took me
  most of the afternoon to figure out how to install beagle on my laptop. Now
  running.

### Update (20140505)

* analysis on diagrid still running (at least one of them), still issues with
  underflow.
* analyses still running on CIPRES, killed and deleted (also showing underflow)
* resubmitted the analysis on CIPRES without using BEAGLE.

### Update (20140508)

* Analyses submitted on May 5th haven't completed yet, but the log doesn't show
  any issues with underflow.
* Resubmitting all the analyses from April 28th without using BEAGLE
* Using (20140505) prefix to avoid confusion

### Results of latest batch of analyses

* These are the results of the 20140505 analyses, no marginal likelihood
  estimation, analyses ran only on **noDup** dataset, 1 model of molecular
  evolution per codon position (submitted the wrong template to CIPRES):
* **TODO** ~~check models of molecular evolution as some might actually be the
  old one with one partition per codon.~~ Yes, they are.

   |    Analysis           |  run1  |  run2  | combined
   |-----------------------|--------|--------|---------
   |`noDup_relaxed_coalcst`| barely |   no   |  barely
   |`noDup_relaxed_coalexp`|   no   |   no   |    no
   |`noDup_relaxed_yule`   | barely | barely |   yes
   |`noDup_strict_coalcst` |  no    | barely |   yes
   |`noDup_strict_coalexp` |  no    |  no    |   no
   |`noDup_strict_yule`    |  barely|  XXX   |   XXX    

### conclusions

* The analyses need to be run longer:
  - Not a problem, only takes a couple of hours without BEAGLE on CIPRES for the
    dataset without duplicates. Shouldn't be too much longer for the all dataset.
* I can't do the marginal likelhood estimations on CIPRES, it takes too long if
  I don't use BEAGLE.
  - not much I can do about this, other than running the analyses on HiPerGator,
    but even there, not sure it's going to work. Might be worth submitting 1 job
    to see how long it would take;
* How important are these analyses? I have spent already a lot of time trying to
  get them to work, and I need to refocus on goal:
  > **Goal** Showing that priors on the analysis play a huge role in number of
  species estimated by GMYC and it has nothing to do with phylogenetic
  uncertainty.
  > **Questions for next round of analyses**
  >* What kind of model of molecular evolution is needed? and how does it
  >   influence the results. It seems to really influence mixing, which probably
  >   influences branch length estimation and has consequences for the question
  >   I'm interested in.
  >   * simple (single GTR+G for all codon position)
  >	 * intermediate (use what's suggested by PartitionFinder, but simplify to
  >     not use GTR and replace with HKY+G)
  >  * intermediate (use SDR06 model)
  >	 * complex keep using what's suggested by PartitionFinder.

* These analyses are pretty quick, so let's run one of each on a data that
  didn't perform well (poor mixing) in previous set (`noDup_relaxed_coalexp`)
  and see what happens.

# 20140512 -- BEAST with COI sequences -- models of molecular evolution tests

* submitted ~~5 analyses in a single job on HiPerGator~~ (didn't finish
  overnight as I was hoping) to CIPRES
- `noDup_relaxed_coalexp_3models_estimated` same models as suggested with
  partitionFinder, and estimated base frequencies
- `noDup_relaxed_coalexp_3models_empirical` same as above but with empircal base
  frequencies
- `noDup_relaxed_coalexp_3HKY_esimated` using only HKY models but still with
  Gamma for 1st and 2nd as in partitionFinder
- `noDup_relaxed_coalexp_1GTRG_estimated` a single GTR+G gamma across all
  partitions
- `noDup_relaxed_coalexp_SDRGTR_estimated` using the SDR06 (1,2 + 3 codon
  position, didn't check that my alignment actually conforms to these codon
  positions but should give an idea of it performs regarding mixing) model with
  GTR+G model.

## Results

The models using SDR06, 1 model per partition (with either empirical or
estimated base frequencies) were showing signs of poor mixing on posterior and
prior that seemed correlated with over-parametrization of the model. This
however doesn't seem a problem (nice mixing, and convergence) for the run with a
single GTR+G model or for the model with 3 partitions (HKY+G for positions 1 &
2, but HKY for position 3). ESUs above 350 in these cases. 
 
There are slight differences in the topologies and the posterior probabilities
associated with the nodes, but nothing too crazy. However, some haplotypes with
low support jump around quite a bit and depending on their positions, I see that
they could influence the number of estimated ESUs with GMYC. For now, I'm going
to use HKY+G for each of the first 2 codon positions, and HKY for the last one.

# 20140513 -- BEAST with COI sequences

## Methods

* Based on the results above, I'm going to use HKY+G on the first 2 codon
  positions, and HKY on the third.
* Re-using code from 20140428 to generate list of species in EP+WA+Gala clade,
  constraining monophyly.
* Used same priors as before:
  - tMRCA for EP+WA+Gala: lognormal distribution mean (1.5), sd (0.75) and
  offset (2.5)
  - strict clock and lognormal relaxed clocks with lognormal priors set with
    mean at -4 and SD at 1

## Results
* Nice mixing for all runs. Used these as final. 


**TODO:**
* run diff on input files for both runs to make sure the input file is
exactly identical.



- - - - - - - - - - - 

# Entire complex analyses

## 20131127 entire complex
* The XML file used for the 20131122 analysis wasn't the correct one. There was
  no model of molecular evolution specified, no molecular clock, etc.
* Re-doing everything, using BEAST 2.0.2
* for each model use default values when different from 0, and click on
  estimate. For pInv use 0.2 as starting value, click on estimate
* fix tmca for WA+Gala+EP to 5MYA and 1SD in real space with lognormal prior
  (conservative given that Galapagos and WA are sister species)
* yule prior

## 20131204 entire complex
* The inital run (20131127) is really bad (bad mixing and poor estimations)
  caused probably by sane amount of over-parameterization
* The clockRate estimates were out of the roof (1e25) with the uniform prior
* ESS really bad for COI\_2 site estimates
  * Modifications:
  - changed clockRate prior to Normal [0,1] truncated on [0,+Inf[
  - changed COI\_2 to HKY+G

## 20131211 entire complex
* It's a little better, but still the clockRate estimates and the mutation rate
  estimates are all over the place. Checked "fix mean mutation rate" under the
  "Site Model". (It's also what they do in the tutorial).
* Used lognormal prior on clockRates with mean -2 and SD 1.25, truncated on [0,100]
  arbitrary truncation as I apparently couldn't specify inifinty having had
  saved the file previously with another value for upper limit on the prior
  distribution.

## 20131226 entire complex
* Almost all ESS values > 200 except for prior, and the mutation rates that are
  very low.
* To see if it influence anything, redoing the analysis by fixing the mutation
  rate of 16S to 1, and unchecking "fix mean mutation rate". If that doesn't do
  anything, I'll restart the previous analysis for another week.

## 20140424 entire complex

### Methods
* I might be setting myself up for failure but...
* In site model:
  - used partitions as before, estimating them all but use fixed mean rate:
    - p1     GTR+I+G  (16S + 16Sc + ATP62)
    - p3     HKY+I+G  (ATP61 + COI1)
	- p4     HKY+I+G  (ATP63 + c0775)
	- p5     GTR+I+G  (c0036)
	- p6     HKY+G    (COI2 + H3a2)
	- H3a3   HKY+I    (H3a3 + COI3)
	- p8     HKY      (H3a1 + LSU)
	- p9     HKY+G    (ITS)
* all clocks estimated and set to relaxed lognormal, with exponential prior with
  mean = 1
* set tMRCA prior for EP+WA+Gala:
  - mean 0.55
  - SD 1.5
  - offset 1.705
  - gives: 2.5% at 3.5, Median at 5.14 and 97.5% at 36.2

### Results
* After a week, the run went for 119e6 generations. Convergence almost reached
  (ESS barely > 100 for prior and posterior, good for likelihood at 695) at 45e6
  generations.
* ESS for Tree likelihood for c0036 is low (only 37)
* ESS for TreeHeight, YuleModel and birthRate below 200
* interaction between alpha and pInv for partition 1, 5, and to a lesser extent
  4. Partition p3 is good though. Also needs to remove pInv for H3a_3 (replaced
     with Gamma)
* Looked quickly into ClockstaR to figure out clocks, but it seems that it needs
  quite a bit of prep work to get it working, I'm going to reduce the number of
  clocks to:
  - one for mt loci
  - one for nuc loci (c0036, c0775 and H3a)
  - one for nuc rDNA (ITS + LSU)
  
## 20140502 -- BEAST entire complex

### Methods

* Partition scheme
    - p1     GTR+G    (16S + 16Sc + ATP62)
    - p3     HKY+I+G  (ATP61 + COI1)
	- p4     HKY+G    (ATP63 + c0775)
	- p5     GTR+G    (c0036)
	- p6     HKY+G    (COI2 + H3a2)
	- H3a3   HKY+G    (H3a3 + COI3)
	- p8     HKY      (H3a1 + LSU)
	- p9     HKY+G    (ITS)
* Clocks lognormal relaxed for:
    - mitochondrial loci
	- coding nuclear loci
	- ribosomal nuclear loci
	- all using lognormal with mean -4, SD 1
* Set uniform prior at 5 MY for divergence of EP+WA+Gala group, will use
  divergence of mediterranean clade to see how good that is.

### Side test

While analyzing the data on my laptop it seemed that BEAST didn't use all of the
cores, and that I was using the options `-threads` and `-beast_instances`
wrongly. I read this
[thread](https://groups.google.com/forum/#!topic/beast-users/y9nQ56zygYQ) and
given that I have exactly 8 partitions, I wondered what to do. I set up a really
short run (100,000 generations, sampled every 100) using the full XML file. In
one case used `-beast_instances 1 -threads 8`, in the other `-beast_instances 8
-threads 8`. I used the same seed (10101) and made sure to run both in the same
call as there might be differences in CPU speed depending on the instances used
on HiPerGator. At the end of the run the likelihoods of both runs were
identical, and the second option (specifying 8 instances) was a little faster
(towards the end getting 1h13m/M samples against 1h30m/M samples; total time 552
seconds against 459 seconds). Didn't keep the results.

```
[note added on 2014-06-07 -- a message from Andrew Rambaut on BEAST mailing list
suggests that for best performance it's best to leave -threads unspecified (or
to use -threads -1) which lets BEAST decides the size of thread pool depending
on the dataset (typically it's set to be equal the number of partitions).]
```

### Results after a week

Terrible mixing. Values are all over the place with mutation rate values
completely out of control. Fixing the mean values is not quite enough. Trying to
run it again, this time fixing the substitution rate for p1 (16S+16Sc+ATP6_2) at
1 and estimating all others.

# 20140513 -- BEAST entire complex

## Methods

??? (see `20140519`)

## After 55 hours (73e6 + generations)
Still terrible mixing. Parameters with very low ESS values associated with whole
model (posterior, prior), shape of tree (TreeHeight, Yule model, birthrate),
some mutation rates (p5, p8, p9) and the molecular clocks.

# 20140515 -- BEAST entire complex

## Modifications:
- reduce number of clocks from 3 to 1
- remove partition 9 and merge with partition 8, making it a HKY+G

## After 48 hours -- (62e6 generations)

Still issues with the mixing. Problems seem to be coming from relaxed clocks
parameters, that in turn (or in parallel) affect the values for TreeHeight,
YuleModel and birthRate. A few of the treeLikelihoods are also pretty low.

# 20140519 -- BEAST entire complex

## Modifications:
- to test if the issue with the clock is caused by too much heterogeneity in the
  rates that can't be accommodated by a single clock, use 3 strict clocks
  - one for mtDNA
  - one for protein coding nuc
  - one for rDNA
- removed outgroup (S0213), suggested in BEAST book.

## After 70 hours (74e6 generations)
Overall looks pretty good. Low/correlated ESS values for mutation rate on partition 8 and clock rate for rDNA.
Will finish this analysis and use it for this version of the manuscript. Will however need:
- to do a second run
- will try to use relaxed clock

## After 168 hours (185e6 generations)
- Looks good, all ESU above 200
- downloaded it and trying to summarize it (haven't succeeded yet)
- submitted `run2` limited number of generations to 150e6 and named file with
  prefix `20140526`.

## Summary of methods

### Partitions used

 Partition | Loci               | Model
 ----------|--------------------|-------
 p1        | 16S, 16Sc, ATP6\_2 | GTR+G
 p3        | ATP6\_1, COI\_1    | HKY+G+I
 p4        | ATP6\_3, c0775     | HKY+G
 p5        | c0036              | GTR+G
 p6        | COI\_2, H3a\_2     | HKY+G
 p8        | H3a\_1, ITS, LSU   | HKY
 H3a_3     | H3a\_3, COI\_3     | HKY+G

### Clocks used

* 3 strict clocks
- mitochondrial markers
- protein coding nuclear markers
- rDNA nuclear markers

### Priors

* Yule model
* Exponential priors on all strict clock rates
* Uniform prior set at 5 MYA for divergence of WA+Gala+E species group

## Tree generation

There were too many trees in the output and it crashed `treeannotator`. I had to
use `logcombiner` to reduce the number of trees using a larger burnin (30%).

I re-used a previous command and removed more of the trees when using
`treeannotator` as the option `-burnin 10` was already in there.

```{bash}
./logcombiner -log ~/Documents/Impatiens/20140519.allImpatiens/20140519.impTreeAll.trees \
  -o ~/Documents/Impatiens/20140519.allImpatiens/20140519.impTreeAll_sub.trees -b 30
./treeannotator -heights ca -burnin 10 ~/Documents/Impatiens/20140519.allImpatiens/20140519.impTreeAll_sub.trees \
  ~/Documents/Impatiens/20140519.allImpatiens/20140519.impTree.nex
```

## Results

Tree looks good. 2nd run give similar values for all parameters. Still need to
generate tree. Tree from run1 temporarily used in manuscript. See `20140605` for
follow up.

## Results (combining with 2nd run)

### Step 1 -- removing extra log

Because the first run was not constrained in term of number of generation,
truncating file to 185e6 generations to match length of second run.

```{bash}
grep -n 185000000 20140519.allImpatiens.log 
mv 20140519.allImpatiens.log 20140519.allImpatiens_full.log
head -18842 20140519.allImpatiens_full.log > 20140519.allImpatiens.log
grep -n 185000000 20140519.impTreeAll.trees 
head -18925 20140519.impTreeAll.trees > 20140519.impTree_185e6.trees
echo "END;" >> 20140519.impTree_185e6.trees
```

### Step 2 -- combining the logs

Just bought more RAM, no need to be as conservative with the burnin!

```{bash}
./logcombiner -log ~/Documents/Impatiens/20140519.allImpatiens/run1/20140519.allImpatiens.log -b 20 \
  -log ~/Documents/Impatiens/20140519.allImpatiens/run2/20140519.allImpatiens.log -b 20\
  -o ~/Documents/Impatiens/20140519.allImpatiens/combined_burnin20.log
```

All good, all ESS very high (lowest 440), most are above 1000.

### Step 3 -- combining the trees

I had to do a 20% burnin on the final tree because of memory issues...

```{bash}
./logcombiner -log ~/Documents/Impatiens/20140519.allImpatiens/run1/20140519.impTree_185e6.trees -b 20 \
  -log ~/Documents/Impatiens/20140519.allImpatiens/run2/20140519.impTreeAll.trees -b 20 \
  -o ~/Documents/Impatiens/20140519.allImpatiens/combined_burnin20.trees.nex
./treeannotator -heights ca -burnin 20 \
  ~/Documents/Impatiens/20140519.allImpatiens/combined_burnin20.trees.nex \
  ~/Documents/Impatiens/20140519.allImpatiens/allimpatiens_strict.tree.nex
```

# 20140522 -- RAxML entire complex

Trying to use RAxML to see if I get the same topology.

## Methods
- data file: 20130924 (from all_impatiens, partition finder)
- use partition full partition scheme proposed by partitionFinder (but using GTR+G for all)
`raxmlHPC-HYBRID -s infile -n result -m GTRGAMMA -x 12345 -N 500 -q part -p 12345 -f a`
- set time limit for the job to 5 hours as it uses 64 CPUs...

## Results
- looks fine.
- seems that I get lower support than with Bayesian (might be caused by difference in how BS works)
- from memory slight differences in topology (where the Mediterranean falls for example)

# 20140605 -- BEAST full complex

## Rationale

To test effect of using strict clock on analysis, redoing the analysis using
relaxed clocks for the 3 clocks.

## Methods

- Keeping the subst rate for 16S fixed at 1, letting the other rates being
estimated relatively (as for `20140519`).
- Priors for clocks all set as exponentional distributions with mean at 1.0.
- set run length at 185e6 generations.

## Results

Poor overall mixing mostly driven by poor mixing on the clocks.

# 20140612 -- BEAST full complex

## Rationale

Reducing the number of clocks to see the impact on mixing. First, use strict
clock on rDNA but keep the other two as log-normal.


- - - - - - - - - 

# Old analyses I didn't keep

## 20130513.StarBeast ##
* `20130507.impatiens.nex` for data file
* use partition scheme from partfinder
* 1 clock model for each locus (so 16Sc and 16S combined into 16S_clock)
* tree partitions: mtTree (16S, 16Sc, ATP6, COI), rDNATree (ITS+LSU), c0036Tree
  (c0036), c0775Tree (c0775)
* species sets: impatiens (all except for outgroup S213); EP+Gala+WA
* partition models as in partfinder, except for c0036 which calls for SYM+I+G
  but used GTR+I+G; base frequency set to "estimated" except where partfinder
  calls for "all equal"
* clocks: lognormal relaxed (uncorrelated) for 16S and ITS; strict clock for all other loci;
  according to Lessios COI mutation rate is 0.037 per site /Myr so fix COI clock with this rate
* all tree priors: yule process, piecewise linear & constant root, random starting tree
* priors:
  - tmrca(EP+Gala+WA): lognormal log(mean) 1.5, log(SD) 1
  - {16S,ATP6,c0036,c0775,ITS}\_clock.ucld.mean: lognormal log(mean) -1; log(SD)
    1, initial = 1
  - {H3a,LSU}_clock.ucld.mean: lognormal log(mean) -5; log(SD) 1, initial 1
* 50e6 generation; log and echo every 5e3
* file names: 20130513.imp_starBeast

## 20130515.Beast 

Given the results from 20130513.StarBeast, with low ESSs for many parameters, I
decided to first run the analysis with the estimation of a single tree to see if
I could get better estimates on some parameters of the partition, to help with
the convergence of the *BEAST run.

* use 20130507.impatiens.nex for data file
* partition scheme from partfinder
* 1 clock for each locus (16Sc and 16S combined into 16S\_clock, ITS and LSU
  combined into rDNA\_clock)
* one tree partition (imp\_Tree)
* 2 taxon sets:
  - impatiens: all but outgroup S213
  - EP+Gala+WA
* same partition models as in partfinder, except for c0036 (p5) set up to
  GTR+I+G; base frequencies set to "estimated" except where partfinder calls for
  "all equal"
* clocks: all strict clocks except for 16S and rDNA set to lognormal relaxed;
  all set to "estimate" except for COI_clock set to 0.037
* prior on the tree: speciation: yule process
* priors not set to default:
  - tmrca(EP+Gala+WA): lognormal [1.5, 1]
  - {16S,ATP6,c0036,c0775,H3a,rDNA}\_clock.clock.rate set to lognormal[-1,1]
  - {16S,rDNA}_clock.ucld.mean set to exponential[1]
  - yule.birthRate: uniform [0,100]
* 50e6 generations; log & echo every 5e3
* !!!file names: 20130507.impatiens!!!!


## 20130517.Beast

Inspection of the log files for 20130515 analysis showed relatively good ESS for
most parameters given that the total run was less than 20e6 generations. The
parameters for the lognormal clock were however very low (<20) for most parameters
so edited MANUALLY XML file from 20130515 to change clocks for 16S and rDNA to
strict clock. Set number of generations to 50e6 and run time to 55 hours.


## 20130519.Beast

Inspection of the log files for 20130517 look relatively good even after 30e6+
generations. The low likelihoods were for pInv and alpha for the 16Sc partition.
I therefore removed the pInv parameter to only have GTR+Gamma. The 55 hours
timelimit was not enough. The run aborted at 31.246e6 generations, so increase
wall time at 100 hours.



## 20130524.Beast

Re-run of 20130519 with modification of seed number to check convergence. Results
similar. Combined run provides all ESS > 300. Created file 20130519+24.impatiens.trees
that contains trees sampled by both runs (with both 5000 trees removed as burnin),
and 20130519+24.impatiens.nex that will probably end up being the final tree for the paper.

## 20130531.starBeast

Goal: reduce number of parameters compared to last run; try to see the values of
the parameters estimated to see if some could be put together. Given that it
seems that some parameters (*treeLikelihood) are estimated on all defined
partitions even if they are pooled together for substitution models for
instance, I edited the original nexus files to group partitions that were grouped
together by partition_finder: the 3 codon positions for c0036 and c0775; named
nexus file: 20130531.impatiens.nex.

It doesn't look like the frequencies are really close to each others. Last trial
didn't use strict clock, so going to start there.

* use 20130531.impatiens.nex for data
* use partition scheme for substitution models as indicated by partition finder
* 1 clock for each locus (16Sc and 16S combined into 16S\_clock, ITS and LSU
  combined into rDNA\_clock)
* Trees:
  - mtTree for 16Sc, 16S, ATP6, COI
  - c0036
  - c0775
  - H3a
  - rDNATree for ITS, LSU
* all base frequencies set to "empirical" except where partition finder suggests
  "all equal"
* all strict clocks set to estimate, except for COI which is fixed at 0.037
* all species tree prior set to "yule", population models set to "piecewise
  linear and constant root", ploidy for mtTree set to "mitochondrial"
* Set walltime to 160 hours
* 100e6 generations, log and echo every 1e5*

## 20130604.starBeast

Given the poor mixing on the previous run, attempting to do the run simply on
ESU1 complex: ESU1, WPac, ESU3, gracilis, Hawaii.

* use 20130604.impatiensESU1cplx.nex as input file
* use partition scheme for substitution models as specified in partition finder
  trees:
  - mtTree for 16Sc, 16S, ATP6, COI
  - c0036
  - c0775
  - H3a
  - rDNATree for ITS, LSU
* 1 clock for each locus (16Sc and 16S combined into 16S\_clock, ITS and LSU
  combined into rDNA\\_clock)
* models of molecular evolution as in partition finder, but as in 20130519,
  removed Inv to have only GTR+G for 16Sc (p2).
* set all clocks to strict and estimated, except for COI with fixed rate of
  0.037
* trees,
  - prior with yule, and piecewise lineary & constant root
  - ploidy type for mtTree set to mitochondrial
* Priors all default, except clock.rate set to exponential [.25]
* 100e6 generations, log and echo every 1e5
* wall time to 80 hours (time it took last time with full dataset, so this run
  should last substantially less)

## 20130606.starBeast
Things look pretty good for the last run. Runnning the chain for twice as long.
Changed seed number.

In folder 20130606.starBeast
1. create postburn tree files, to remove the burnin from species.trees. Ran:

    ./logcombiner -trees -burnin 10000000 ~/Documents/Impatiens/20130606.impatiens\_starBeast/20130604.StarBEAST.species.trees ~/Documents/Impatiens/20130606.impatiens\_starBeast/20130604.StarBEAST.species.postburn.trees

and

    ./logcombiner -trees -burnin 16000000 ~/Documents/Impatiens/20130606.impatiens\_starBeast/20130606.StarBEAST.species.trees ~/Documents/Impatiens/20130606.impatiens\_starBeast/20130606.StarBEAST.species.postburn.trees
	
1. combining both sets of trees into one file:

    ./logcombiner -trees ~/Documents/Impatiens/20130606.impatiens\_starBeast/20130604.StarBEAST.species.postburn.trees ~/Documents/Impatiens/20130606.impatiens\_starBeast/20130606.StarBEAST.species.postburn.trees ~/Documents/Impatiens/20130604+06.StarBeast.species.trees

* same thing this log files:

    ./logcombiner -burnin 10000000 ~/Documents/Impatiens/20130606.impatiens\_starBeast/20130604.StarBEAST.log ~/Documents/Impatiens/20130606.impatiens\_starBeast/20130604.StarBEAST.postburn.log

    ./logcombiner -burnin 16000000 ~/Documents/Impatiens/20130606.impatiens\_starBeast/20130606.StarBEAST.log ~/Documents/Impatiens/20130606.impatiens\_starBeast/20130606.StarBEAST.postburn.log

* combining both sets:

    ./logcombiner ~/Documents/Impatiens/20130606.impatiens\_starBeast/20130604.StarBEAST.postburn.log ~/Documents/Impatiens/20130606.impatiens\_starBeast/20130606.StarBEAST.postburn.log ~/Documents/Impatiens/20130606.impatiens\_starBeast/20130604+06.StarBEAST.log

* analysis in tracer of this combined run reveal low ESS for some parameters,
  need to run chain longer. Set up new analysis 20130618.starBeast (see below)

## 20130618.starBeast
* exact same run, just changed seed, and walltime to 120 hours.
  
