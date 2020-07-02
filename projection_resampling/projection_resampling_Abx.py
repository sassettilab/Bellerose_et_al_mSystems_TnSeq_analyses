# copyright 2020, Thomas R. Ioerger

import sys,math,numpy
import statsmodels.stats.multitest

#########################################################

class Spreadsheet:
  def __init__(self,filename): # ,key="Id"
    self.keys,self.data,self.headers = [],[],None
    self.rowhash,self.colhash = {},{}
    for line in open(filename,'rU'): # also handle mac files with x0d vs x0a
      if len(line)==0 or line[0]=='#': continue
      w = line.rstrip().split('\t') 
      if self.headers==None: 
        self.headers = w
        for i,h in enumerate(self.headers): 
          if h in self.colhash: print "error: '%s' appears multiple times in first row of '%s' - column headers must be unique" % (h,filename); sys.exit(-1)
          self.colhash[h] = i
        continue
      self.data.append(w)
      key = w[0]
      if key=="": continue
      if key in self.rowhash: print "error: '%s' appears multiple times in first column of '%s' - keys must be unique" % (key,filename); sys.exit(-1)
      self.rowhash[key] = len(self.keys)
      self.keys.append(key) # check if unique?
    self.ncols = max([len(row) for row in self.data])
    self.nrows = len(self.keys)
  def getrow(self,r):
    if r in self.rowhash: r = self.rowhash[r] 
    # check that it is an integer (in range)?
    return self.data[r]
  def getcol(self,c):
    if c in self.colhash: c = self.colhash[c] # otherwise, assume it is an integer
    return [r[c] for r in self.data] # what if not all same length?
  def get(self,r,c):
    if c in self.colhash: c = self.colhash[c] # otherwise, assume it is an integer
    if r in self.rowhash: r = self.rowhash[r]
    return self.data[r][c]


#########################################################
# adapted from resampling() function from Transit, https://transit.readthedocs.io/en/latest/
# in src/pytransit/stat_tools.py

def F_mean_diff_flat(*args, **kwargs):
    A = args[0]
    B = args[1]
    return numpy.mean(B) - numpy.mean(A)

def F_shuffle_flat(*args, **kwargs):
    X = args[0]
    return numpy.random.permutation(X)

def resampling_transit(data1, data2, S=10000, testFunc=F_mean_diff_flat,
            permFunc=F_shuffle_flat, adaptive=False, PC=1,n1=None,n2=None):
    """Does a permutation test on two sets of data.

    Performs the resampling / permutation test given two sets of data using a
    function defining the test statistic and a function defining how to permute
    the data.

    Args:
        ar: List or numpy array with the first set of observations.
        data2: List or numpy array with the second set of observations.
        S: Number of permutation tests (or samples) to obtain.
        testFunc: Function defining the desired test statistic. Should accept
                two lists as arguments. Default is difference in means between
                the observations.
        permFunc: Function defining the way to permute the data. Should accept
                one argument, the combined set of data. Default is random
                shuffle.
        adaptive: Cuts-off resampling early depending on significance.
        n1,n2 = equivalent sample sizes for null distn; len(data1) by default, but user could specify subsets

    Returns:
        Tuple with described values
            - test_obs -- Test statistic of observation.
            - mean1 -- Arithmetic mean of first set of data.
            - mean2 -- Arithmetic mean of second set of data.
            - log2FC -- Normalized log2FC the means.
            - pval_ltail -- Lower tail p-value.
            - pval_utail -- Upper tail p-value.
            - pval_2tail -- Two-tailed p-value.
            - test_sample -- List of samples of the test statistic.
    """

    # - Check input has some data
    assert len(data1) > 0, "Data1 cannot be empty"
    assert len(data2) > 0, "Data2 cannot be empty"

    #FAKE = numpy.array([[1,1,1]])[-1,:]
    #data1 = numpy.concatenate((data1,FAKE))
    #data2 = numpy.concatenate((data2,FAKE))

    count_ltail = 0
    count_utail = 0
    count_2tail = 0

    test_list = []

    # Calculate basic statistics for the input data:
    if n1==None: n1 = len(data1)
    if n2==None: n2 = len(data2)

    mean1 = 0
    if n1 > 0:
        mean1 = numpy.mean(data1)
        #mean1 = numpy.sum(data1)/float(n1)
    mean2 = 0
    if n2 > 0:
        mean2 = numpy.mean(data2)
        #mean2 = numpy.sum(data2)/float(n2)

    if PC>0: log2FC = math.log((mean2+PC)/(mean1+PC),2) # as of 3/5/20
    else:
      # Only adjust log2FC if one of the means is zero
      if mean1 > 0 and mean2 > 0: log2FC = math.log((mean2)/(mean1),2)
      else: log2FC = math.log((mean2+1.0)/(mean1+1.0),2)

 
    # Get stats and info 
    nTAs = 0
    try: test_obs = testFunc(data1, data2)
    except Exception as e:
            print("")
            print("!"*100)
            print("Error: Could not apply test function to input data!")
            print("data1", data1)
            print("data2", data2)
            print("")
            print("\t%s" % e)
            print("!"*100)
            print("")
            return None

    # construct pooled list of counts
    perm = numpy.zeros(len(data1)+len(data2))
    perm[:len(data1)] = data1
    perm[len(data1):] = data2

    count_ltail = 0
    count_utail = 0
    count_2tail = 0
    test_list = []
    s_performed = 0
    for s in range(S):
        if len(perm) >0:
            perm = permFunc(perm)
            test_sample = testFunc(perm[:n1], perm[n1:n1+n2]) # use equivalent sizes for null distribution
        else:
            test_sample = 0

        test_list.append(test_sample)
        if test_sample <= test_obs: count_ltail+=1
        if test_sample >= test_obs: count_utail+=1
        if abs(test_sample) >= abs(test_obs): count_2tail+=1

        s_performed+=1
        if adaptive:
            if s_performed == round(S*0.01) or s_performed == round(S*0.1) or s_performed == round(S*1):
                    if count_2tail >= round(S*0.01*0.10):
                        break



    pval_ltail = count_ltail/float(s_performed)
    pval_utail = count_utail/float(s_performed)
    pval_2tail = count_2tail/float(s_performed)

    return (test_obs, mean1, mean2, log2FC, pval_ltail, pval_utail,  pval_2tail, test_list)

#########################################################

# returns list of (start,end,Rv,gene,dir,descr)
# convert to 0-based indexing of coords

def read_genes(fname,descriptions=False):
  genes = []
  for line in open(fname):
    w = line.split('\t')
    data = [int(w[1])-1,int(w[2])-1,w[8],w[7],w[3]]
    if descriptions==True: data.append(w[0])
    genes.append(data)
  return genes

def probdist(vals):
  #temp = [math.exp(alpha*x) for x in vals]
  temp = [x*x for x in vals]
  tot = sum(temp)
  normed = [x/tot for x in temp]
  return normed

#########################################################
# main
#########################################################

#### read input files

if len(sys.argv)<6:
  print "usage: python projection_resampling_Abx.py <counts.txt> <samples_metadata.txt> <genome.prot_table> <loadings.txt> <genes_list.txt> [-v]"
  sys.exit()

VERBOSE = False
if "-v" in sys.argv: VERBOSE = True

print "# command: python",' '.join(sys.argv)

data = Spreadsheet(sys.argv[1])
coords = [int(x) for x in data.getcol(0)]
colnames = data.headers[3:]
counts = []

samples = Spreadsheet(sys.argv[2])

for row in range(data.nrows):
  counts.append([int(float(x)) for x in data.data[row][3:]]) # user is responsible for normalizing input counts
colconds = [samples.get(x,"Condition") for x in colnames]

genes = read_genes(sys.argv[3])

skip = 1
conditions,mat = [],[]
varimax_headers = None
for line in open(sys.argv[4]):
  w = line.rstrip().split('\t')
  if skip>=1: skip -= 1; varimax_headers = w; continue
  conditions.append(w[0])
  mat.append([float(x) for x in w[1:]])

loadings = [] # columns of mat, i.e. transpose
ndims,dims = len(mat[0]),[]
for i in range(ndims):
  loadings.append([x[i] for x in mat])

if VERBOSE:
 #print "rotation matrix:"
 #for i,cond in enumerate(conditions):
 #  print "%-20s" % cond,
 #  for j in range(ndims): print "%8.5f" % mat[i][j],
 #  print 
 print conditions
 print "component loadings:"
 for i in range(ndims):
   print "V%s" % (1+i),
   for x in loadings[i]: print "%8.5f" % x,
   print 
 print "sampling weights:"
 for i in range(ndims):
   print "V%s" % (1+i),
   for x in probdist(loadings[i]): print "%8.5f" % x,
   print

anova = {}
skip = 1
for line in open(sys.argv[5]): # temp_anova.txt
  if skip>0: skip -= 1; continue
  w = line.rstrip().split('\t')
  if w[-1]=='nan': padj = 1
  else: padj = float(w[-1])  
  anova[w[0]] = padj


### do projection and resampling

results = [] # list of [rv,dim,log2FC,pval]
gene_data = []
for (start,end,rv,gene,strand) in genes:
  if rv not in anova: 
    if VERBOSE: sys.stderr.write("warning: %s not found in %s\n" % (rv,sys.argv[5]))
    continue
  if anova[rv]>=0.05: continue
  sys.stderr.write("%s\n" % rv)
  localCounts = {} # for each condition, observations over all TA sites and all replicates pooled
  localCountsPerSite = {} # list of counts over replicates, indexed by (cond,coord)
  for cond in conditions: localCounts[cond] = []
  sites = []
  for i in range(len(coords)):
    coord = coords[i]
    if coord>=start and coord<=end:
      sites.append(coord)
      for j in range(len(counts[i])):
        cond,cnt = colconds[j],counts[i][j]
        localCounts[cond].append(cnt)
        if (cond,coord) not in localCountsPerSite: localCountsPerSite[(cond,coord)] = []
        localCountsPerSite[(cond,coord)].append(cnt)
  condMeansPerSite = {}
  for cond in conditions: 
    for site in sites:
      condMeansPerSite[(cond,site)] = numpy.mean(localCountsPerSite[(cond,site)])

  vals = [rv,gene,len(sites)]
  means = []
  if VERBOSE: print '\t'.join([str(x) for x in vals])
  for cond in conditions:
    mean = numpy.mean(localCounts[cond])
    means.append(numpy.round(mean,1))
    if VERBOSE: print "%-20s %5.1f" % (cond,mean),localCounts[cond]

  # do projection of counts onto varimax dimensions

  totObs = sum([len(x) for x in localCounts.values()])
  n2 = int(totObs/float(ndims)) # equivalent sample size for individual PC or Vmax dim (tot/k)
  n1 = totObs-n2 # equivalent sample size for all the other PCs or dims (tot*(k-1)/k)

  nTA,ncond = len(sites),len(conditions)
  totReps = totObs/float(nTA)
  equivSampSize = int(totObs/float(ndims))
  if VERBOSE: print "nTA=%s, totObs=%s, totReps=%s, ncond=%s, meanReps=%0.1f, ndims=%s, virtReps=%0.1f, equivSampSize=%s" % (nTA,totObs,totReps,ncond,totReps/float(ncond),ndims,equivSampSize/float(nTA),equivSampSize)

  # take weighted average of counts at each TA site for each PC or Varimax dim
  projectedcounts = []
  SCALE = len(conditions)# /float(ndims)
  for i in range(ndims):
    probs = probdist(loadings[i]) # note: loadings has been transposed from file
    newcounts = []
    for site in sites:
      for k,cond in enumerate(conditions):
        newcounts.append(probs[k]*condMeansPerSite[(cond,site)]*SCALE) 
    projectedcounts.append(newcounts)
  #grandmean = numpy.mean([numpy.mean(x) for x in projectedcounts])
  projmeans = [numpy.mean(x) for x in projectedcounts]
  projmeans.sort()
  temp = len(projmeans)
  median = (projmeans[temp/2]+projmeans[temp/2+1])/2.0

  # do resampling of each condition against all the others

  if VERBOSE: print "projection:"
  vals = [rv,gene,len(sites)]+means
  gene_data.append(vals)
  for i in range(ndims):    
    data2 = numpy.array(projectedcounts[i]) # counts for the condition of interest
    data1 = numpy.concatenate([projectedcounts[x] for x in range(len(projectedcounts)) if x!=i]) # all others as ref

    # args based on pytransit/src/resampling.py: data1, data2, S=self.samples, testFunc=stat_tools.F_mean_diff_dict, permFunc=stat_tools.F_shuffle_dict_libraries, adaptive=self.adaptive, lib_str1=self.ctrl_lib_str, lib_str2=self.exp_lib_str)
    # test_obj is test statistic (diff of means); testlist is samples from null distn
    # log2FC represents log2(mean2/mean1) with pseudo-counts of 1

    (test_obs, mean1, mean2, log2FC, pval_ltail, pval_utail,  pval_2tail, testlist) =  resampling_transit(data1, data2, S=10000, adaptive=True) 
    dim = "V%s" % (i+1)
    log2FC = math.log((mean2+1)/(median+1),2) # redefine based on median of means, instead of Vdim vs others from resampling()
    results.append([rv,dim,mean2,numpy.round(log2FC,1),pval_2tail]) 

    vals2 = ["#resamp","V%s" % (i+1)]+[numpy.round(x,3) for x in [mean2,mean1,log2FC]]+[pval_2tail]
    if VERBOSE: print '\t'.join([str(x) for x in vals2])
  if VERBOSE: print "median of projected counts on Vdims:",round(median,1)


### print out results

pvals = [x[-1] for x in results]
(rejected,qvals) = statsmodels.stats.multitest.fdrcorrection(pvals) # compute adjusted p-values (q-vals)

vals = "ORF gene numTAs".split()+conditions+["%s_mean" % x for x in varimax_headers]+["%s_LFC" % x for x in varimax_headers]+["%s_padj" % x for x in varimax_headers]
print '\t'.join(vals)

for vals in gene_data:
  rv = vals[0]
  lfcs,means,rvqvals = [],[],[]
  for i in range(len(results)):
    if results[i][0]==rv: 
      means.append(results[i][2])
      lfcs.append(results[i][3])
      rvqvals.append(numpy.round(qvals[i],6))
  vals += [round(x,1) for x in means]+[round(x,3) for x in lfcs]+rvqvals
  print '\t'.join([str(x) for x in vals])
  
