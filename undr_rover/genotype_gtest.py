import math
from collections import defaultdict
from memoized import memoized
import mpmath
from scipy.stats import chi2

DEFAULT_READ_ERROR = 1.0/500.0
DIPLOID_GENOTYPES = ['AA', 'AT', 'AC', 'AG', 'TT', 'TC', 'TG', 'GG', 'GC', 'CC']
HAPLOID_GENOTYPES = ['A', 'T', 'G', 'C']
DNA_BASES = ['A', 'T', 'G', 'C']
DEGREES_OF_FREEDOM = len(DNA_BASES) - 1
 

# Below computes the Chi2 CDF to larger precision than scipy, however it is likely
# to be a bit slower and we probably don't care about the very low probabilities
# which are pretty close to zero anyway. So we opt to use chi2.cdf instead
'''
@memoized
def mpmathcdf(g, df, dps=10):
    mpmath.mp.dps = dps
     
    x,k = mpmath.mpf(g), mpmath.mpf(df)
    cdf = mpmath.gammainc(k/2, 0, x/2, regularized=True)
     
    # floating point precision insufficient, use more precision
    if cdf == 1.0:
        if dps > 4000:
            return cdf # give up after a certain point
        else:
            cdf = mpmathcdf(g,df,dps*2)
    return cdf

def g_test_to_probability(g,df):
    assert g >= 0, g
    return float(1-mpmathcdf(g,df))
'''

def g_test_to_probability(g_score, degrees_of_freedom):
    return 1.0 - chi2.cdf(g_score, degrees_of_freedom)

#from http://blog.mcbryan.co.uk/2013/07/the-g-test-and-python.html

def flnf(f):
    return f * math.log(f) if f > 0.5 else 0

def gtest(a):
    rowtotals = []
    coltotals = [0]*len(a[0])
    cells = []
    for row in a:
        rowtotals.append(sum(row))
        for index,col in enumerate(row):
            coltotals[index] += col
            cells.append(col)
    return 2 * (math.fsum([flnf(x) for x in cells]) + flnf(sum(cells)) - (math.fsum([flnf(x) for x in rowtotals]) + math.fsum([flnf(x) for x in coltotals])))

class Genotype(object):
    def __init__(self, pileup, genotypes, prob_read_error=DEFAULT_READ_ERROR):
        self.prob_read_error = prob_read_error
        self.pileup = pileup
        self.genotypes = genotypes 

    def snv(self):
        observed_dist = defaultdict(int)
        num_observations = len(self.pileup)
    
        for base in self.pileup:
            observed_dist[base] += 1

        best_genotype = None
        max_probability = 0.0

        #scores = []

        for gt in self.genotypes:
            expect_dist = self.expected_distribution(gt, num_observations)
            contingency = [[expect_dist[dna_base], observed_dist[dna_base]] for dna_base in DNA_BASES]
            g_test_score = gtest(contingency)
            probability = g_test_to_probability(g_test_score, DEGREES_OF_FREEDOM)
            #scores.append((g_test_score, probability, gt))
            if probability > max_probability: 
                best_genotype = gt
                max_probability = probability

        if best_genotype is not None:
            #return best_genotype, max_probability, scores
            return best_genotype, max_probability
        else:
            return None

    def expected_distribution(self, genotype, num_observations):
        raise NotImplementedError("expected_distribution called on base class")
    

class Diploid(Genotype):

    def __init__(self, pileup, prob_read_error=DEFAULT_READ_ERROR):
        super(Diploid, self).__init__(pileup, DIPLOID_GENOTYPES, prob_read_error)
        self.homozygous_prob = 1.0 - self.prob_read_error
        self.heterozygous_prob = (3.0 - (2.0 * self.prob_read_error)) / 6.0
        self.error_prob = self.prob_read_error / 3.0

    '''
    Compute the expected frequency distribution of observations in a pilup of a given
    length for a given genotype.
    
    P(X|XX) = 1 - e
    P(X|XY) = (3 - 2e) / 6 
    P(Y|YX) = (3 - 2e) / 6
    P(X|YZ) = e/3
    '''
    def expected_distribution(self, genotype, num_observations):
        gt_1, gt_2 = genotype[0], genotype[1]
        dist = {}
        for base in DNA_BASES:
            if base == gt_1 and base == gt_2:
                dist[base] = num_observations * self.homozygous_prob
            elif base == gt_1 or base == gt_2:
                dist[base] = num_observations * self.heterozygous_prob
            else:
                dist[base] = num_observations * self.error_prob
        return dist

class Haploid(Genotype):

    def __init__(self, pileup, prob_read_error=DEFAULT_READ_ERROR):
        super(Haploid, self).__init__(pileup, HAPLOID_GENOTYPES, prob_read_error)

    '''
    Compute the expected frequency distribution of observations in a pilup of a given
    length for a given genotype.
    
    P(X|X) = 1 - e
    P(X|Y) = e 
    '''
    def expected_distribution(self, genotype, num_observations):
        dist = {}
        for base in DNA_BASES:
            if base == genotype:
                dist[base] = num_observations * (1.0 - self.prob_read_error) 
            else: 
                dist[base] = num_observations * self.prob_read_error
        return dist

def test(pileup, type):
    print(pileup)
    result = type(pileup).snv()
    if result is not None:
        gt, prob, scores  = result
        print("{} {}".format(gt, prob))
        for score, prob, g  in sorted(scores):
            print("{} {} {}".format(g, score, prob))
    else:
        print("No genotype available")

def test_diploid(pileup):
    return test(pileup, Diploid)

def test_haploid(pileup):
    return test(pileup, Haploid)

def run_tests():
    test_diploid("A")
    test_diploid("AT")
    test_diploid("ATGC")
    test_diploid("AA")
    test_diploid("AT")
    test_diploid("AAAA")
    test_diploid("AATTG")
    test_diploid(50 * "A" + 30 * "T")
    test_diploid(50 * "A" + 10 * "T")
    test_diploid(50 * "A" + 5 * "T")
    test_diploid(50 * "A" + 1 * "T")
    test_diploid(500 * "A" + 500 * "G" + 10 * "C")
    
    test_haploid("A")
    test_haploid("AT")
    test_haploid("ATGC")
    test_haploid("AA")
    test_haploid("AT")
    test_haploid("AAAA")
    test_haploid("AATTG")
    test_haploid(50 * "A" + 30 * "T")
    test_haploid(50 * "A" + 10 * "T")
    test_haploid(50 * "A" + 5 * "T")
    test_haploid(50 * "A" + 1 * "T")
    test_haploid(500 * "A" + 500 * "G" + 10 * "C")
