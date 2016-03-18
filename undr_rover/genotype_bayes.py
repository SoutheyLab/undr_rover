from math import log, exp

'''
Genotyping from variants.

Algorithm based on:

    SNP detection for massively parallel whole-genome resequencing
    Ruiqiang Li, et al.
    http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2694485/

    There are 10 possible genotypes for diploid cell:

    GT = {AA, AT, AC, AG, TT, TC, TG, GG, GC, CC}

    probability of GT_i given some pileup data

    P(GT_i | D) = P(D | GT_i) * P(GT_i)
                  ---------------------
                          P(D)

    The genotype with the highest posterior probability P(GT_i | D)
    is taken as the consensus. 

    Assume uniform distribution of genotypes:

    P(GT_i) = 0.1

    P(D) = Sum(n = 1..10, P(GT_n) * P(D | GT_n))

    Liklihood:

    For a set of n observed alleles in pileup: D = {d_1, d_2 .. d_n}

    P(D | GT_i) = multinomial.probability(D, GT_i)

    count(b, D) = |{x, x <- D, x = b}|
    counts(D) = <count(A, D), count(T, D), count(G, D), count(C, D)>

    prob(b, GT) = P(b | GT)
    probs(GT) = <prob(A, GT), prob(T, GT), prob(G, GT), prob(C, GT)>

    PMF is the probability mass function of the multinomial distribution

    multinomial.probability(D, GT_i) = PMF(counts(D), probs(GT_i), n)

    PMF(<c1, c2 .. ck>, <p1, p2, .. pk>, n) = (n!/c1!c2!..ck!) * p1^c1 * p2^c2 .. pk^ck

    Notice that the PMF is based on the product of:
        - (n!/c1!c2!..ck!)
        - p1^c1 * p2^c2 .. pk^ck
    The first component is independent of the genotype, therefore we can factor it
    out. Let's write the PMF as M*N

    We observe that Bayes' rule can be written as:

    P(GT_i | D) = P(D | GT_i) * P(GT_i)
                  ---------------------
                  Sum(n = 1..10, P(GT_n) x P(D | GT_n))

                = M * N * P(GT_i) 
                  ---------------
                  (Mp + Mp2 + ... Mp10)

   So the Ms cancel out.

   Therefore it is okay to compute:

   P(D | GT_i) = p1^c1 * p2^c2 .. pk^ck
               = P(d_1|GT_i) * P(d_2|GT_i) * ... * P(d_n|GT_i)

   However, this value can become VERY small for long pileups, with a risk of
   underflow. So we should compute it in log space instead:

               = log(P(d_1|GT_i) + log(P(d_2|GT_i)) + ... + log(P(d_n|GT_i))

    We use a simplified single base conditional probabilty compared to the
    paper (diploid case):

    Assuming an input base error rate of E:

    P(d=X | GT=XX) = 1 - E
    P(d=X | GT=XY) = (3 - 2e) / 6
    P(d=X | GT=YX) = (3 - 2e) / 6
    P(d=X | GT=YZ) = e/3

    We observe that the "normalising factor" P(D) is independent of the genotype.
    It is also difficult to compute in log-space because it is a summation
    of terms. We can simply drop the computation of this term and
    assume that the result is now a score rather than a probability:

    P(GT | D) ~ P(G) * P(D|GT)

    in log space:

    log(P(GT | D)) ~ log(P(G)) + log(P(D | GT))

    we compute the above score for each genotype and take the maximum as the result.

'''

DEFAULT_READ_ERROR = 1.0/500.0
UNIFORM_GENOTYPE_PROBABILITY = 0.1
LOG_UNIFORM_GENOTYPE_PROBABILITY = log(UNIFORM_GENOTYPE_PROBABILITY)
DIPLOID_GENOTYPES = ['AA', 'AT', 'AC', 'AG', 'TT', 'TC', 'TG', 'GG', 'GC', 'CC']

class Genotype(object):

    def __init__(self, pileup, prob_read_error=DEFAULT_READ_ERROR):
        self.prob_read_error = prob_read_error
        self.pileup = pileup

    # assumes all bases in pileup are uppercase and only one of A,T,G,C
    def snv(self):
        scores = []
    
        for gt in self.genotypes:
            log_score = self.log_prob_gt_given_pileup(gt) 
            scores.append((log_score, gt))
    
        scores.sort()

        normaliser = 0.0

        for log_score, gt in scores:
            normaliser += exp(log_score)

        if len(scores) > 0: 
            best_log_score, best_gt = scores[-1]
            best_score = exp(best_log_score) 
            if normaliser != 0:
                best_prob = best_score / normaliser
            else:
                best_prob = None
            return best_gt, best_prob 
        else:
            return None
    
    def log_prob_gt_given_pileup(self, genotype):
        '''P(gt|d) = P(d|gt) * P(g)
           log(P(gt|d)) = log(P(d|gt)) + log(P(g))
        '''
        return self.log_prob_pileup_given_gt(genotype) + LOG_UNIFORM_GENOTYPE_PROBABILITY
    
    def log_prob_pileup_given_gt(self, genotype):
        '''P(d|gt)      = P(d_1|gt) * P(d_2|gt) ... * P(d_n|gt)
           log(P(d|gt)) = log(P(d_1|gt)) + log(P(d_2|gt)) + ... + log(P(d_n|gt))
        '''
        result = 0.0 
        for dna_base in self.pileup:
            result += self.log_prob_base_given_gt(dna_base, genotype)
        return result

    def log_prob_base_given_gt(self, dna_base, genotype):
        raise NotImplementedError("log_prob_base_given_gt called on base class")

class Diploid(Genotype):

    def __init__(self, pileup, prob_read_error=DEFAULT_READ_ERROR):
        super(Diploid, self).__init__(pileup, prob_read_error)
        self.genotypes = DIPLOID_GENOTYPES
    
    def log_prob_base_given_gt(self, dna_base, genotype):
        '''log(P(d|gt))
           
           P(X|XX) = 1 - e
           P(X|XY) = ((1 - e) / 2) + (e/3)/2
                   = (3 - 2e) / 6
           P(X|YX) = (3 - 2e) / 6
           P(X|YZ) = e/3
        '''
        genotype_1 = genotype[0]
        genotype_2 = genotype[1]
    
        if dna_base == genotype_1 and dna_base == genotype_2:
            result = 1.0 - self.prob_read_error
        elif dna_base == genotype_1 or dna_base == genotype_2:
            result = (3.0 - (2.0 * self.prob_read_error)) / 6.0 
        else:
            result = self.prob_read_error / 3.0
    
        return log(result)
