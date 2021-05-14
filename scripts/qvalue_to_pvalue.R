
# adapted from: https://stats.stackexchange.com/questions/51070/how-can-i-convert-a-q-value-distribution-to-a-p-value-distribution
# use to transform ASD FDR qvalue to pvalues

convert.qval.pval = function(qvalues) {
        # you need to know the estimate of pi0 used to create the q-value
        # that's the maximum q-value (or very, very close to it)
        pi0 = max(qvalues)
    # compute m0, the estimated number of true nulls
    m0 = length(qvalues) * pi0
        # then you multiply each q-value by the proportion of true nulls
        # expected to be under it (the inverse of how you get there from
        # the p-value):
        return(qvalues * rank(qvalues) / m0)
}
