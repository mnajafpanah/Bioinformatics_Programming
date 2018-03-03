
# Gibbs Sampler Motif Finder

Gibbs sampling is the basis behind a general class of algorithms that is a type of local search.  This is one of the algorithm which can find the best conserved motifs among a huge
number of sequences. In the following you can see the pseudocode for the corresponding algorithm which is implemented in this project:
# # # # # # # # # #
    GibbsSampler( Dna , k, t, N)
      randomly select k-mer Motifs = (Motif1, ..., Motif-t) in each string from  Dna BestMotifs <- Motifs
        for j = 1 to N
            i = RANDOM(t)
            Profile <- profile matrix formed from all strings in Motifs except for Motif-i Motif-i <- Profile randomly          
            generated k-mer in the i- th  sequence
            if Score(Motifs) < Score(BestMorifs)
                    BestMorifs <- Motifs
    return BestMorifs
