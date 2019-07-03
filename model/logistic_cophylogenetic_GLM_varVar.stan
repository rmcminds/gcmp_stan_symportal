data {
    int NSamples;
    int NObs;
    int NMicrobeNodes;
    int NMicrobeTips;
    int NFactors;
    int NSubPerFactor[NFactors];
    int NEffects;
    int NHostNodes;
    int NHostTips;
    int present[NObs];
    int sampleNames[NObs];
    int microbeTipNames[NObs];
    real<lower=0> aveStDPriorExpect;
    real<lower=0> aveStDMetaPriorExpect;
    matrix[NEffects, sum(NSubPerFactor)] subfactLevelMat;
    matrix[NSamples, NEffects + NHostNodes + 1] modelMat;
    int NSumTo0;
    matrix[NSumTo0, NEffects] baseLevelMat;
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestorsT;
    matrix[NMicrobeNodes + 1, NMicrobeTips] microbeTipAncestorsT;
    matrix[NHostNodes, NHostNodes] hostAncestors;
    matrix[NHostTips, NHostNodes] hostTipAncestors;
    vector[NHostNodes] hostEdges;
    row_vector[NMicrobeNodes] microbeEdges;
}
transformed data {
    int NSubfactorGammas = 0;
    int NSubfactors = sum(NSubPerFactor);
    for(i in 1:NFactors) {
        if(NSubPerFactor[i] > 1) {
            NSubfactorGammas += NSubPerFactor[i] - 1;
        }
    }
}
parameters {
    real<lower=0> aveStD;
    vector<lower=0>[2 * NSubfactorGammas] subfactPropsRaw;
    simplex[2 * NFactors + 3] stDProps;
    vector<lower=0>[6] stDsMeta;
    simplex[2] codivVsCophyVarProps[3];
    row_vector[NMicrobeNodes] phyloLogVarMultPrev;
    row_vector[NMicrobeNodes] phyloLogVarDivPrev;
    vector[NHostNodes] phyloLogVarMultADiv;
    vector[NHostNodes] phyloLogVarDivADiv;
    matrix[NHostNodes, NMicrobeNodes] phyloLogVarMultRaw;
    matrix[NHostNodes, NMicrobeNodes] phyloLogVarCodivRaw;
    matrix[NEffects + NHostNodes + 1, NMicrobeNodes + 1] rawMicrobeNodeEffects;
}
transformed parameters {
    simplex[2 * NSubfactors + 3] subfactProps;
    vector<lower=0>[2 * NSubfactors + 3] scales;
    row_vector[NMicrobeNodes] microbeRateShifts;
    row_vector<lower=0>[NMicrobeNodes] microbeRates;
    row_vector[NMicrobeNodes] microbeDivergenceVariance;
    row_vector<lower=0>[NMicrobeNodes] microbeDivergence;
    row_vector<lower=0>[NMicrobeNodes] microbeScales;
    vector[NHostNodes] hostRateShifts;
    vector<lower=0>[NHostNodes] hostRates;
    vector[NHostNodes] hostDivergenceVariance;
    vector<lower=0>[NHostNodes] hostDivergence;
    vector<lower=0>[NHostNodes] hostScales;
    matrix[NHostNodes, NMicrobeNodes] cophyloRateShifts;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] cophyloRates;
    matrix[NHostNodes, NMicrobeNodes] coDivergenceVariance;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] coDivergence;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] coPhyloScales;
    matrix[NEffects + NHostNodes + 1, NMicrobeNodes + 1] scaledMicrobeNodeEffects;
    real dirichSubFact_lpdf = 0;
    {
        int rawStart = 1;
        int normStart = 1;
        for (i in 1:NFactors) {
            if(NSubPerFactor[i] > 1) {
                real sum_gamma = 1 + sum(segment(subfactPropsRaw,
                                                 rawStart,
                                                 NSubPerFactor[i] - 1));
                subfactProps[normStart:(normStart - 1 + NSubPerFactor[i])]
                    = append_row(1, segment(subfactPropsRaw, rawStart, NSubPerFactor[i] - 1))
                      / sum_gamma;
                dirichSubFact_lpdf += -NSubPerFactor[i] * log(sum_gamma)
                                      + dirichlet_lpdf(subfactProps[normStart:(normStart - 1 + NSubPerFactor[i])] | rep_vector(1, NSubPerFactor[i]));
                subfactProps[normStart:(normStart - 1 + NSubPerFactor[i])]
                    = subfactProps[normStart:(normStart - 1 + NSubPerFactor[i])]
                      * stDProps[i];
                sum_gamma = 1 + sum(segment(subfactPropsRaw,
                                            NSubfactorGammas + rawStart,
                                            NSubPerFactor[i] - 1));
                subfactProps[(NSubfactors + normStart):(NSubfactors + normStart - 1 + NSubPerFactor[i])]
                    = append_row(1, segment(subfactPropsRaw, NSubfactorGammas + rawStart, NSubPerFactor[i] - 1))
                      / sum_gamma;
                dirichSubFact_lpdf += -NSubPerFactor[i] * log(sum_gamma)
                                      + dirichlet_lpdf(subfactProps[(NSubfactors + normStart):(NSubfactors + normStart - 1 + NSubPerFactor[i])] | rep_vector(1, NSubPerFactor[i]));
                subfactProps[(NSubfactors + normStart):(NSubfactors + normStart - 1 + NSubPerFactor[i])]
                    = subfactProps[(NSubfactors + normStart):(NSubfactors + normStart - 1 + NSubPerFactor[i])]
                      * stDProps[NFactors + i];
                rawStart += NSubPerFactor[i] - 1;
            } else {
                subfactProps[normStart]
                    = stDProps[i];
                subfactProps[NSubfactors + normStart]
                    = stDProps[NFactors + i];
            }
            normStart += NSubPerFactor[i];
        }
    }
    subfactProps[(2 * NSubfactors + 1):(2 * NSubfactors + 3)]
        = stDProps[(2 * NFactors + 1):(2 * NFactors + 3)];
    scales
        = sqrt((2 * NSubfactors + 3) * subfactProps)
          * aveStD;
    microbeRateShifts
        = sqrt(microbeEdges)
          .* phyloLogVarMultPrev
          * stDsMeta[1];
    microbeRates
        = microbeEdges
          .* exp(microbeRateShifts
                 * microbeAncestorsT);
    microbeRates
        = codivVsCophyVarProps[1,1]
          * microbeRates
          / mean(microbeRates * microbeTipAncestorsT[2:,]); //this normalization and the 5 other instances are a major source of computation that could potentially be eliminated to both speed up the model and improve identifiability.
// sum(mR * mTA) = constant //(where constant is the length of the vector so mean would be 1)
// make the parameter the vector equal to mR * mTA, then divide by mTA to get mR
// this vector would need to be constrained to both sum to the constant AND have every element greater than zero (simplex and add len(vector)-1)
    microbeDivergenceVariance
        = sqrt(microbeEdges)
          .* phyloLogVarDivPrev
          * stDsMeta[2];
    microbeDivergence
        = exp(microbeDivergenceVariance
              * microbeAncestorsT);
    microbeDivergence
        = codivVsCophyVarProps[1,2]
          * microbeDivergence
          / mean(microbeDivergence * microbeTipAncestorsT[2:,]);
    microbeScales
        = sqrt(microbeRates + microbeDivergence);
    hostRateShifts
        = sqrt(hostEdges)
          .* phyloLogVarMultADiv
          * stDsMeta[3];
    hostRates
        = hostEdges
          .* exp(hostAncestors
                 * hostRateShifts);
    hostRates
        = codivVsCophyVarProps[2,1]
          * hostRates
          / mean(hostTipAncestors * hostRates);
    hostDivergenceVariance
        = sqrt(hostEdges)
          .* phyloLogVarDivADiv
          * stDsMeta[4];
    hostDivergence
        = exp(hostAncestors
              * hostDivergenceVariance);
    hostDivergence
        = codivVsCophyVarProps[2,2]
          * hostDivergence
          / mean(hostTipAncestors * hostDivergence);
    hostScales
        = scales[2 * NSubfactors + 1]
          * sqrt(hostRates + hostDivergence);
    cophyloRateShifts
        = sqrt(hostEdges * microbeEdges)
          .* phyloLogVarMultRaw
          * stDsMeta[5];
    cophyloRates
        = hostRates * microbeRates
          .* exp(hostAncestors
                 * cophyloRateShifts
                 * microbeAncestorsT);
    cophyloRates
        = codivVsCophyVarProps[3,1]
          * cophyloRates
          / mean(hostTipAncestors
                 * cophyloRates
                 * microbeTipAncestorsT[2:,]);
    coDivergenceVariance
        = sqrt(hostEdges * microbeEdges)
          .* phyloLogVarCodivRaw
          * stDsMeta[6];
    coDivergence
        = exp(hostAncestors
              * coDivergenceVariance
              * microbeAncestorsT);
    coDivergence
        = codivVsCophyVarProps[3,2]
          * coDivergence
          / mean(hostTipAncestors
                 * coDivergence
                 * microbeTipAncestorsT[2:,]);
    coPhyloScales
        = scales[2 * NSubfactors + 2]
          * sqrt(cophyloRates + coDivergence);
    scaledMicrobeNodeEffects
        = append_col(
                append_row(1.0,
                           append_row(subfactLevelMat * segment(scales, 1, NSubfactors),
                                      hostScales)),
                append_row(
                    append_row(
                        scales[2 * NSubfactors + 3],
                        subfactLevelMat * segment(scales, NSubfactors + 1, NSubfactors))
                    * microbeScales,
                    coPhyloScales))
          .* rawMicrobeNodeEffects;
}
model {
    matrix[NSamples, NMicrobeTips] sampleTipEffects;
    vector[NObs] logit_ratios;
    target += exponential_lpdf(aveStD | 1.0 / aveStDPriorExpect);
    target += dirichSubFact_lpdf;
    target += dirichlet_lpdf(stDProps | rep_vector(1, 2 * NFactors + 3));
    for (i in 1:3)
        target += dirichlet_lpdf(codivVsCophyVarProps[i] | rep_vector(1,2));
    target += exponential_lpdf(stDsMeta | 1.0 / aveStDMetaPriorExpect);
    target += std_normal_lpdf(phyloLogVarMultPrev);
    target += std_normal_lpdf(phyloLogVarMultADiv);
    target += std_normal_lpdf(to_vector(phyloLogVarMultRaw));
    target += std_normal_lpdf(phyloLogVarDivPrev);
    target += std_normal_lpdf(phyloLogVarDivADiv);
    target += std_normal_lpdf(to_vector(phyloLogVarCodivRaw));
    target += std_normal_lpdf(to_vector(rawMicrobeNodeEffects)[2:]);
    target += logistic_lpdf(rawMicrobeNodeEffects[1,1] | 0,1);
    target += std_normal_lpdf(to_vector(baseLevelMat * rawMicrobeNodeEffects[2:(NEffects + 1),]));
    sampleTipEffects = modelMat * (scaledMicrobeNodeEffects * microbeTipAncestorsT);
    for (n in 1:NObs)
        logit_ratios[n] = sampleTipEffects[sampleNames[n], microbeTipNames[n]];
    target += bernoulli_logit_lpmf(present | logit_ratios);
}
generated quantities {
    matrix[NSumTo0, NMicrobeNodes + 1] baseLevelEffects;
    baseLevelEffects
        = baseLevelMat * scaledMicrobeNodeEffects[2:(NEffects + 1),];
}
