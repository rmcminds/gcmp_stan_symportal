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
    row_vector[NHostNodes] PDescHost
        = rep_row_vector(1.0 / NHostTips, NHostTips)
          * hostTipAncestors;
    vector[NHostNodes] hostDivEdges
        = rep_vector(1.0 / sum(PDescHost), NHostNodes);
    vector[NMicrobeNodes] PDescMicrobe
        = microbeTipAncestorsT[2:,]
          * rep_vector(1.0 / NMicrobeTips, NMicrobeTips);
    row_vector[NMicrobeNodes] microbeDivEdges
        = rep_row_vector(1.0 / sum(PDescMicrobe), NMicrobeNodes);
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
    real<lower=0> aveStDMeta;
    simplex[2] hostDivVsTime;
    simplex[2] microbeDivVsTime;
    simplex[3] metaVarProps;
    row_vector[NMicrobeNodes] phyloLogVarMultPrev;
    vector[NHostNodes] phyloLogVarMultADiv;
    matrix[NHostNodes, NMicrobeNodes] phyloLogVarMultRaw;
    matrix[NEffects + NHostNodes + 1, NMicrobeNodes + 1] rawMicrobeNodeEffects;
}
transformed parameters {
    simplex[2 * NSubfactors + 3] subfactProps;
    vector<lower=0>[2 * NSubfactors + 3] scales;
    vector<lower=0>[3] metaScales;
    row_vector<lower=0>[NMicrobeNodes] microbeDivPlusTime;
    row_vector<lower=0>[NMicrobeNodes] microbeVarRaw;
    row_vector<lower=0>[NMicrobeNodes] microbeScales;
    vector<lower=0>[NHostNodes] hostDivPlusTime;
    vector<lower=0>[NHostNodes] hostVarRaw;
    vector<lower=0>[NHostNodes] hostScales;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] phyloVarRaw;
    matrix<lower=0>[NHostNodes, NMicrobeNodes] phyloScales;
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
    metaScales
        = sqrt(3 * metaVarProps)
          * aveStDMeta;
    microbeDivPlusTime
        = microbeEdges * microbeDivVsTime[1]
          + microbeDivEdges * microbeDivVsTime[2];
    microbeVarRaw
        = microbeDivPlusTime
          .* exp((microbeDivPlusTime
                  .* phyloLogVarMultPrev
                  * metaScales[1])
                 * microbeAncestorsT);
    microbeScales
        = sqrt(microbeVarRaw
               / (microbeVarRaw * PDescMicrobe));
    hostDivPlusTime
        = hostEdges * hostDivVsTime[1]
          + hostDivEdges * hostDivVsTime[2];
    hostVarRaw
        = hostDivPlusTime
          .* exp(hostAncestors
                 * (hostDivPlusTime
                    .* phyloLogVarMultADiv
                    * metaScales[2]));
    hostScales
        = scales[2 * NSubfactors + 1]
          * sqrt(hostVarRaw
                 / (PDescHost * hostVarRaw));
    phyloVarRaw
        = exp(hostAncestors
              * (phyloLogVarMultRaw
                 * metaScales[3])
              * microbeAncestorsT)
          .* (hostVarRaw * microbeVarRaw);
    phyloScales
        = scales[2 * NSubfactors + 2]
          * sqrt(phyloVarRaw
                 / (PDescHost
                    * phyloVarRaw
                    * PDescMicrobe));
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
                    phyloScales))
          .* rawMicrobeNodeEffects;
}
model {
    matrix[NSamples, NMicrobeTips] sampleTipEffects;
    vector[NObs] logit_ratios;
    target += exponential_lpdf(aveStD | 1.0 / aveStDPriorExpect);
    target += dirichSubFact_lpdf;
    target += dirichlet_lpdf(stDProps | rep_vector(1, 2 * NFactors + 3));
    target += exponential_lpdf(aveStDMeta | 1.0 / aveStDMetaPriorExpect);
    target += dirichlet_lpdf(hostDivVsTime | rep_vector(1, 2));
    target += dirichlet_lpdf(microbeDivVsTime | rep_vector(1, 2));
    target += dirichlet_lpdf(metaVarProps | rep_vector(1, 3));
    target += std_normal_lpdf(phyloLogVarMultPrev);
    target += std_normal_lpdf(phyloLogVarMultADiv);
    target += std_normal_lpdf(to_vector(phyloLogVarMultRaw));
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
