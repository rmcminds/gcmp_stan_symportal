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
    vector[NHostNodes] logPDescHost
        = log(PDescHost)';
    vector[NHostNodes] hostDivEdges
        = rep_vector(1.0 / sum(PDescHost), NHostNodes);
    vector[NMicrobeNodes] PDescMicrobe
        = microbeTipAncestorsT[2:,]
          * rep_vector(1.0 / NMicrobeTips, NMicrobeTips);
    row_vector[NMicrobeNodes] logPDescMicrobe
        = log(PDescMicrobe)';
    matrix[NHostNodes, NMicrobeNodes] logPDescBoth
        = rep_matrix(logPDescMicrobe, NHostNodes) + rep_matrix(logPDescHost, NMicrobeNodes);
    row_vector[NMicrobeNodes] microbeDivEdges
        = rep_row_vector(1.0 / sum(PDescMicrobe), NMicrobeNodes);
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestorsTCont
        = microbeAncestorsT;
    matrix[NHostNodes, NHostNodes] hostAncestorsCont
        = hostAncestors;
    int NSubfactorGammas = 0;
    int NSubfactors = sum(NSubPerFactor);
    for(i in 1:NFactors) {
        if(NSubPerFactor[i] > 1) {
            NSubfactorGammas += NSubPerFactor[i] - 1;
        }
    }
    for(i in 1:NMicrobeNodes) {
        microbeAncestorsTCont[i,i] = 1.0 - exp(-microbeAncestorsTCont[i,i]);
    }
    for(i in 1:NHostNodes) {
        hostAncestorsCont[i,i] = 1.0 - exp(-hostAncestorsCont[i,i]);
    }
}
parameters {
    real<lower=0> aveStD;
    vector<lower=0>[2 * NSubfactorGammas] subfactPropsRaw;
    simplex[2 * NFactors + 3] stDProps;
    real<lower=0> aveStDMeta;
    simplex[2] hostDivVsTime;
    simplex[2] microbeDivVsTime;
    simplex[NFactors + 3] metaVarProps;
    vector<lower=0>[NSubfactorGammas] subfactMetaPropsRaw;
    row_vector[NMicrobeNodes] phyloLogVarMultPrev;
    vector[NHostNodes] phyloLogVarMultADiv;
    matrix[NHostNodes, NMicrobeNodes] phyloLogVarMultRaw;
    matrix[NSubfactors, NMicrobeNodes] phyloLogVarMultFacts;
    matrix[NEffects + NHostNodes + 1, NMicrobeNodes + 1] rawMicrobeNodeEffects;
}
transformed parameters {
    simplex[2 * NSubfactors + 3] subfactProps;
    vector<lower=0>[2 * NSubfactors + 3] scales;
    simplex[NSubfactors + 3] subfactMetaProps;
    vector<lower=0>[NSubfactors + 3] metaScales;
    row_vector<lower=0>[NMicrobeNodes] microbeScales;
    vector<lower=0>[NHostNodes] hostScales;
    matrix<lower=0>[NSubfactors, NMicrobeNodes] factScales;
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
                dirichSubFact_lpdf += lmultiply(-NSubPerFactor[i], sum_gamma)
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
                dirichSubFact_lpdf += lmultiply(-NSubPerFactor[i], sum_gamma)
                                      + dirichlet_lpdf(subfactProps[(NSubfactors + normStart):(NSubfactors + normStart - 1 + NSubPerFactor[i])] | rep_vector(1, NSubPerFactor[i]));
                subfactProps[(NSubfactors + normStart):(NSubfactors + normStart - 1 + NSubPerFactor[i])]
                    = subfactProps[(NSubfactors + normStart):(NSubfactors + normStart - 1 + NSubPerFactor[i])]
                      * stDProps[NFactors + i];
                sum_gamma = 1 + sum(segment(subfactMetaPropsRaw,
                                            rawStart,
                                            NSubPerFactor[i] - 1));
                subfactMetaProps[normStart:(normStart - 1 + NSubPerFactor[i])]
                    = append_row(1, segment(subfactMetaPropsRaw, rawStart, NSubPerFactor[i] - 1))
                      / sum_gamma;
                dirichSubFact_lpdf += lmultiply(-NSubPerFactor[i], sum_gamma)
                                      + dirichlet_lpdf(subfactMetaProps[normStart:(normStart - 1 + NSubPerFactor[i])] | rep_vector(1, NSubPerFactor[i]));
                subfactMetaProps[normStart:(normStart - 1 + NSubPerFactor[i])]
                    = subfactMetaProps[normStart:(normStart - 1 + NSubPerFactor[i])]
                      * metaVarProps[i];
                rawStart += NSubPerFactor[i] - 1;
            } else {
                subfactProps[normStart]
                    = stDProps[i];
                subfactProps[NSubfactors + normStart]
                    = stDProps[NFactors + i];
                subfactMetaProps[normStart]
                    = metaVarProps[i];
            }
            normStart += NSubPerFactor[i];
        }
    }
    subfactProps[(2 * NSubfactors + 1):(2 * NSubfactors + 3)]
        = stDProps[(2 * NFactors + 1):(2 * NFactors + 3)];
    subfactMetaProps[(NSubfactors + 1):(NSubfactors + 3)]
        = metaVarProps[(NFactors + 1):(NFactors + 3)];
    scales
        = sqrt((2 * NSubfactors + 3) * subfactProps)
          * aveStD;
    metaScales
        = sqrt((NSubfactors + 3) * subfactMetaProps)
          * aveStDMeta;
    {
        row_vector[NMicrobeNodes] microbeDivPlusTime;
        vector[NHostNodes] hostDivPlusTime;
        matrix[NHostNodes, NMicrobeNodes] bothDivPlusTime;
        row_vector[NMicrobeNodes] logMicrobeVarRaw;
        vector[NHostNodes] logHostVarRaw;
        matrix[NHostNodes, NMicrobeNodes] logPhyloVarRaw;
        matrix[NSubfactors, NMicrobeNodes] logFactVarRaw;
        matrix[NHostNodes, NMicrobeNodes] phyloScales;
        microbeDivPlusTime
            = microbeEdges * microbeDivVsTime[1]
              + microbeDivEdges * microbeDivVsTime[2];
        hostDivPlusTime
            = hostEdges * hostDivVsTime[1]
              + hostDivEdges * hostDivVsTime[2];
        bothDivPlusTime
            = hostDivPlusTime
              * microbeDivPlusTime
              * NMicrobeTips
              * NHostTips
              / (NHostTips + NMicrobeTips);
        logMicrobeVarRaw
            = log(microbeDivPlusTime)
              + (sqrt(microbeDivPlusTime)
                 .* phyloLogVarMultPrev
                 * metaScales[NSubfactors + 1])
                * microbeAncestorsTCont;
        logHostVarRaw
            = log(hostDivPlusTime)
              + hostAncestorsCont
                * (sqrt(hostDivPlusTime)
                   .* phyloLogVarMultADiv
                   * metaScales[NSubfactors + 2]);
        logPhyloVarRaw
            = hostAncestorsCont
                  * (sqrt(bothDivPlusTime)
                     .* phyloLogVarMultRaw
                     * metaScales[NSubfactors + 3])
                  * microbeAncestorsTCont
              + rep_matrix(logHostVarRaw, NMicrobeNodes)
              + rep_matrix(logMicrobeVarRaw, NHostNodes);
        logFactVarRaw
            = metaScales[1:NSubfactors]
              * sqrt(microbeDivPlusTime)
              .* phyloLogVarMultFacts
              * microbeAncestorsTCont
              + rep_matrix(logMicrobeVarRaw, NSubfactors);
        microbeScales
            = scales[2 * NSubfactors + 3]
              * exp((logMicrobeVarRaw
                     - log_sum_exp(logPDescMicrobe + logMicrobeVarRaw))
                    * 0.5);
        hostScales
            = scales[2 * NSubfactors + 1]
              * exp((logHostVarRaw
                     - log_sum_exp(logPDescHost + logHostVarRaw))
                    * 0.5);
        phyloScales
            = scales[2 * NSubfactors + 2]
              * exp((logPhyloVarRaw
                     - log_sum_exp(logPDescBoth + logPhyloVarRaw))
                    * 0.5);
        for(f in 1:NSubfactors) {
            factScales[f,]
                = scales[NSubfactors + f]
                  * exp((logFactVarRaw[f,]
                        - log_sum_exp(logPDescMicrobe + logFactVarRaw[f,]))
                        * 0.5);
        }
        scaledMicrobeNodeEffects[1,1]
            = rawMicrobeNodeEffects[1,1]; //intercept
        scaledMicrobeNodeEffects[1,2:]
            = microbeScales .* rawMicrobeNodeEffects[1,2:]; //microbe prevalence
        scaledMicrobeNodeEffects[2:(NEffects + 1),1]
            = subfactLevelMat
              * segment(scales, 1, NSubfactors)
              .* rawMicrobeNodeEffects[2:(NEffects + 1),1]; //samplewise factor alpha diversity
        scaledMicrobeNodeEffects[(NEffects + 2):,1]
            = hostScales .* rawMicrobeNodeEffects[(NEffects + 2):,1]; //host alpha diversity
        scaledMicrobeNodeEffects[2:(NEffects + 1),2:]
            = subfactLevelMat
              * factScales
              .* rawMicrobeNodeEffects[2:(NEffects + 1),2:]; //samplewise factor microbe interactions
        scaledMicrobeNodeEffects[(NEffects + 2):,2:]
            = phyloScales .* rawMicrobeNodeEffects[(NEffects + 2):,2:]; //host microbe interactions
    }
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
    target += dirichlet_lpdf(metaVarProps | rep_vector(1, NFactors + 3));
    target += std_normal_lpdf(phyloLogVarMultPrev);
    target += std_normal_lpdf(phyloLogVarMultADiv);
    target += std_normal_lpdf(to_vector(phyloLogVarMultRaw));
    target += std_normal_lpdf(to_vector(phyloLogVarMultFacts));
    target += std_normal_lpdf(to_vector(rawMicrobeNodeEffects)[2:]);
    target += logistic_lpdf(rawMicrobeNodeEffects[1,1] | 0,1);
    target += std_normal_lpdf(to_vector(baseLevelMat * rawMicrobeNodeEffects[2:(NEffects + 1),]));
    sampleTipEffects = modelMat * (scaledMicrobeNodeEffects * microbeTipAncestorsT);
    for (n in 1:NObs)
        logit_ratios[n] = sampleTipEffects[sampleNames[n], microbeTipNames[n]];
    target += bernoulli_logit_lpmf(present | logit_ratios);
}
generated quantities {
    row_vector[NMicrobeNodes] microbeNewEdges
        = square(microbeScales / scales[2 * NSubfactors + 3]);
    vector[NHostNodes] hostNewEdges
        = square(hostScales / scales[2 * NSubfactors + 1]);
    matrix[NSumTo0, NMicrobeNodes + 1] baseLevelEffects
        = baseLevelMat * scaledMicrobeNodeEffects[2:(NEffects + 1),];
}
