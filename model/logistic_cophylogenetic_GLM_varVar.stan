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
    vector[NHostNodes] logHostEdges
        = log(hostEdges);
    row_vector[NMicrobeNodes] logMicrobeEdges
        = log(microbeEdges);
    vector[NHostNodes] sqrtHostEdges
        = sqrt(hostEdges);
    row_vector[NMicrobeNodes] sqrtMicrobeEdges
        = sqrt(microbeEdges);
    vector[NHostNodes] hostIntegral
        = sqrtHostEdges * 0.5; // constant growth of log variance over time of branch means average log variance along branch is 1/2 of the final effect size
    row_vector[NMicrobeNodes] microbeIntegral
        = sqrtMicrobeEdges * 0.5;
    matrix[NHostNodes, NHostNodes] hostAncestorsScaled
        = diag_post_multiply(hostAncestors, sqrtHostEdges);
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestorsTScaled
        = diag_pre_multiply(sqrtMicrobeEdges, microbeAncestorsT);
    vector[NHostNodes] logPDescHost
        = log(rep_row_vector(1.0 / NHostTips, NHostTips)
              * hostTipAncestors)';
    row_vector[NMicrobeNodes] logPDescMicrobe
        = log(microbeTipAncestorsT[2:,]
              * rep_vector(1.0 / NMicrobeTips, NMicrobeTips))';
    matrix[NHostNodes, NMicrobeNodes] logPDescBoth
        = rep_matrix(logPDescMicrobe, NHostNodes) + rep_matrix(logPDescHost, NMicrobeNodes);
    real tipPerEdgeHost
        = NHostTips / sum(hostTipAncestors);
    real tipPerEdgeMicrobe
        = NMicrobeTips / sum(microbeTipAncestorsT[2:,]);
    real logTipPerEdgeHost
        = log(tipPerEdgeHost);
    real logTipPerEdgeMicrobe
        = log(tipPerEdgeMicrobe);
    real sqrtTipPerEdgeHost
        = sqrt(tipPerEdgeHost);
    real sqrtTipPerEdgeMicrobe
        = sqrt(tipPerEdgeMicrobe);
    matrix[NHostNodes, NHostNodes] hostAncestorsDivScaled
        = hostAncestors * sqrtTipPerEdgeHost;
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestorsTDivScaled
        = microbeAncestorsT * sqrtTipPerEdgeMicrobe;
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestorsTScaledStart
        = microbeAncestorsTScaled;
    matrix[NHostNodes, NHostNodes] hostAncestorsScaledStart
        = hostAncestorsScaled;
    int NSubfactorGammas = 0;
    int NSubfactors = sum(NSubPerFactor);
    for(i in 1:NMicrobeNodes) {
        microbeAncestorsTScaledStart[i,i] = 0;
    }
    for(i in 1:NHostNodes) {
        hostAncestorsScaledStart[i,i] = 0;
    }
    for(i in 1:NFactors) {
        if(NSubPerFactor[i] > 1) {
            NSubfactorGammas += NSubPerFactor[i] - 1;
        }
    }
}
parameters {
    simplex[2] hostDivVsTime;
    simplex[2] microbeDivVsTime;
    real<lower=0> aveStD;
    vector<lower=0>[2 * NSubfactorGammas] subfactPropsRaw;
    simplex[2 * NFactors + 3] stDProps;
    real<lower=0> aveStDMeta;
    simplex[NFactors + 3] metaVarProps;
    vector<lower=0>[NSubfactorGammas] subfactMetaPropsRaw;
    row_vector[NMicrobeNodes] phyloLogVarMultPrev;
    vector[NHostNodes] phyloLogVarMultADiv;
    matrix[NHostNodes, NMicrobeNodes] phyloLogVarMultRaw;
    matrix[NSubfactors, NMicrobeNodes] phyloLogVarMultFacts;
    matrix[NEffects + NHostNodes + 1, NMicrobeNodes + 1] rawMicrobeNodeEffects;
    vector<lower=-1, upper=1>[3] varEffectCor;
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
    vector<lower=0, upper=1>[3] varEffectChol2 = sqrt(1 - square(varEffectCor));
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
        = sqrt((NHostNodes + NSubfactors + 1) * (NMicrobeNodes + 1) * subfactProps)
          * aveStD;
    metaScales
        = sqrt(((NHostNodes + 1) * (NMicrobeNodes + 1) + NSubfactors * NMicrobeNodes) * subfactMetaProps)
          * aveStDMeta;
    {
        vector[2] logHostDivVsTime;
        vector[2] logMicrobeDivVsTime;
        row_vector[NMicrobeNodes] microbeIntegralScaled;
        vector[NHostNodes] hostIntegralScaled;
        matrix[NMicrobeNodes, NMicrobeNodes] microbeStart;
        matrix[NHostNodes, NHostNodes] hostStart;
        row_vector[NMicrobeNodes] logMicrobeVarRaw;
        row_vector[NMicrobeNodes] microbeMetaVar;
        vector[NHostNodes] logHostVarRaw;
        vector[NHostNodes] hostMetaVar;
        matrix[NHostNodes, NMicrobeNodes] logPhyloVarRaw;
        matrix[NHostNodes, NMicrobeNodes] phyloMetaVar;
        matrix[NHostNodes, NMicrobeNodes] phyloMetaVarP1;
        matrix[NHostNodes, NMicrobeNodes] phyloMetaVarP1Scaled;
        matrix[NHostNodes, NMicrobeNodes] phyloMetaVarP2Scaled;
        matrix[NSubfactors, NMicrobeNodes] logFactVarRaw;
        matrix[NSubfactors, NMicrobeNodes] factMetaVar;
        row_vector[NMicrobeNodes] logMicrobeTimeVarRaw;
        vector[NHostNodes] logHostTimeVarRaw;
        matrix[NHostNodes, NMicrobeNodes] logPhyloTimeVarRaw;
        matrix[NSubfactors, NMicrobeNodes] logFactTimeVarRaw;
        row_vector[NMicrobeNodes] logMicrobeDivVarRaw;
        vector[NHostNodes] logHostDivVarRaw;
        matrix[NHostNodes, NMicrobeNodes] logPhyloDivTimeVarRaw;
        matrix[NHostNodes, NMicrobeNodes] logPhyloTimeDivVarRaw;
        matrix[NHostNodes, NMicrobeNodes] logPhyloDivVarRaw;
        matrix[NSubfactors, NMicrobeNodes] logFactDivVarRaw;
        matrix[NHostNodes, NMicrobeNodes] phyloScales;
        logHostDivVsTime
            = log(hostDivVsTime);
        logMicrobeDivVsTime
            = log(microbeDivVsTime);
        microbeIntegralScaled
            = microbeIntegral * sqrt(microbeDivVsTime[1]);
        hostIntegralScaled
            = hostIntegral * sqrt(hostDivVsTime[1]);
        microbeStart
            = microbeAncestorsTScaledStart
              * sqrt(microbeDivVsTime[1])
              + microbeAncestorsTDivScaled
                * sqrt(microbeDivVsTime[2]);
        hostStart
            = hostAncestorsScaledStart
              * sqrt(hostDivVsTime[1])
              + hostAncestorsDivScaled
                * sqrt(hostDivVsTime[2]);
        logMicrobeVarRaw
            = phyloLogVarMultPrev * metaScales[NSubfactors + 1];
        microbeMetaVar
            = logMicrobeVarRaw
              * microbeStart;
        logMicrobeTimeVarRaw
            = logMicrobeEdges
              + microbeMetaVar
              + logMicrobeVarRaw
                .* microbeIntegralScaled; // includes instantaneous changes in variance and gradual changes of self and ancestors
        logMicrobeTimeVarRaw
            = logMicrobeTimeVarRaw
              - log_sum_exp(logPDescMicrobe + logMicrobeTimeVarRaw);
        logMicrobeDivVarRaw
            = logTipPerEdgeMicrobe
              + microbeMetaVar; // includes instantaneous changes in variance and only ancestral gradual changes
        logMicrobeDivVarRaw
            = logMicrobeDivVarRaw
              - log_sum_exp(logPDescMicrobe + logMicrobeDivVarRaw);
        logHostVarRaw
            = phyloLogVarMultADiv * metaScales[NSubfactors + 2];
        hostMetaVar
            = hostStart
              * logHostVarRaw;
        logHostTimeVarRaw
            = logHostEdges
              + hostMetaVar
              + logHostVarRaw
                .* hostIntegralScaled;
        logHostTimeVarRaw
            = logHostTimeVarRaw
              - log_sum_exp(logPDescHost + logHostTimeVarRaw);
        logHostDivVarRaw
            = logTipPerEdgeHost
              + hostMetaVar;
        logHostDivVarRaw
            = logHostDivVarRaw
              - log_sum_exp(logPDescHost + logHostDivVarRaw);
        logPhyloVarRaw
            = phyloLogVarMultRaw * metaScales[NSubfactors + 3];
        phyloMetaVarP1
            = hostStart
              * logPhyloVarRaw;
        phyloMetaVarP1Scaled
            = diag_post_multiply(phyloMetaVarP1, microbeIntegralScaled)
              + rep_matrix(logMicrobeTimeVarRaw, NHostNodes);
        phyloMetaVarP2Scaled
            = diag_pre_multiply(hostIntegralScaled, logPhyloVarRaw * microbeStart)
              + rep_matrix(logHostTimeVarRaw, NMicrobeNodes);
        phyloMetaVar
            = phyloMetaVarP1
              * microbeStart;
        logPhyloTimeVarRaw
            = phyloMetaVar
              + phyloMetaVarP1Scaled
              + phyloMetaVarP2Scaled
              + hostIntegralScaled * microbeIntegralScaled;
        logPhyloTimeVarRaw
            = logPhyloTimeVarRaw
              - log_sum_exp(logPDescBoth + logPhyloTimeVarRaw);
        logPhyloTimeDivVarRaw
            = rep_matrix(logHostDivVarRaw, NMicrobeNodes)
              + phyloMetaVar
              + phyloMetaVarP1Scaled;
        logPhyloTimeDivVarRaw
            = logPhyloTimeDivVarRaw
              - log_sum_exp(logPDescBoth + logPhyloTimeDivVarRaw);
        logPhyloDivTimeVarRaw
            = rep_matrix(logMicrobeDivVarRaw, NHostNodes)
              + phyloMetaVar
              + phyloMetaVarP2Scaled;
        logPhyloDivTimeVarRaw
            = logPhyloDivTimeVarRaw
              - log_sum_exp(logPDescBoth + logPhyloDivTimeVarRaw);
        logPhyloDivVarRaw
            = rep_matrix(logHostDivVarRaw, NMicrobeNodes)
              + rep_matrix(logMicrobeDivVarRaw, NHostNodes)
              + phyloMetaVar;
        logPhyloDivVarRaw
            = logPhyloDivVarRaw
              - log_sum_exp(logPDescBoth + logPhyloDivVarRaw);
        logFactVarRaw
            = diag_pre_multiply(metaScales[1:NSubfactors], phyloLogVarMultFacts);
        factMetaVar
            = logFactVarRaw
              * microbeStart;
        logFactTimeVarRaw
            = rep_matrix(logMicrobeTimeVarRaw, NSubfactors)
              + factMetaVar
              + diag_post_multiply(logFactVarRaw,
                                   microbeIntegralScaled);
        logFactDivVarRaw
            = rep_matrix(logMicrobeDivVarRaw, NSubfactors)
              + factMetaVar;
        microbeScales
            = scales[2 * NSubfactors + 3]
              * sqrt(exp(logMicrobeTimeVarRaw + logMicrobeDivVsTime[1])
                     + exp(logMicrobeDivVarRaw + logMicrobeDivVsTime[2]));
        hostScales
            = scales[2 * NSubfactors + 1]
              * sqrt(exp(logHostTimeVarRaw + logHostDivVsTime[1])
                     + exp(logHostDivVarRaw + logHostDivVsTime[2]));
        phyloScales
            = scales[2 * NSubfactors + 2]
              * sqrt(exp(logPhyloTimeVarRaw + logMicrobeDivVsTime[1] + logHostDivVsTime[1])
                     + exp(logPhyloTimeDivVarRaw + logMicrobeDivVsTime[1] + logHostDivVsTime[2])
                     + exp(logPhyloDivTimeVarRaw + logMicrobeDivVsTime[2] + logHostDivVsTime[1])
                     + exp(logPhyloDivVarRaw + logMicrobeDivVsTime[2] + logHostDivVsTime[2]));
        for(f in 1:NSubfactors) {
            factScales[f,]
                = scales[NSubfactors + f]
                  * sqrt(exp(logFactTimeVarRaw[f,]
                             - log_sum_exp(logPDescMicrobe + logFactTimeVarRaw[f,])
                             + logMicrobeDivVsTime[1])
                         + exp(logFactDivVarRaw[f,]
                               - log_sum_exp(logPDescMicrobe + logFactDivVarRaw[f,])
                               + logMicrobeDivVsTime[2]));
        }
        scaledMicrobeNodeEffects[1,1]
            = rawMicrobeNodeEffects[1,1]; //intercept
        scaledMicrobeNodeEffects[1,2:]
            = microbeScales .* (varEffectChol2[1] * rawMicrobeNodeEffects[1,2:]
                                + varEffectCor[1] * phyloLogVarMultPrev); //microbe prevalence
        scaledMicrobeNodeEffects[2:(NEffects + 1),1]
            = subfactLevelMat
              * segment(scales, 1, NSubfactors)
              .* rawMicrobeNodeEffects[2:(NEffects + 1),1]; //samplewise factor alpha diversity
        scaledMicrobeNodeEffects[(NEffects + 2):,1]
            = hostScales .* (varEffectChol2[2] * rawMicrobeNodeEffects[(NEffects + 2):,1]
                             + varEffectCor[2] * phyloLogVarMultADiv); //host alpha diversity
        scaledMicrobeNodeEffects[2:(NEffects + 1),2:]
            = subfactLevelMat
              * factScales
              .* rawMicrobeNodeEffects[2:(NEffects + 1),2:]; //samplewise factor microbe interactions
        scaledMicrobeNodeEffects[(NEffects + 2):,2:]
            = phyloScales .* (varEffectChol2[3] * rawMicrobeNodeEffects[(NEffects + 2):,2:]
                              + varEffectCor[3] * phyloLogVarMultRaw); //host microbe interactions
    }
}
model {
    matrix[NSamples, NMicrobeTips] sampleTipEffects;
    vector[NObs] logit_ratios;
    target += exponential_lpdf(aveStD | 1.0 / aveStDPriorExpect);
    target += dirichSubFact_lpdf;
    target += dirichlet_lpdf(stDProps | rep_vector(1, 2 * NFactors + 3));
    target += exponential_lpdf(aveStDMeta | 1.0 / aveStDMetaPriorExpect);
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
