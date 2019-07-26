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
    matrix[NSamples, NEffects + NHostTips + 1] modelMat;
    int NSumTo0;
    matrix[NSumTo0, NEffects] baseLevelMat;
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestors;
    matrix[NHostNodes, NHostNodes] hostAncestors;
    vector[NHostNodes] hostEdges;
    row_vector[NMicrobeNodes] microbeEdges;
}
transformed data {
    cov_matrix[NMicrobeNodes] microbeAncestorsLLInv
        = inverse_spd(tcrossprod(microbeAncestors));
    cov_matrix[NHostNodes] hostAncestorsLLInv
        = inverse_spd(tcrossprod(hostAncestors));
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestorsCont
        = microbeAncestors;
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestorsContInv;
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestorsContXInvT;
    matrix[NHostNodes, NHostNodes] hostAncestorsCont
        = hostAncestors;
    matrix[NHostNodes, NHostNodes] hostAncestorsContInv;
    matrix[NHostNodes, NHostNodes] hostAncestorsContXInv;
    row_vector[NHostNodes] PDescHost
        = rep_row_vector(1.0 / NHostTips, NHostTips)
          * hostAncestors[1:NHostTips,];
    vector[NHostNodes] logPDescHost
        = log(PDescHost)';
    row_vector[NMicrobeNodes] PDescMicrobe
        = rep_row_vector(1.0 / NMicrobeTips, NMicrobeTips)
          * microbeAncestors[1:NMicrobeTips,];
    row_vector[NMicrobeNodes] logPDescMicrobe
        = log(PDescMicrobe);
    matrix[NHostNodes, NMicrobeNodes] logPDescBoth
        = rep_matrix(logPDescMicrobe, NHostNodes) + rep_matrix(logPDescHost, NMicrobeNodes);
    vector[NHostNodes] logHostEdges
        = log(hostEdges);
    row_vector[NMicrobeNodes] logMicrobeEdges
        = log(microbeEdges);
    vector[NHostNodes] sqrtHostEdges
        = sqrt(hostEdges);
    row_vector[NMicrobeNodes] sqrtMicrobeEdges
        = sqrt(microbeEdges);
    row_vector[NMicrobeNodes] rawMicrobeVarScales;
    vector[NHostNodes] rawHostVarScales;
    matrix[NHostNodes, NMicrobeNodes] outerEdges;
    int NSubfactorGammas
        = 0;
    int NSubfactors
        = sum(NSubPerFactor);
    for(i in 1:NMicrobeNodes) {
        microbeAncestorsCont[i,i] = (1.0 - exp(-1));
    }
    for(i in 1:NHostNodes) {
        hostAncestorsCont[i,i] = (1.0 - exp(-1));
    }
    for(i in 1:NFactors) {
        if(NSubPerFactor[i] > 1) {
            NSubfactorGammas += NSubPerFactor[i] - 1;
        }
    }
    rawMicrobeVarScales
        = sqrtMicrobeEdges * microbeAncestorsCont';
    rawHostVarScales
        = hostAncestorsCont * sqrtHostEdges;
    outerEdges
        = hostAncestorsCont
          * sqrt(hostEdges
                 * microbeEdges
                 * NMicrobeTips
                 * NHostTips
                 / (NHostTips + NMicrobeTips))
          * microbeAncestorsCont';
    microbeAncestorsContInv
        = inverse_spd(tcrossprod(microbeAncestorsCont));
    hostAncestorsContInv
        = inverse_spd(tcrossprod(hostAncestorsCont));
    microbeAncestorsContXInvT
        = (microbeAncestorsCont / microbeAncestors)';
    hostAncestorsContXInv
        = hostAncestorsCont / hostAncestors;
}
parameters {
    real<lower=0> aveStDRaw;
    vector<lower=0>[2 * NSubfactorGammas] subfactPropsRaw;
    simplex[2 * NFactors + 3] stDProps;
    real<lower=0> aveStDMetaRaw;
    simplex[NFactors + 3] metaVarProps;
    vector<lower=0>[NSubfactorGammas] subfactMetaPropsRaw;
    row_vector[NMicrobeNodes] scaledPhyloLogVarMultPrev;
    vector[NHostNodes] scaledPhyloLogVarMultADiv;
    matrix[NHostNodes, NMicrobeNodes] scaledPhyloLogVarMultRaw;
    matrix[NSubfactors, NMicrobeNodes]  scaledPhyloLogVarMultFacts;
    real intercept;
    matrix[NHostNodes, NMicrobeNodes] scaledHostMicrobeSpecificity;
    matrix[NEffects, NMicrobeNodes] scaledMicrobeEffectSpecificity;
    vector[NHostNodes] scaledHostADiv;
    row_vector[NMicrobeNodes] scaledMicrobePrevalence;
    vector[NEffects] scaledEffectsADiv;
    real<lower=-1, upper=1> varEffectCor;
}
transformed parameters {
    real<lower=0> aveStD = aveStDRaw * aveStDPriorExpect;
    real<lower=0> aveStDMeta = aveStDMetaRaw * aveStDMetaPriorExpect;
    simplex[2 * NSubfactors + 3] subfactProps;
    vector<lower=0>[2 * NSubfactors + 3] scales;
    simplex[NSubfactors + 3] subfactMetaProps;
    vector<lower=0>[NSubfactors + 3] metaScales;
    row_vector[NMicrobeNodes] rawPhyloLogVarMultPrev;
    vector[NHostNodes] rawPhyloLogVarMultADiv;
    matrix[NHostNodes, NMicrobeNodes] rawPhyloLogVarMultRaw;
    matrix[NSubfactors, NMicrobeNodes]  rawPhyloLogVarMultFacts;
    row_vector<lower=0>[NMicrobeNodes] microbeScales;
    vector<lower=0>[NHostNodes] hostScales;
    matrix<lower=0>[NSubfactors, NMicrobeNodes] factScales;
    matrix[NHostNodes, NMicrobeNodes] rawHostMicrobeSpecificity;
    matrix[NEffects, NMicrobeNodes] rawMicrobeEffectSpecificity;
    vector[NHostNodes] rawHostADiv;
    row_vector[NMicrobeNodes] rawMicrobePrevalence;
    vector[NEffects] rawEffectsADiv;
    real<lower=0, upper=1> varEffectChol2 = sqrt(1 - varEffectCor^2);
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
        = sqrt((2 * NFactors + 3) * subfactProps)
          * aveStD;
    metaScales
        = sqrt((NFactors + 3) * subfactMetaProps)
          * aveStDMeta;
    {
        row_vector[NMicrobeNodes] logMicrobeVarRaw;
        vector[NHostNodes] logHostVarRaw;
        matrix[NHostNodes, NMicrobeNodes] logPhyloVarRaw;
        matrix[NSubfactors, NMicrobeNodes] logFactVarRaw;
        matrix[NHostNodes, NMicrobeNodes] phyloScales;
        logMicrobeVarRaw
            = scaledPhyloLogVarMultPrev + logMicrobeEdges;
        logHostVarRaw
            = scaledPhyloLogVarMultADiv + logHostEdges;
        logPhyloVarRaw
            = scaledPhyloLogVarMultRaw
              + rep_matrix(logHostVarRaw, NMicrobeNodes)
              + rep_matrix(logMicrobeVarRaw, NHostNodes);
        logFactVarRaw
            = scaledPhyloLogVarMultFacts + rep_matrix(logMicrobeVarRaw, NSubfactors);
        microbeScales
            = scales[2 * NSubfactors + 3]
              * exp((logMicrobeVarRaw
                     - log_sum_exp(logPDescMicrobe + logMicrobeVarRaw))
                    * 0.5)
              * microbeAncestors';
        hostScales
            = scales[2 * NSubfactors + 1]
              * hostAncestors
              * exp((logHostVarRaw
                     - log_sum_exp(logPDescHost + logHostVarRaw))
                    * 0.5);
        phyloScales
            = scales[2 * NSubfactors + 2]
              * hostAncestors
              * exp((logPhyloVarRaw
                     - log_sum_exp(logPDescBoth + logPhyloVarRaw))
                    * 0.5)
              * microbeAncestors';
        for(f in 1:NSubfactors) {
            factScales[f,]
                = scales[NSubfactors + f]
                  * exp((logFactVarRaw[f,]
                        - log_sum_exp(logPDescMicrobe + logFactVarRaw[f,]))
                        * 0.5)
                  * microbeAncestors';
        }
        rawEffectsADiv
            = scaledEffectsADiv ./ (subfactLevelMat * segment(scales, 1, NSubfactors));
        rawMicrobePrevalence
            = scaledMicrobePrevalence ./ microbeScales;
        rawHostADiv
            = scaledHostADiv ./ hostScales;
        rawMicrobeEffectSpecificity
            = scaledMicrobeEffectSpecificity ./ (subfactLevelMat * factScales);
        rawHostMicrobeSpecificity
            = scaledHostMicrobeSpecificity ./ phyloScales;
        rawPhyloLogVarMultFacts
            = scaledPhyloLogVarMultFacts ./ rep_matrix(rawMicrobeVarScales, NSubfactors) ./ rep_matrix(metaScales[1:NSubfactors], NMicrobeNodes);
        rawPhyloLogVarMultPrev
            = (scaledPhyloLogVarMultPrev ./ rawMicrobeVarScales / metaScales[NSubfactors + 1]
               - varEffectCor * rawMicrobePrevalence * microbeAncestorsContXInvT)
              / varEffectChol2;
        rawPhyloLogVarMultADiv
            = (scaledPhyloLogVarMultADiv ./ rawHostVarScales / metaScales[NSubfactors + 2]
               - varEffectCor * hostAncestorsContXInv * rawHostADiv)
              / varEffectChol2;
        rawPhyloLogVarMultRaw
            = (scaledPhyloLogVarMultRaw ./ outerEdges / metaScales[NSubfactors + 3]
               - varEffectCor * hostAncestorsContXInv * rawHostMicrobeSpecificity * microbeAncestorsContXInvT)
              / varEffectChol2;
    }
}
model {
    matrix[NEffects + NHostTips + 1, NMicrobeTips] scaledMicrobeNodeEffects;
    matrix[NSamples, NMicrobeTips] sampleTipEffects;
    vector[NObs] logit_ratios;
    target += exponential_lpdf(aveStDRaw | 1.0);
    target += dirichSubFact_lpdf;
    target += dirichlet_lpdf(stDProps | rep_vector(1, 2 * NFactors + 3));
    target += exponential_lpdf(aveStDMetaRaw | 1.0);
    target += dirichlet_lpdf(metaVarProps | rep_vector(1, NFactors + 3));
    target += (- (NSubfactors + 1) * NMicrobeNodes * log(2*pi())
               - trace_quad_form(microbeAncestorsContInv,
                                 append_row(rawPhyloLogVarMultPrev,
                                            rawPhyloLogVarMultFacts)'))
              * 0.5; //matrix normal lpdf for log rate changes in microbe prevalence and factors
    target += multi_normal_cholesky_lpdf(rawPhyloLogVarMultADiv |
                                         rep_vector(0,NHostNodes),
                                         hostAncestorsCont);
    target += (- NHostNodes * NMicrobeNodes * log(2*pi())
               - trace_gen_quad_form(microbeAncestorsContInv,
                                     hostAncestorsContInv,
                                     rawPhyloLogVarMultRaw))
              * 0.5; //matrix normal lpdf for log rate changes in host-microbe interactions
    target += logistic_lpdf(intercept | 0,1);
    target += std_normal_lpdf(rawEffectsADiv);
    target += std_normal_lpdf(baseLevelMat * rawEffectsADiv);
    target += (- (NEffects + NSumTo0 + 1) * NMicrobeNodes * log(2*pi())
               - trace_quad_form(microbeAncestorsLLInv,
                                 append_row(rawMicrobePrevalence,
                                   append_row(rawMicrobeEffectSpecificity,
                                              baseLevelMat * rawMicrobeEffectSpecificity))'))
              * 0.5; //matrix normal lpdf for microbe prevalence and effect-microbe interactions
    target += multi_normal_cholesky_lpdf(rawHostADiv |
                                         rep_vector(0,NHostNodes),
                                         hostAncestors);
    target += (- NHostNodes * NMicrobeNodes * log(2*pi())
               - trace_gen_quad_form(microbeAncestorsLLInv,
                                     hostAncestorsLLInv,
                                     rawHostMicrobeSpecificity))
              * 0.5; //matrix normal lpdf for host-microbe interactions
    scaledMicrobeNodeEffects
        = rep_matrix(append_row(intercept,
                                append_row(scaledEffectsADiv,
                                           scaledHostADiv[1:NHostTips])),
                     NMicrobeTips); // intercept and adiv effects
    scaledMicrobeNodeEffects[1,]
        += scaledMicrobePrevalence[1:NMicrobeTips]; // microbe prevalence
    scaledMicrobeNodeEffects[2:(NEffects + 1),]
        += scaledMicrobeEffectSpecificity[,1:NMicrobeTips]; // samplewise factor microbe interactions
    scaledMicrobeNodeEffects[(NEffects + 2):,]
        += scaledHostMicrobeSpecificity[1:NHostTips, 1:NMicrobeTips]; // host microbe interactions
    sampleTipEffects = modelMat * scaledMicrobeNodeEffects;
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
        = baseLevelMat * append_col(scaledEffectsADiv,
                                    scaledMicrobeEffectSpecificity);
}
