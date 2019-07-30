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
    matrix[NSamples, NEffects + NHostTips] modelMat;
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
    vector<lower=0, upper=pi()/2>[2 * NFactors + 3] superScales_unif;
    vector<lower=0, upper=pi()/2>[NFactors + 3] superMetaScales_unif;
    vector<lower=0>[2 * NSubfactorGammas] subScalesRaw;
    vector<lower=0>[NSubfactorGammas] subMetaScalesRaw;
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
    vector<lower=0>[2 * NFactors + 3] superScales = tan(superScales_unif);
    vector<lower=0>[NFactors + 3] superMetaScales = tan(superMetaScales_unif);
    vector<lower=0>[2 * NSubfactors + 3] scales;
    vector<lower=0>[NSubfactors + 3] metaScales;
    row_vector[NMicrobeNodes] rawPhyloLogVarMultPrev;
    vector[NHostNodes] rawPhyloLogVarMultADiv;
    matrix[NHostNodes, NMicrobeNodes] rawPhyloLogVarMultRaw;
    matrix[NSubfactors, NMicrobeNodes]  rawPhyloLogVarMultFacts;
    row_vector[NMicrobeNodes] microbeNewEdges;
    vector[NHostNodes] hostNewEdges;
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
                real sum_gamma = 1 + sum(segment(subScalesRaw,
                                                 rawStart,
                                                 NSubPerFactor[i] - 1));
                scales[normStart:(normStart - 1 + NSubPerFactor[i])]
                    = append_row(1, segment(subScalesRaw, rawStart, NSubPerFactor[i] - 1))
                      / sum_gamma;
                dirichSubFact_lpdf += lmultiply(-NSubPerFactor[i], sum_gamma)
                                      + dirichlet_lpdf(scales[normStart:(normStart - 1 + NSubPerFactor[i])] | rep_vector(1, NSubPerFactor[i]));
                scales[normStart:(normStart - 1 + NSubPerFactor[i])]
                    = scales[normStart:(normStart - 1 + NSubPerFactor[i])]
                      * superScales[i];
                sum_gamma = 1 + sum(segment(subScalesRaw,
                                            NSubfactorGammas + rawStart,
                                            NSubPerFactor[i] - 1));
                scales[(NSubfactors + normStart):(NSubfactors + normStart - 1 + NSubPerFactor[i])]
                    = append_row(1, segment(subScalesRaw, NSubfactorGammas + rawStart, NSubPerFactor[i] - 1))
                      / sum_gamma;
                dirichSubFact_lpdf += lmultiply(-NSubPerFactor[i], sum_gamma)
                                      + dirichlet_lpdf(scales[(NSubfactors + normStart):(NSubfactors + normStart - 1 + NSubPerFactor[i])] | rep_vector(1, NSubPerFactor[i]));
                scales[(NSubfactors + normStart):(NSubfactors + normStart - 1 + NSubPerFactor[i])]
                    = scales[(NSubfactors + normStart):(NSubfactors + normStart - 1 + NSubPerFactor[i])]
                      * superScales[NFactors + i];
                sum_gamma = 1 + sum(segment(subMetaScalesRaw,
                                            rawStart,
                                            NSubPerFactor[i] - 1));
                metaScales[normStart:(normStart - 1 + NSubPerFactor[i])]
                    = append_row(1, segment(subMetaScalesRaw, rawStart, NSubPerFactor[i] - 1))
                      / sum_gamma;
                dirichSubFact_lpdf += lmultiply(-NSubPerFactor[i], sum_gamma)
                                      + dirichlet_lpdf(metaScales[normStart:(normStart - 1 + NSubPerFactor[i])] | rep_vector(1, NSubPerFactor[i]));
                metaScales[normStart:(normStart - 1 + NSubPerFactor[i])]
                    = metaScales[normStart:(normStart - 1 + NSubPerFactor[i])]
                      * superMetaScales[i];
                rawStart += NSubPerFactor[i] - 1;
            } else {
                scales[normStart]
                    = superScales[i];
                scales[NSubfactors + normStart]
                    = superScales[NFactors + i];
                metaScales[normStart]
                    = superMetaScales[i];
            }
            normStart += NSubPerFactor[i];
        }
    }
    scales[(2 * NSubfactors + 1):(2 * NSubfactors + 3)]
        = superScales[(2 * NFactors + 1):(2 * NFactors + 3)];
    scales
        = scales
          * aveStDPriorExpect;
    metaScales[(NSubfactors + 1):(NSubfactors + 3)]
        = superScales[(NFactors + 1):(NFactors + 3)];
    metaScales
        = metaScales
          * aveStDMetaPriorExpect;
    {
        row_vector[NMicrobeNodes] logMicrobeVarRaw;
        vector[NHostNodes] logHostVarRaw;
        matrix[NHostNodes, NMicrobeNodes] logPhyloVarRaw;
        matrix[NSubfactors, NMicrobeNodes] logFactVarRaw;
        row_vector[NMicrobeNodes] microbeScales;
        vector[NHostNodes] hostScales;
        matrix[NHostNodes, NMicrobeNodes] phyloScales;
        matrix[NSubfactors, NMicrobeNodes] factScales;
        logMicrobeVarRaw
            = scaledPhyloLogVarMultPrev + logMicrobeEdges;
        microbeNewEdges
            = exp(logMicrobeVarRaw
                  - log_sum_exp(logPDescMicrobe + logMicrobeVarRaw));
        logHostVarRaw
            = scaledPhyloLogVarMultADiv + logHostEdges;
        hostNewEdges
            = exp(logHostVarRaw
                  - log_sum_exp(logPDescHost + logHostVarRaw));
        logPhyloVarRaw
            = scaledPhyloLogVarMultRaw
              + rep_matrix(logHostVarRaw, NMicrobeNodes)
              + rep_matrix(logMicrobeVarRaw, NHostNodes);
        logFactVarRaw
            = scaledPhyloLogVarMultFacts + rep_matrix(logMicrobeVarRaw, NSubfactors);
        microbeScales
            = scales[2 * NSubfactors + 1]
              * sqrt(microbeNewEdges
                     * microbeAncestors');
        hostScales
            = scales[2 * NSubfactors + 2]
              * sqrt(hostAncestors
                     * hostNewEdges);
        phyloScales
            = scales[2 * NSubfactors + 3]
              * sqrt(hostAncestors
                     * exp(logPhyloVarRaw
                           - log_sum_exp(logPDescBoth + logPhyloVarRaw))
                     * microbeAncestors');
        for(f in 1:NSubfactors) {
            factScales[f,]
                = scales[NSubfactors + f]
                  * sqrt(exp(logFactVarRaw[f,]
                             - log_sum_exp(logPDescMicrobe + logFactVarRaw[f,]))
                         * microbeAncestors');
        }
        rawEffectsADiv
            = scaledEffectsADiv ./ (subfactLevelMat * segment(scales, 1, NSubfactors));
        rawHostADiv
            = scaledHostADiv ./ hostScales;
        rawMicrobePrevalence
            = scaledMicrobePrevalence ./ microbeScales;
        rawMicrobeEffectSpecificity
            = (scaledMicrobeEffectSpecificity - rep_matrix(scaledEffectsADiv, NMicrobeNodes)) ./ (subfactLevelMat * factScales);
        rawHostMicrobeSpecificity
            = (scaledHostMicrobeSpecificity - rep_matrix(scaledHostADiv, NMicrobeNodes) - rep_matrix(scaledMicrobePrevalence, NHostNodes) - intercept) ./ phyloScales;
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
    matrix[NSamples, NMicrobeTips] sampleTipEffects;
    vector[NObs] logit_ratios;
    target += dirichSubFact_lpdf;
    target += (- (NSubfactors + 1) * NMicrobeNodes * log(2*pi())
               - trace_quad_form(microbeAncestorsContInv,
                                 append_row(rawPhyloLogVarMultPrev,
                                            rawPhyloLogVarMultFacts)'))
              * 0.5; //matrix normal lpdf for log rate changes in microbe prevalence and factors (http://rpubs.com/mshvarts/matnorm-stan)
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
    sampleTipEffects
        = modelMat * append_row(scaledMicrobeEffectSpecificity[,1:NMicrobeTips], scaledHostMicrobeSpecificity[1:NHostTips, 1:NMicrobeTips]);
    for (n in 1:NObs)
        logit_ratios[n] = sampleTipEffects[sampleNames[n], microbeTipNames[n]];
    target += bernoulli_logit_lpmf(present | logit_ratios);
}
generated quantities {
    vector[2 * NFactors + 3] stDProps = superScales / sum(superScales);
    vector[2 * NSubfactors + 3] subfactProps = scales / sum(superScales);
    vector[NFactors + 3] metaVarProps = superMetaScales / sum(superMetaScales);
    matrix[NSumTo0, NMicrobeNodes + 1] baseLevelEffects
        = baseLevelMat * append_col(scaledEffectsADiv,
                                    scaledMicrobeEffectSpecificity - rep_matrix(scaledEffectsADiv, NMicrobeNodes));
}
