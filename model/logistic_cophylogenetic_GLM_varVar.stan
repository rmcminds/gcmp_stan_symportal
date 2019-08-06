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
    matrix[NMicrobeNodes,NMicrobeNodes] microbeAncestorsInvT
        = inverse(microbeAncestors)';
    matrix[NHostNodes,NHostNodes] hostAncestorsInv
        = inverse(hostAncestors);
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestorsCont
        = microbeAncestors;
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestorsContInvT;
    matrix[NHostNodes, NHostNodes] hostAncestorsCont
        = hostAncestors;
    matrix[NHostNodes, NHostNodes] hostAncestorsContInv;
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
    row_vector[NMicrobeNodes] rawMicrobeVarScales
        = sqrt(microbeEdges);
    vector[NHostNodes] rawHostVarScales
        = sqrt(hostEdges);
    matrix[NHostNodes, NMicrobeNodes] rawHostMicrobeVarScales
        = sqrt(hostEdges
          * microbeEdges);
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
    microbeAncestorsContInvT
        = inverse(microbeAncestorsCont)';
    hostAncestorsContInv
        = inverse(hostAncestorsCont);
}
parameters {
    real<lower=0> globalScaleRaw;
    real<lower=0> globalMetaScaleRaw;
    vector<lower=0>[2 * NFactors + 3] superScalesRaw;
    vector<lower=0>[NFactors + 3] superMetaScalesRaw;
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
    vector<lower=0>[2 * NFactors + 3] superScales
        = aveStDPriorExpect * superScalesRaw;
    vector<lower=0>[NFactors + 3] superMetaScales
        = aveStDMetaPriorExpect * superMetaScalesRaw;
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
                    = sqrt(scales[normStart:(normStart - 1 + NSubPerFactor[i])])
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
                    = sqrt(scales[(NSubfactors + normStart):(NSubfactors + normStart - 1 + NSubPerFactor[i])])
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
                    = sqrt(metaScales[normStart:(normStart - 1 + NSubPerFactor[i])])
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
    metaScales[(NSubfactors + 1):(NSubfactors + 3)]
        = superMetaScales[(NFactors + 1):(NFactors + 3)];
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
              * sqrt(microbeNewEdges);
        hostScales
            = scales[2 * NSubfactors + 2]
              * sqrt(hostNewEdges);
        phyloScales
            = scales[2 * NSubfactors + 3]
              * sqrt(exp(logPhyloVarRaw
                         - log_sum_exp(logPDescBoth + logPhyloVarRaw)));
        for(f in 1:NSubfactors) {
            factScales[f,]
                = scales[NSubfactors + f]
                  * sqrt(exp(logFactVarRaw[f,]
                             - log_sum_exp(logPDescMicrobe + logFactVarRaw[f,])));
        }
        rawEffectsADiv
            = scaledEffectsADiv ./ (subfactLevelMat * segment(scales, 1, NSubfactors));
        rawHostADiv
            = (hostAncestorsInv * scaledHostADiv) ./ hostScales;
        rawMicrobePrevalence
            = (scaledMicrobePrevalence * microbeAncestorsInvT) ./ microbeScales;
        rawMicrobeEffectSpecificity
            = ((scaledMicrobeEffectSpecificity - rep_matrix(scaledEffectsADiv, NMicrobeNodes)) * microbeAncestorsInvT) ./ (subfactLevelMat * factScales);
        rawHostMicrobeSpecificity
            = (hostAncestorsInv * (scaledHostMicrobeSpecificity - rep_matrix(scaledHostADiv, NMicrobeNodes) - rep_matrix(scaledMicrobePrevalence, NHostNodes) - intercept) * microbeAncestorsInvT) ./ phyloScales;
        rawPhyloLogVarMultFacts
            = (scaledPhyloLogVarMultFacts * microbeAncestorsContInvT) ./ rep_matrix(rawMicrobeVarScales, NSubfactors) ./ rep_matrix(metaScales[1:NSubfactors], NMicrobeNodes);
        rawPhyloLogVarMultPrev
            = ((scaledPhyloLogVarMultPrev * microbeAncestorsContInvT) ./ rawMicrobeVarScales / metaScales[NSubfactors + 1]
               - varEffectCor * rawMicrobePrevalence)
              / varEffectChol2;
        rawPhyloLogVarMultADiv
            = ((hostAncestorsContInv * scaledPhyloLogVarMultADiv) ./ rawHostVarScales / metaScales[NSubfactors + 2]
               - varEffectCor * rawHostADiv)
              / varEffectChol2;
        rawPhyloLogVarMultRaw
            = ((hostAncestorsContInv * scaledPhyloLogVarMultRaw * microbeAncestorsContInvT) ./ rawHostMicrobeVarScales / metaScales[NSubfactors + 3]
               - varEffectCor * rawHostMicrobeSpecificity)
              / varEffectChol2;
    }
}
model {
    matrix[NSamples, NMicrobeTips] sampleTipEffects;
    vector[NObs] logit_ratios;
    target += student_t_lpdf(globalScaleRaw | 25,0,1);
    target += student_t_lpdf(globalMetaScaleRaw | 25,0,1);
    target += student_t_lpdf(superScalesRaw | 25,0, globalScaleRaw);
    target += student_t_lpdf(superMetaScalesRaw | 25,0, globalMetaScaleRaw);
    target += dirichSubFact_lpdf;
    target += logistic_lpdf(intercept | 0,1);
    target += std_normal_lpdf(rawEffectsADiv);
    target += std_normal_lpdf(baseLevelMat * rawEffectsADiv);
    target += std_normal_lpdf(rawHostADiv);
    target += std_normal_lpdf(rawMicrobePrevalence);
    target += std_normal_lpdf(to_vector(rawMicrobeEffectSpecificity));
    target += std_normal_lpdf(to_vector(baseLevelMat * rawMicrobeEffectSpecificity));
    target += std_normal_lpdf(to_vector(rawHostMicrobeSpecificity));
    target += std_normal_lpdf(to_vector(rawPhyloLogVarMultFacts));
    target += std_normal_lpdf(rawPhyloLogVarMultPrev);
    target += std_normal_lpdf(rawPhyloLogVarMultADiv);
    target += std_normal_lpdf(to_vector(rawPhyloLogVarMultRaw));
    sampleTipEffects
        = modelMat * append_row(scaledMicrobeEffectSpecificity[,1:NMicrobeTips], scaledHostMicrobeSpecificity[1:NHostTips, 1:NMicrobeTips]);
    for (n in 1:NObs)
        logit_ratios[n] = sampleTipEffects[sampleNames[n], microbeTipNames[n]];
    target += bernoulli_logit_lpmf(present | logit_ratios);
}
generated quantities {
    vector[2 * NFactors + 3] stDProps
        = square(superScales) / sum(square(superScales));
    vector[2 * NSubfactors + 3] subfactProps
        = square(scales) / sum(square(superScales));
    vector[NFactors + 3] metaVarProps
        = square(superMetaScales) / sum(square(superMetaScales));
    matrix[NSumTo0, NMicrobeNodes + 1] baseLevelEffects
        = baseLevelMat * append_col(scaledEffectsADiv,
                                    scaledMicrobeEffectSpecificity - rep_matrix(scaledEffectsADiv, NMicrobeNodes));
}
