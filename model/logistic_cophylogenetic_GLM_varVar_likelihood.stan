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
    int idxMicrobe[NMicrobeNodes];
    int UNidxMicrobe[NMicrobeNodes];
    int idxHost[NHostNodes];
    int UNidxHost[NHostNodes];
    row_vector[NHostNodes] PDescHost
        = rep_row_vector(1.0 / NHostTips, NHostTips)
          * hostTipAncestors;
    vector[NHostNodes] logPDescHost
        = log(PDescHost)';
    vector[NMicrobeNodes] PDescMicrobe
        = microbeTipAncestorsT[2:,]
          * rep_vector(1.0 / NMicrobeTips, NMicrobeTips);
    row_vector[NMicrobeNodes] logPDescMicrobe
        = log(PDescMicrobe)';
    matrix[NHostNodes, NMicrobeNodes] logPDescBoth
        = rep_matrix(logPDescMicrobe, NHostNodes) + rep_matrix(logPDescHost, NMicrobeNodes);
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestorsTCont
        = microbeAncestorsT;
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestorsCont;
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestors;
    matrix[NMicrobeNodes, NMicrobeNodes] microbeAncestorsTContInv;
    matrix[NHostNodes, NHostNodes] hostAncestorsCont
        = hostAncestors;
    matrix[NHostNodes, NHostNodes] hostAncestorsContInv;
    vector[NHostNodes] logHostEdges = log(hostEdges);
    row_vector[NMicrobeNodes] logMicrobeEdges = log(microbeEdges);
    vector[NHostNodes] sqrtHostEdges = sqrt(hostEdges);
    row_vector[NMicrobeNodes] sqrtMicrobeEdges = sqrt(microbeEdges);
    matrix[NHostNodes, NMicrobeNodes] outerEdges
        = sqrt(hostEdges
               * microbeEdges
               * NMicrobeTips
               * NHostTips
               / (NHostTips + NMicrobeTips));
    int NSubfactorGammas = 0;
    int NSubfactors = sum(NSubPerFactor);
    for(i in 1:NMicrobeNodes) {
        microbeAncestorsTCont[i,i] = (1.0 - exp(-microbeEdges[i])) / sqrt(microbeEdges[i]);
    }
    microbeAncestorsCont = microbeAncestorsTCont';
    microbeAncestors = microbeAncestorsT';
    microbeAncestorsTContInv
        = inverse(microbeAncestorsTCont);
    for(i in 1:NHostNodes) {
        hostAncestorsCont[i,i] = (1.0 - exp(-hostEdges[i])) / sqrt(hostEdges[i]);
    }
    hostAncestorsContInv
        = inverse(hostAncestorsCont);
    for(i in 1:NFactors) {
        if(NSubPerFactor[i] > 1) {
            NSubfactorGammas += NSubPerFactor[i] - 1;
        }
    }
    for(i in 1:(NMicrobeNodes - NMicrobeTips)) {
        idxMicrobe[i] = i + NMicrobeTips;
        UNidxMicrobe[i + NMicrobeTips] = i;
    }
    for(i in (NMicrobeNodes - NMicrobeTips + 1):NMicrobeNodes) {
        idxMicrobe[i] = i - NMicrobeNodes + NMicrobeTips;
        UNidxMicrobe[i - NMicrobeNodes + NMicrobeTips] = i;
    }
    for(i in 1:(NHostNodes - NHostTips)) {
        idxHost[i] = i + NHostTips;
        UNidxHost[i + NHostTips] = i;
    }
    for(i in (NHostNodes - NHostTips + 1):NHostNodes) {
        idxHost[i] = i - NHostNodes + NHostTips;
        UNidxHost[i - NHostNodes + NHostTips] = i;
    }
}
parameters {
    real<lower=0> aveStDRaw;
    vector<lower=0>[2 * NSubfactorGammas] subfactPropsRaw;
    simplex[2 * NFactors + 3] stDProps;
    real<lower=0> aveStDMetaRaw;
    simplex[NFactors + 3] metaVarProps;
    vector<lower=0>[NSubfactorGammas] subfactMetaPropsRaw;
    simplex[NMicrobeNodes] microbeVarRaw;
    simplex[NHostNodes] hostVarRaw;
    simplex[NHostNodes * NMicrobeNodes] phyloVarRaw;
    simplex[NMicrobeNodes] factVarRaw[NSubfactors];
    matrix[NEffects, NMicrobeNodes] scaledMicrobeEffectSpecificity;
    matrix[NHostNodes, NMicrobeNodes] scaledHostMicrobeSpecificity;
    row_vector[NMicrobeNodes] scaledMicrobePrevalence;
    vector[NHostNodes] scaledHostADiv;
    vector[NEffects] scaledEffectsADiv;
    real intercept;
    real<lower=-1, upper=1> varEffectCor;
}
transformed parameters {
    simplex[2 * NSubfactors + 3] subfactProps;
    vector[2 * NSubfactors + 3] logScales;
    simplex[NSubfactors + 3] subfactMetaProps;
    vector<lower=0>[NSubfactors + 3] metaScales;
    row_vector<lower=0>[NMicrobeNodes] microbeScales;
    vector<lower=0>[NHostNodes] hostScales;
    matrix[NHostNodes, NMicrobeNodes] phyloScales;
    matrix<lower=0>[NSubfactors, NMicrobeNodes] factScales;
    row_vector[NMicrobeNodes] phyloLogVarMultPrev;
    vector[NHostNodes] phyloLogVarMultADiv;
    vector[NHostNodes] rawHostADiv;
    row_vector[NMicrobeNodes] rawMicrobePrevalence;
    vector[NEffects] rawEffectsADiv;
    matrix[NHostNodes, NMicrobeNodes] phyloLogVarMultRaw;
    matrix[NSubfactors, NMicrobeNodes] phyloLogVarMultFacts;
    matrix[NHostNodes, NMicrobeNodes] rawHostMicrobeSpecificity;
    matrix[NEffects, NMicrobeNodes] rawMicrobeEffectSpecificity;
    real<lower=0> aveStD = aveStDRaw * aveStDPriorExpect;
    real<lower=0> aveStDMeta = aveStDMetaRaw * aveStDMetaPriorExpect;
    real<lower=0, upper=1> varEffectChol2 = sqrt(1 - varEffectCor^2);
    real rawVarJacobians = 0;
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
    logScales
        = log(sqrt((2 * NSubfactors + 3) * subfactProps)
              * aveStD);
    metaScales
        = sqrt((NSubfactors + 3) * subfactMetaProps)
          * aveStDMeta;
    {
        row_vector[NMicrobeNodes] logMicrobeVarRaw = log(microbeVarRaw');
        row_vector[NMicrobeNodes] logMicrobeVar = logMicrobeVarRaw - logPDescMicrobe;
        vector[NHostNodes] logHostVarRaw = log(hostVarRaw);
        vector[NHostNodes] logHostVar = logHostVarRaw - logPDescHost;
        matrix[NHostNodes, NMicrobeNodes] logPhyloVarRaw = to_matrix(log(phyloVarRaw), NHostNodes, NMicrobeNodes);
        matrix[NHostNodes, NMicrobeNodes] logPhyloVar = logPhyloVarRaw - logPDescBoth;
        row_vector[NMicrobeNodes] logFactVarRaw[NSubfactors];
        row_vector[NMicrobeNodes] logFactVar[NSubfactors];
        rawVarJacobians += -sum(logMicrobeVarRaw) - sum(logHostVarRaw) - sum(logPhyloVarRaw);
        for(f in 1:NSubfactors) {
            logFactVarRaw[f] = log(factVarRaw[f]');
            logFactVar[f] = logFactVarRaw[f] - logPDescMicrobe;
            rawVarJacobians += -sum(logFactVarRaw[f]);
        }
        microbeScales
            = exp(logScales[2 * NSubfactors + 3] + logMicrobeVar * 0.5);
        hostScales
            = exp(logScales[2 * NSubfactors + 1] + logHostVar * 0.5);
        phyloScales
            = exp(logScales[2 * NSubfactors + 2] + logPhyloVar * 0.5);
        for(f in 1:NSubfactors) {
            factScales[f,]
                = exp(logScales[NSubfactors + f] + logFactVar[f,] * 0.5);
        }
        rawMicrobePrevalence
            = scaledMicrobePrevalence ./ microbeScales;
        rawEffectsADiv
            = scaledEffectsADiv ./ (subfactLevelMat * segment(exp(logScales), 1, NSubfactors));
        rawHostADiv
            = scaledHostADiv ./ hostScales;
        rawMicrobeEffectSpecificity
            = scaledMicrobeEffectSpecificity ./ (subfactLevelMat * factScales);
        rawHostMicrobeSpecificity
            = scaledHostMicrobeSpecificity ./ phyloScales;
        phyloLogVarMultPrev
            = (mdivide_left_tri_low(microbeAncestorsCont[idxMicrobe,idxMicrobe],
                                    (logMicrobeVar - logMicrobeEdges)[idxMicrobe]')'[UNidxMicrobe]
               / metaScales[NSubfactors + 1]
               ./ sqrtMicrobeEdges
               - (varEffectCor * rawMicrobePrevalence))
              / varEffectChol2;
        phyloLogVarMultADiv
            = (mdivide_left_tri_low(hostAncestorsCont[idxHost,idxHost],
                                    (logHostVar - logHostEdges)[idxHost])[UNidxHost]
               / metaScales[NSubfactors + 2]
               ./ sqrtHostEdges
               - (varEffectCor * rawHostADiv))
              / varEffectChol2;
        phyloLogVarMultRaw
            = (mdivide_left_tri_low(hostAncestorsCont[idxHost,idxHost],
                                    mdivide_left_tri_low(microbeAncestorsCont[idxMicrobe,idxMicrobe],
                                                         (logPhyloVarRaw
                                                          - rep_matrix(logHostVar, NMicrobeNodes)
                                                          - rep_matrix(logMicrobeVar, NHostNodes))[idxHost,idxMicrobe]')')[UNidxHost,UNidxMicrobe]
               / metaScales[NSubfactors + 3]
               ./ outerEdges
               - (varEffectCor * rawHostMicrobeSpecificity))
              / varEffectChol2;
        for(f in 1:NSubfactors) {
            phyloLogVarMultFacts[f,]
                = mdivide_left_tri_low(microbeAncestorsCont[idxMicrobe,idxMicrobe],
                                       (logFactVar[f] - logMicrobeVar)[idxMicrobe]')'[UNidxMicrobe]
                  ./ (metaScales[f]
                      * sqrtMicrobeEdges);
        }
    }
}
model {
    matrix[NEffects + NHostNodes + 1, NMicrobeTips] scaledMicrobeNodeEffects;
    matrix[NSamples, NMicrobeTips] sampleTipEffects;
    vector[NObs] logit_ratios;
/*
    target += rawVarJacobians;
    target += exponential_lpdf(aveStDRaw | 1.0);
    target += dirichSubFact_lpdf;
    target += dirichlet_lpdf(stDProps | rep_vector(1, 2 * NFactors + 3));
    target += exponential_lpdf(aveStDMetaRaw | 1.0);
    target += dirichlet_lpdf(metaVarProps | rep_vector(1, NFactors + 3));
    target += std_normal_lpdf(phyloLogVarMultPrev);
    target += std_normal_lpdf(phyloLogVarMultADiv);
    target += std_normal_lpdf(to_vector(phyloLogVarMultRaw));
    target += std_normal_lpdf(to_vector(phyloLogVarMultFacts));
    for(e in 1:NEffects) {
        target += multi_normal_cholesky_lpdf(rawMicrobeEffectSpecificity[e,idxMicrobe] |
                                             rep_vector(0,NMicrobeNodes),
                                             microbeAncestors[idxMicrobe,idxMicrobe]);
    }
    for(s in 1:NSumTo0) {
        target += multi_normal_cholesky_lpdf(baseLevelMat[s,] * rawMicrobeEffectSpecificity[,idxMicrobe] |
                                             rep_vector(0,NMicrobeNodes),
                                             microbeAncestors[idxMicrobe,idxMicrobe]);
    }
    for(h in 1:NHostNodes) {
        target += multi_normal_cholesky_lpdf(rawHostMicrobeSpecificity[h,idxMicrobe] |
                                             rep_vector(0,NMicrobeNodes),
                                             microbeAncestors[idxMicrobe,idxMicrobe]);
    }
    target += multi_normal_cholesky_lpdf(rawMicrobePrevalence[idxMicrobe] |
                                         rep_vector(0,NMicrobeNodes),
                                         microbeAncestors[idxMicrobe,idxMicrobe]);
    target += std_normal_lpdf(rawHostADiv);
    target += std_normal_lpdf(rawEffectsADiv);
    target += logistic_lpdf(intercept | 0,1);
    target += uniform_lpdf(varEffectCor | -1,1);
*/
    scaledMicrobeNodeEffects[1,]
        = scaledMicrobePrevalence[1:NMicrobeTips]; //microbe prevalence
    scaledMicrobeNodeEffects[2:(NEffects + 1),]
        = scaledMicrobeEffectSpecificity[,1:NMicrobeTips]; //samplewise factor microbe interactions
    scaledMicrobeNodeEffects[(NEffects + 2):,]
        = scaledHostMicrobeSpecificity[,1:NMicrobeTips]; //host microbe interactions
    sampleTipEffects = rep_matrix(modelMat
                                  * append_row(intercept,
                                               append_row(scaledEffectsADiv,
                                                          scaledHostADiv)),
                                  NMicrobeTips)
                       + modelMat * scaledMicrobeNodeEffects;
    for (n in 1:NObs)
        logit_ratios[n] = sampleTipEffects[sampleNames[n], microbeTipNames[n]];
    target += bernoulli_logit_lpmf(present | logit_ratios);
}
generated quantities {
    row_vector[NMicrobeNodes] microbeNewEdges
        = square(microbeScales / exp(logScales[2 * NSubfactors + 3]));
    vector[NHostNodes] hostNewEdges
        = square(hostScales / exp(logScales[2 * NSubfactors + 1]));
    matrix[NSumTo0, NMicrobeNodes + 1] baseLevelEffects
        = baseLevelMat * append_col(scaledEffectsADiv,
                                    scaledMicrobeEffectSpecificity);
}
