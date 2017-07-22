library(dplyr)
library(tidyr)
library(ggplot2)

# Read-in belowground
{
  # Soil C, N data
  soilCN <- read.table("/home/pradeep/Workspace/MSc/Project/soc-fragments/Pradeep_CoorgSoilCN.csv",
                       as.is=T, header=T, sep=",")
  
  # Clean-up col names, drop needless cols
  soilCN <- dplyr::select(soilCN, c(SAMPLE, MASS, C_PERC, N_PERC))
  soilCN <- soilCN %>% rename(POINT=SAMPLE) %>% rename(SOIL.C=C_PERC) %>%
    rename(SOIL.N=N_PERC) %>% rename(SOIL.WT=MASS)
  
  # Remove "Blank", "Soil 502-309" standard sample rows
  soilCN <- filter(soilCN, grepl("CT|FR", POINT))
  
  # Copy CT4_PT0 row to CT5_PT0 row since 
  # they are the same point
  t <- soilCN[grep("CT4_PT0",soilCN$POINT), ]
  t$POINT[1] <- "CT5_PT0"
  soilCN <- bind_rows(soilCN, t)
  
  # Calculate "adjusted" soil C by subtracting out value of "dummy" point
  # collected at bare patch inside/outside the site
  soilCN <- separate(soilCN, col="POINT",
                     into=c("SITE.ID", "POINT.TYPE", "DUMMY", "POINT.IN.SITE", "POINT.NO"),
                     sep=c(3,4,6,7), remove=F)
  soilCN <- dplyr::select(soilCN, -DUMMY)
  soilCN$SITE.ID <- as.factor(soilCN$SITE.ID)
  soilCN <- soilCN %>% group_by(SITE.ID) %>%
    mutate(ADJ.SOIL.C=SOIL.C-SOIL.C[POINT.IN.SITE=="0"])
  
  # Group points according to sites
  soilCN <- arrange(soilCN, POINT)
  
  # Soil SIR data
  soilSIR <- read.table("/home/pradeep/Workspace/MSc/Project/soc-fragments/soilSIR.csv",
                        as.is=T, header=T, sep=",")
  
  concNaOH <- 2 # in N (normality)
  
  # Titre value corresponding to substrate-induced microbial respiration alone
  # Calculated as (total sample titre - "blank" titre)
  avgBlanks <- soilSIR %>% group_by(SET) %>% filter(grepl("B", SAMPLE)) %>% 
    summarise(AVG.BLANKS=mean(TITRE))
  soilSIRS1 <- soilSIR %>% filter(SET == "S1") %>% 
    mutate(ADJ.TITRE = -TITRE+avgBlanks$AVG.BLANKS[avgBlanks$SET=="S1"])
  soilSIRS2 <- soilSIR %>% filter(SET == "S2") %>% 
    mutate(ADJ.TITRE = -TITRE+avgBlanks$AVG.BLANKS[avgBlanks$SET=="S2"])
  soilSIR <- bind_rows(soilSIRS1, soilSIRS2)
  
  # Ignore treated samples that respired respired less than the "blank"s
  soilSIR$ADJ.TITRE[soilSIR$ADJ.TITRE<0] = NA  
  soilSIR <- filter(soilSIR, !grepl("B", SAMPLE))
  
  # SIR <- Compute C in CO2 released by microbes due to substrate treatment.
  # Computed in units:  micro gram C-CO2 per hour, per gram dry weight soil
  #
  # --             --   --                                      --
  #| Tadj_ml * HCl_N | |                    12_g                  |
  #| --------------- |*| ---------------------------------------- |*10^6_mugPerg
  #|     NaOH_N      | |  1000_mlPerLitre * 26_hr * 10_gDryWtSoil |
  # --             --   --                                      --
  #
  soilSIR <- mutate(soilSIR, C.RELEASED = 1000000*
                      (ADJ.TITRE*(HCL_N/concNaOH))*12/(1000*26*10))
  # Ignoring outliers
  soilSIR$C.RELEASED[soilSIR$C.RELEASED>25] = NA
  soilSIR$C.RELEASED[soilSIR$C.RELEASED<0.18] = NA
  
  soilSIR <- soilSIR %>% rename(POINT=SAMPLE) %>%
    dplyr::select(POINT, C.RELEASED)
  
  # Merging soil C, N, SIR into one dataframe
  belowGround <- left_join(soilCN, soilSIR, by="POINT")
  belowGround <- rename(belowGround, SOIL.SIR=C.RELEASED)
}

# Read-in aboveground data
{
  # Litter C, N data
  litterCN <- read.table("/home/pradeep/Workspace/MSc/Project/soc-fragments/Pradeep_CoorgLitterCN.csv",
                         as.is=T, header=T, sep=",")
  litterCN <- dplyr::select(litterCN, c(SAMPLE, C_PERC, N_PERC))
  litterCN <- litterCN %>% rename(POINT=SAMPLE) %>% rename(LITTER.C=C_PERC) %>%
    rename(LITTER.N=N_PERC)
  
  # Remove "Blank", "EDTA*" standard sample rows
  litterCN <- filter(litterCN, grepl("CT|FR", POINT))
  
  # Group points according to sites
  litterCN <- arrange(litterCN, POINT)
  
  # Calculate C/N ratio
  litterCN <- mutate(litterCN, LITTER.CN.RATIO=LITTER.C/LITTER.N)
  
  # Litter P data
  litterP <- read.table("/home/pradeep/Workspace/MSc/Project/soc-fragments/ICP_calculated_P_pradeep.csv",
                        as.is=T, header=T, sep=",")
  litterP <- dplyr::select(litterP, c(SAMPLE, P_PERC))
  litterP <- litterP %>% rename(POINT=SAMPLE) %>% rename(LITTER.P=P_PERC) %>%
    filter(grepl("CT|FR", POINT))
  
  # Litter weights
  litterWts <- read.table("/home/pradeep/Workspace/MSc/Project/soc-fragments/litterProperties.csv",
                          as.is=T, header=T, sep=",")
  litterWts$AVGWT <- rowMeans(litterWts[ ,c("DRYWT1", "DRYWT2", "DRYWT3")],
                              na.rm=TRUE)
  litterAvgWt <- litterWts %>% dplyr::select(c(POINT, AVGWT)) %>%
    rename(LITTER.WT = AVGWT)
  # Standardise weights to per m^2
  #              litter wt from 0.3m * 0.3m quadrat
  # Litter wt. = ----------------------------------  per m^2
  #                           (0.3)^2
  litterAvgWt$LITTER.WT = litterAvgWt$LITTER.WT/(0.3*0.3)
  
  # Consolidate litter C, N, P, weight
  litter <- left_join(litterCN, litterP, by="POINT")
  litter <- left_join(litter, litterAvgWt, by="POINT")
  
  # Tree ommunity structure data
  treeDistrib <- read.table("/home/pradeep/Workspace/MSc/Project/soc-fragments/treeInfluencePotentialData.csv",
                            as.is=T, header=T, sep=",")
  
  # Effective DBH : sqrt(dbh1^2 + dbh2^2 + dbh3^2 + dbh4^2)
  gbhs <- dplyr::select(treeDistrib, c(GBH1, GBH2, GBH3, GBH4))
  treeDistrib$DBH <- (gbhs/pi)^2 %>% rowSums(na.rm=T) %>% sqrt()
  
  treeDistrib <- treeDistrib %>% separate(col="POINT",
                                          into=c("SITE.TYPE", "T"), sep=2,
                                          remove=F) %>%
    dplyr::select(-T) %>%
    separate(col="POINT",
             into=c("SITE.ID", "PT.TYPE", "T", "PT.IN.SITE", "TT"),
             sep=c(3,4,6,7), remove=F) %>%
    dplyr::select(-T, -TT)
  
  # Ignore trees with, DBH < 10 cm (0.1 m),
  #             point type = "D" (from "degraded" parts of sites (fragments))
  #          point in site = "0" (from "blank" points (outside sites))
  treeDistrib <- treeDistrib %>% filter(DBH>=0.1) %>% filter(PT.TYPE!="D") %>%
    filter(PT.IN.SITE==1) %>%
    dplyr::select(POINT, SITE.TYPE, SITE.ID, DISTANCE, DBH)
  
  # Compute basal area
  treeDistrib$TreeBasArea <- (pi*(treeDistrib$DBH)^2)/4
  # Estimate of site-level basal area (m^2 per ha), calculated at each point
  # (sum of basal area per tree)_m^2
  # -------------------------------- * 10,000_m^2PerHa
  #           pi*10_m^2 
  basArea <- treeDistrib %>% group_by(POINT) %>% 
    summarise(BAS.AREA=sum(TreeBasArea, na.rm=T)*10000/(pi*10^2))
  # Estimate of site-level tree density (per ha), calculated at each point
  #         no. of trees
  #         ------------ * 10,000_m^2PerHa
  #          pi*10_m^2 
  treeDens <- treeDistrib %>% group_by(POINT) %>% 
    summarise(TREE.DENS=n()*10000/(pi*10^2))
  
  # Proportion of large trees ... 
  binBounds = c(0.1, 0.15, 0.25, 0.45, 0.85, 2)
  treeDistrib$SIZE.CLASS[(binBounds[1]<treeDistrib$DBH) &
                           (treeDistrib$DBH<=binBounds[2])] = "S1"
  treeDistrib$SIZE.CLASS[(binBounds[2]<treeDistrib$DBH) &
                           (treeDistrib$DBH<=binBounds[3])] = "S2"
  treeDistrib$SIZE.CLASS[(binBounds[3]<treeDistrib$DBH) &
                           (treeDistrib$DBH<=binBounds[4])] = "S3"
  treeDistrib$SIZE.CLASS[(binBounds[4]<treeDistrib$DBH) &
                           (treeDistrib$DBH<=binBounds[5])] = "S4"
  treeDistrib$SIZE.CLASS[(binBounds[5]<treeDistrib$DBH) &
                           (treeDistrib$DBH<=binBounds[6])] = "S5"
  treeDistrib$SIZE.CLASS <- as.ordered((treeDistrib$SIZE.CLASS))
  
  largeTreeDens <- treeDistrib %>% mutate(IS.LARGE = DBH>0.7) %>%
    group_by(POINT) %>%
    summarise(NO.LARGE.TREES = sum(IS.LARGE)*10000/(pi*10^2), 
              PROP.LARGE.TREES = sum(IS.LARGE)/n())
  largeBasalAreaFrac <- treeDistrib %>% mutate(IS.LARGE = DBH>0.7) %>%
    group_by(POINT) %>%
    summarise(LARGE.TREE.BA=
                sum(TreeBasArea[DBH>0.7],na.rm=T)*10000/(pi*10^2),
              LARGE.TREE.BA.FRAC=
                sum(TreeBasArea[DBH>0.7])/sum(TreeBasArea))
  largeTreeDens <- left_join(largeTreeDens, largeBasalAreaFrac, by="POINT")
  
  # Compute tree density by size class
  treeDensBySizeClass <- treeDistrib %>% group_by(POINT, SIZE.CLASS) %>% 
    summarise(TREE.DENS=n()*10000/(pi*10^2))
  treeDensBySizeClass <- treeDensBySizeClass %>% group_by(POINT) %>% 
    mutate(PROP.STEMS=TREE.DENS/sum(TREE.DENS))
  treeDensBySizeClass <- treeDensBySizeClass %>% 
    separate(col="POINT", into=c("SITE.TYPE", "T"), sep=2,
             remove=F) %>%
    dplyr::select(-T) %>% 
    separate(col="POINT", into=c("SITE.ID", "T"), sep=3,
             remove=F) %>% 
    dplyr::select(-T)
  treeDensMeans <- treeDensBySizeClass %>% group_by(SITE.TYPE, SIZE.CLASS) %>%
    summarise(MEAN=mean(TREE.DENS),
              SD=sd(TREE.DENS),
              SE=sd(TREE.DENS)/sqrt(n()))
  treeStemPropsMeans <- treeDensBySizeClass %>%
    group_by(SITE.TYPE, SIZE.CLASS) %>%
    summarise(MEAN=mean(PROP.STEMS),
              SD=sd(PROP.STEMS),
              SE=sd(PROP.STEMS)/sqrt(n()))
  
  aboveGround <- left_join(litter, basArea, by="POINT")
  aboveGround <- left_join(aboveGround, treeDens, by="POINT")
  
}

# Anand's species-wise leaf litter quality data
{
  # Replicated species-wise litter quality data
  anandSpWiseLitter <- read.table("/home/pradeep/Workspace/MSc/Project/scratchPad/Anand_specieswise_litter_CN.csv",
                                  as.is=T, header=T, sep=",")
  anandSpWiseLitter <- anandSpWiseLitter %>% mutate (SP.CTON = C../N.) 
  # Average litter quality for each species
  speciesWiseLitterQual <- anandSpWiseLitter %>% 
    group_by(species) %>% summarise(CN.RATIO=mean(SP.CTON))
  
  # Plot-level tree-id and basal area measurements. 
  anandSpDominance <- read.table("/home/pradeep/Workspace/MSc/Project/scratchPad/Anand_TreeCommunity_SpeciesDominance_plotData.csv",
                                 as.is=T, header=T, sep=",")
  anandSpDominance <- unite(anandSpDominance, SITE.PLOT.ID, site.name, plot.id,
                            remove=F)
  # Fractional contribution of each tree to its plot-level basal area
  anandSpDominance <- anandSpDominance %>% group_by(SITE.PLOT.ID) %>% 
    arrange(desc(basal.area_sq.cm)) %>% 
    mutate(BAS.AREA.FRAC = basal.area_sq.cm/sum(basal.area_sq.cm)) %>%
    mutate(BAS.AREA.CUMFRAC = cumsum(BAS.AREA.FRAC))
  # Biggest trees that contribute to top 70% basal area at plot level
  anand70PercDom <- anandSpDominance %>%
    dplyr::select(-index, -date, -family, -diameter_cm, -height_m, -Remarks) %>%
    filter(BAS.AREA.CUMFRAC <= 0.7)
  
  # LUT for tree id and code used
  anandSpeciesCodes <- read.table("/home/pradeep/Workspace/MSc/Project/scratchPad/Anand_TreeCommunity_Families.csv",
                                  as.is=T, header=T, sep=",")
  anandSpeciesCodes <- anandSpeciesCodes %>% dplyr::select(code, combine) %>%
    rename(species=combine)
  
  # Retain dominant species for which leaf litter qual data is available
  anand70PercDom <- inner_join(anandSpeciesCodes, anand70PercDom, by="species")
  anand70PercDom <- right_join(speciesWiseLitterQual, anand70PercDom,
                               by= c("species"="code"))
  
  # Select sites where I have sampled
  anand70PercDomMySites <- filter (anand70PercDom, 
                                   site.name=="Kokka" |
                                     site.name=="Ruduraguppae" |
                                     site.name=="Arji" |
                                     site.name=="Arpattu")
  
  # Compute community-averaged C/N.
  #                w1*cn1 + w2*cn2 + ... + wN*cnN
  # avg litt C/N = ------------------------------
  #                      w1 + w2 + ... wN
  # wN : fractional basal area of Nth tree
  # cnN: leaf litter C/N of Nth tree
  # (w1 + w2 + ... wN = 0.7, since only those trees included that contribute to
  # top 70% of basal area)
  commAvgdLittQualMySites <- anand70PercDomMySites %>%
    group_by(SITE.PLOT.ID) %>%
    summarise(SITE.AVG.CN.RATIO =
                sum(CN.RATIO*BAS.AREA.FRAC,na.rm=T)/0.7,
              SITE.TYPE = first(type))
  noLargeTreesMySitesAnandData <- anand70PercDomMySites %>% group_by(type) %>%
    summarise(NO.LRG.TREES=n())
  commAvgLittQual <- commAvgdLittQualMySites %>% group_by(SITE.TYPE) %>% 
    summarise(LITT.CTON.MEAN= mean(SITE.AVG.CN.RATIO, na.rm=T),
              LITT.CTON.SE=
                sd(SITE.AVG.CN.RATIO,na.rm=T)/sqrt(n()))
  commAvgLittQual <- bind_cols(data.frame(NAME=c("Comm. Avg.", "Comm. Avg.")),
                               commAvgLittQual)
}

# Processing and plotting ...
{
  # merge aboveground & belowground data point-wise
  cCycle <- belowGround %>% 
    dplyr::select(-SITE.ID, -POINT.TYPE, -POINT.IN.SITE, 
                  -POINT.NO, -SOIL.WT) %>%
    left_join(aboveGround, ., by="POINT")
  cCycle <- cCycle %>% 
    separate(col="SITE.ID", into=c("SITE.TYPE", "SITE.NO"), 
             sep=2, remove=F) %>%
    dplyr::select(-SITE.NO)
  cCycle$SITE.TYPE <- as.factor(cCycle$SITE.TYPE)
  
  cCycle <- left_join(cCycle, largeTreeDens, by="POINT")
  
  #####################################
  #            Figure 2               #
  #####################################
  
  # Comparing average soil C and SIR (response) values between FR and CT
  soilcMeansPerSite <- cCycle %>% group_by(SITE.ID) %>% 
    summarise(SITE.MEAN = mean(SOIL.C),
              SITE.TYPE=first(SITE.TYPE)) 
  t.test(soilcMeansPerSite$SITE.MEAN~soilcMeansPerSite$SITE.TYPE)
  soilcMeans <- soilcMeansPerSite %>% group_by(SITE.TYPE) %>%
    summarise(MEAN=mean(SITE.MEAN),
              SE=sd(SITE.MEAN)/sqrt(n()))
  ggplot(cCycle, aes(y=SOIL.C, x=SITE.TYPE, fill=SITE.TYPE)) + theme_classic() +
    geom_boxplot() + guides(fill=F) +
    stat_summary(fun.y=mean, geom="point", shape=18, size=6) +
    scale_x_discrete(labels=
                       c("Contiguous \n (N = 6 sites)", "Fragmented \n (N = 6 sites)")) +    
    ylab("Soil %C \n") +
    theme(legend.title=element_blank(),
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 18),
          #axis.title.y = element_text(angle = 0, margin=margin(0,20,0,0)),
          axis.text = element_text(size = 18)) +
    scale_fill_grey(start=0.5)
  #ggsave("figs/SoilC_CTvsFRCompare_YlabFix.png")
  
  #####################################
  #            Figure 3               #
  #####################################
  
  sirMeansPerSite <- cCycle %>% group_by(SITE.ID) %>% 
    summarise(SITE.MEAN = mean(SOIL.SIR, na.rm=T),
              SITE.TYPE = first(SITE.TYPE))
  t.test(sirMeansPerSite$SITE.MEAN~sirMeansPerSite$SITE.TYPE)
  sirMeans <- sirMeansPerSite %>% group_by(SITE.TYPE) %>%
    summarise(MEAN = mean(SITE.MEAN), SE = sd(SITE.MEAN)/sqrt(n()))
  ggplot(cCycle, aes(y=SOIL.SIR, x=SITE.TYPE, fill=SITE.TYPE)) + theme_classic() +
    geom_boxplot() + guides(fill=F) +
    stat_summary(fun.y=mean, geom="point", shape=18, size=6) +
    scale_x_discrete(labels=
                       c("Contiguous \n (N = 6 sites)", "Fragmented \n (N = 6 sites)")) +
    ylab(expression(atop("Soil SIR", (paste(mu, "gm C-", CO[2], " ", gm^{-1}, " ", hr^{-1}))))) +
    theme(legend.title=element_blank(),
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 18),
          #axis.title.y = element_text(angle = 0, margin=margin(0,20,0,0)),
          axis.text = element_text(size = 18)) +
    scale_fill_grey(start=0.5)
  #ggsave("figs/SIR_CTvsFRCompare_YlabFix.png")
  
  #####################################
  #            Figure 6               #
  #####################################
  
  # Comparing average predictor values between FR and CT
  littWtMeansPerSite <- cCycle %>% group_by(SITE.ID) %>% 
    summarise(SITE.MEAN = mean(LITTER.WT),
              SITE.TYPE=first(SITE.TYPE)) 
  t.test(littWtMeansPerSite$SITE.MEAN~littWtMeansPerSite$SITE.TYPE)
  littWtMeans <- littWtMeansPerSite %>% group_by(SITE.TYPE) %>%
    summarise(MEAN=mean(SITE.MEAN), SE=sd(SITE.MEAN)/sqrt(n()))
  ggplot(cCycle, aes(y=LITTER.WT, x=SITE.TYPE, fill=SITE.TYPE)) + theme_classic() +
    geom_boxplot() + guides(fill=F) +
    stat_summary(fun.y=mean, geom="point", shape=18, size=6) +
    scale_x_discrete(labels=
                       c("Contiguous \n (N = 6 sites)", "Fragment \n (N = 6 sites)")) +    
    ylab(expression(atop("Litter weight", (paste(gm, " ", m^{-2}))))) +
    theme(legend.title=element_blank(),
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 24),
          #axis.title.y = element_text(angle = 0, margin=margin(0,20,0,0)),
          axis.text = element_text(size = 18)) +
    scale_fill_grey(start=0.5)
  #ggsave("figs/LittWt_CTvsFRCompare_YlabFix.png")
  
  littCNMeansPerSite <- cCycle %>% group_by(SITE.ID) %>% 
    summarise(SITE.MEAN = mean(LITTER.CN.RATIO),
              SITE.TYPE=first(SITE.TYPE)) 
  t.test(littCNMeansPerSite$SITE.MEAN~littCNMeansPerSite$SITE.TYPE)
  littCNMeans <- littCNMeansPerSite %>% group_by(SITE.TYPE) %>%
    summarise(MEAN=mean(SITE.MEAN), SE=sd(SITE.MEAN)/sqrt(n()))
  ggplot(cCycle, aes(y=LITTER.CN.RATIO, x=SITE.TYPE, fill=SITE.TYPE)) + theme_classic() +
    geom_boxplot() + guides(fill=F) +
    stat_summary(fun.y=mean, geom="point", shape=18, size=6) +
    scale_x_discrete(labels=
                       c("Contiguous \n (N = 6 sites)", "Fragment \n (N = 6 sites)")) +    
    ylab("Litter C/N \n") +
    theme(legend.title=element_blank(),
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 24),
          #axis.title.y = element_text(angle = 0, margin=margin(0,20,0,0)),
          axis.text = element_text(size = 18)) +
    scale_fill_grey(start=0.5)
  #ggsave("figs/LittCN_CTvsFRCompare_YlabFix.png")
  
  basAreaMeansPerSite <- cCycle %>% group_by(SITE.ID) %>% 
    summarise(SITE.MEAN = mean(LITTER.CN.RATIO),
              SITE.TYPE=first(SITE.TYPE)) 
  t.test(basAreaMeansPerSite$SITE.MEAN~basAreaMeansPerSite$SITE.TYPE)
  basAreaMeans <- basAreaMeansPerSite %>% group_by(SITE.TYPE) %>%
    summarise(MEAN=mean(SITE.MEAN), SE=sd(SITE.MEAN)/sqrt(n()))
  ggplot(cCycle, aes(y=BAS.AREA, x=SITE.TYPE, fill=SITE.TYPE)) + theme_classic() +
    geom_boxplot() + guides(fill=F) +
    stat_summary(fun.y=mean, geom="point", shape=18, size=6) +
    scale_x_discrete(labels=
                       c("Contiguous \n (N = 6 sites)", "Fragment \n (N = 6 sites)")) +    
    ylab(expression(atop("Basal area", (paste(m^{2}, " ", ha^{-1}))))) +
    theme(legend.title=element_blank(),
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 24),
          #axis.title.y = element_text(angle = 0, margin=margin(0,20,0,0)),
          axis.text = element_text(size = 18)) +
    scale_fill_grey(start=0.5)
  #ggsave("figs/BasArea_CTvsFRCompare_YlabFix.png")
  
  littPMeansPerSite <- cCycle %>% group_by(SITE.ID) %>% 
    summarise(SITE.MEAN = mean(LITTER.P),
              SITE.TYPE=first(SITE.TYPE)) 
  t.test(littPMeansPerSite$SITE.MEAN~littPMeansPerSite$SITE.TYPE)
  littPMeans <- littPMeansPerSite %>% group_by(SITE.TYPE) %>%
    summarise(MEAN=mean(SITE.MEAN), SE=sd(SITE.MEAN)/sqrt(n()))
  ggplot(cCycle, aes(y=LITTER.P, x=SITE.TYPE, fill=SITE.TYPE)) + theme_classic() +
    geom_boxplot() + guides(fill=F) +
    stat_summary(fun.y=mean, geom="point", shape=18, size=6) +
    scale_x_discrete(labels=
                       c("Contiguous \n (N = 6 sites)", "Fragment \n (N = 6 sites)")) +    
    ylab("Litter %P \n") +
    theme(legend.title=element_blank(),
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 24),
          #axis.title.y = element_text(angle = 0, margin=margin(0,20,0,0)),
          axis.text = element_text(size = 18)) +
    scale_fill_grey(start=0.5)
  #ggsave("figs/LittP_CTvsFRCompare_YlabFix.png")
  
  #####################################
  #            Figure 7               #
  #####################################
  
  # To plot SIR and soil C v/s litter C/N and P, but site-wise with SEs
  soilNLittMeans <- cCycle %>% group_by(SITE.ID) %>% 
    summarise(SITE.TYPE=first(SITE.TYPE),
              SOIL.C.MEAN = mean(SOIL.C),
              SOIL.C.SE = sd(SOIL.C)/sqrt(n()),
              SIR.MEAN = mean(SOIL.SIR, na.rm=T),
              SIR.SE = sd(SOIL.SIR, na.rm=T)/sqrt(n()),
              LITTER.CN.MEAN = mean(LITTER.CN.RATIO),
              LITTER.CN.SE = sd(LITTER.CN.RATIO)/sqrt(n()),
              LITTER.P.MEAN = mean(LITTER.P),
              LITTER.P.SE = sd(LITTER.P)/sqrt(n()),
              LITTER.WT.MEAN = mean(LITTER.WT),
              LITTER.WT.SE = sd(LITTER.WT)/sqrt(n()),
              BAS.AREA.MEAN = mean(BAS.AREA),
              BAS.AREA.SE = sd(BAS.AREA)/sqrt(n()))
  
  
  #####################################
  #            Figure 9               #
  #####################################
  
  # Comparing community avg v/s my sample avg litter qual
  sampleAvgLittQual <- cCycle %>% filter(SITE.ID == "CT1" |       # Kokka 1
                                           SITE.ID == "CT2" |       # Kokka 2
                                           SITE.ID == "FR1" |       # Rudraguppe
                                           SITE.ID == "FR3" |       # Arji
                                           SITE.ID == "FR6") %>%    # Arapattu
    group_by(SITE.ID) %>% 
    summarise(MEAN= mean(LITTER.CN.RATIO, na.rm=T),
              SE=sd(LITTER.CN.RATIO, na.rm=T)/sqrt(n()),
              SITE.TYPE=first(SITE.TYPE)) %>%
    group_by(SITE.TYPE) %>%
    summarise(LITT.CTON.MEAN=mean(MEAN), 
              LITT.CTON.SE=sd(MEAN, na.rm=T)/sqrt(n()))
  levels(sampleAvgLittQual$SITE.TYPE)[levels(sampleAvgLittQual$SITE.TYPE)=="FR"] <-
    "FRAGMENT"
  levels(sampleAvgLittQual$SITE.TYPE)[levels(sampleAvgLittQual$SITE.TYPE)=="CT"] <-
    "CONTROL"
  sampleAvgLittQual <- bind_cols(data.frame(NAME=c("Sample Avg.", "Sample Avg.")),
                                 sampleAvgLittQual)
  littQualCompare <- bind_rows(commAvgLittQual, sampleAvgLittQual)
  ggplot(littQualCompare, aes(y=LITT.CTON.MEAN, x=NAME, fill=SITE.TYPE)) +
    theme_classic() +
    geom_bar(position="dodge", stat="identity") +
    geom_errorbar(aes(ymin=LITT.CTON.MEAN-(LITT.CTON.SE),
                      ymax=LITT.CTON.MEAN+(LITT.CTON.SE)), width=.1,
                  position=position_dodge(0.9)) +
    xlab("Estimate type") + ylab("Average litter C/N \n") +
    theme(legend.title=element_blank(),
          legend.text = element_text(size = 15),
          legend.position = c(0.15, 0.9),
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(),
          axis.line.x = element_line(colour='black',size=0.5,linetype='solid'),
          axis.line.y = element_line(colour='black',size=0.5,linetype='solid'),
          axis.title = element_text(size = 24),
          #axis.title.y = element_text(angle = 0, margin=margin(0,20,0,0)),
          axis.text = element_text(size = 18)) +
    scale_x_discrete(labels=
                       c("Community \n weighted-average", "Sample average")) +
    scale_fill_grey(labels=c("Contiguous", "Fragment"), start=0.5)
  #ggsave("figs/LittQualCompare_CommVsSampleAvg_YlabFix.png")
  
  #####################################
  #            Figure 10              #
  #####################################
  # Size class-wise plotting of tree community
  treeDensS1 <- filter(treeDensBySizeClass, SIZE.CLASS=="S1")
  t.test(treeDensS1$PROP.STEMS ~ treeDensS1$SITE.TYPE)
  treeDensS2 <- filter(treeDensBySizeClass, SIZE.CLASS=="S2")
  t.test(treeDensS2$PROP.STEMS ~ treeDensS2$SITE.TYPE)
  treeDensS3 <- filter(treeDensBySizeClass, SIZE.CLASS=="S3")
  t.test(treeDensS3$PROP.STEMS ~ treeDensS3$SITE.TYPE)
  treeDensS4 <- filter(treeDensBySizeClass, SIZE.CLASS=="S4")
  t.test(treeDensS4$PROP.STEMS ~ treeDensS4$SITE.TYPE)
  treeDensS5 <- filter(treeDensBySizeClass, SIZE.CLASS=="S5")
  t.test(treeDensS5$PROP.STEMS ~ treeDensS5$SITE.TYPE)
  
  t.test(treeDensS1$TREE.DENS ~ treeDensS1$SITE.TYPE)
  t.test(treeDensS2$TREE.DENS ~ treeDensS2$SITE.TYPE)
  t.test(treeDensS3$TREE.DENS ~ treeDensS3$SITE.TYPE)
  t.test(treeDensS4$TREE.DENS ~ treeDensS4$SITE.TYPE)
  t.test(treeDensS5$TREE.DENS ~ treeDensS5$SITE.TYPE)
  
  ggplot(treeDensMeans, aes(y=MEAN, x=SIZE.CLASS, fill=SITE.TYPE)) +
    theme_classic() +
    geom_bar(position="dodge", stat="identity") +
    geom_errorbar(aes(ymin=MEAN-SE, ymax=MEAN+SE), width=.2, position=position_dodge(0.9)) +
    scale_x_discrete(labels=c(paste(binBounds[1]*100,binBounds[2]*100, sep=" - "),
                              paste(binBounds[2]*100,binBounds[3]*100, sep=" - "),
                              paste(binBounds[3]*100,binBounds[4]*100, sep=" - "),
                              paste(binBounds[4]*100,binBounds[5]*100, sep=" - "),
                              paste(binBounds[5]*100,binBounds[6]*100, sep=" - "))) +
    xlab("\n DBH (cm) class") + ylab("Stems/ha \n") + ylim(0,200) +
    theme(legend.position=c(0.85, 0.85),
          legend.title=element_blank(),
          legend.text = element_text(size = 15),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 15)) +
    annotate("text", x = treeDensMeans$SIZE.CLASS[1], 
             y = treeDensMeans$MEAN[1] + (1.5*treeDensMeans$SE[1]),
             label = "*", size=8) +
    annotate("text", x = treeDensMeans$SIZE.CLASS[2], 
             y = treeDensMeans$MEAN[2] + (1.5*treeDensMeans$SE[2]),
             label = "*", size=8) +
    annotate("text", x = treeDensMeans$SIZE.CLASS[3], 
             y = treeDensMeans$MEAN[3] + (1.5*treeDensMeans$SE[3]),
             label = "*", size=8) +
    scale_fill_grey(labels=c("Contiguous", "Fragment"), start=0.5)
  #ggsave("figs/TreeSizeClassDistr_StemsPerHa.png")
  
  
  ggplot(treeStemPropsMeans, aes(y=MEAN, x=SIZE.CLASS, fill=SITE.TYPE)) +
    theme_classic() +
    geom_bar(position="dodge", stat="identity") + 
    geom_errorbar(aes(ymin=MEAN-SE, ymax=MEAN+SE), width=.2,
                  position=position_dodge(0.9)) +
    scale_x_discrete(labels=c(paste(binBounds[1]*100,binBounds[2]*100,sep=" - "),
                              paste(binBounds[2]*100,binBounds[3]*100,sep=" - "),
                              paste(binBounds[3]*100,binBounds[4]*100,sep=" - "),
                              paste(binBounds[4]*100,binBounds[5]*100,sep=" - "),
                              paste(binBounds[5]*100,binBounds[6]*100,sep=" - "))) +
    xlab("\n DBH (cm) class \n (Diameter at breast height)") + 
    ylab("Proportion \n of stems") +
    ylim(0,0.4) +
    theme(legend.position=c(0.85, 0.85),
          legend.title=element_blank(),
          legend.text = element_text(size = 15),
          axis.line.x = element_line(colour='black',size=0.5,linetype='solid'),
          axis.line.y = element_line(colour='black',size=0.5,linetype='solid'),
          axis.title = element_text(size = 18),
          axis.title.y = element_text(angle = 0, margin=margin(0,20,0,0)),
          axis.text = element_text(size = 15)) +
    annotate("text", x = treeStemPropsMeans$SIZE.CLASS[9], 
             y = treeStemPropsMeans$MEAN[9] + (1.5*treeStemPropsMeans$SE[9]),
             label = "**", size=8) +
    scale_fill_grey(labels=c("Contiguous", "Fragment"), start=0.5)
  #ggsave("figs/TreeSizeClassDistr_PropStems_YlabFix.png")
  
  # Scaling variables to zero-mean and unit-std.dev
  cCycle.S <- data.frame(POINT = cCycle$POINT,
                         SITE.TYPE=cCycle$SITE.TYPE,
                         SITE.ID=cCycle$SITE.ID,
                         SOIL.SIR.S=scale(cCycle$SOIL.SIR), 
                         SOIL.C.S=scale(cCycle$SOIL.C),
                         LITTER.CN.RATIO.S=scale(cCycle$LITTER.CN.RATIO),
                         LITTER.P.S=scale(cCycle$LITTER.P),
                         LITTER.WT.S=scale(cCycle$LITTER.WT),
                         BAS.AREA.S=scale(cCycle$BAS.AREA))
  
  
  # Linear mixed effects models of the responses
  library(lme4)
  
  library(MuMIn)
  
  S.soilcLittWtFragOnlyLMM <- lmer(SOIL.C.S ~ LITTER.WT.S + (1|SITE.ID),
                                   data=filter(cCycle.S, SITE.TYPE=="FR"),
                                   REML=FALSE)
  r.squaredGLMM(S.soilcLittWtFragOnlyLMM)
  
  ggplot(cCycle.S[cCycle.S$SITE.TYPE == "FR",], aes(y=SOIL.C.S, x=LITTER.WT.S,
                                                    color=SITE.ID)) + 
    theme_classic() + geom_point(size=3) +
    geom_line(aes(y = predict(S.soilcLittWtFragOnlyLMM)), size = 1) +
    xlab(expression(atop("Litter weight", (paste(gm, " ", m^{-2}))))) +
    ylab("Soil %C \n") +
    theme(plot.title = element_text(size=20),
          legend.position = "none",
          #legend.title=element_blank(),
          #legend.text = element_text(size = 15),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18)) #+
  #scale_color_discrete(labels=as.character(soilNLittMeans$SITE.NAME))
  #ggsave("figs/SoilCvsLitterWt_FragInt_ModelFitPlots.png")
}
