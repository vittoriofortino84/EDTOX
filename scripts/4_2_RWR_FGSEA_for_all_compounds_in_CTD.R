
source("functions/pipeline.R")                                                 # pipeline functions
load('outputData/toxdb2gene_final.RData')                                      # pathways file
load("outputData/chem2gene_no_out.RData")                                      # seeds from all compounds
tox2db<-c(go2gene,kegg2gene,reactome2gene,wikiptw2gene,msigdb2gene)            # tox2db is all pathways

library(BiocParallel)

# Drug matrix -------------------------------------------------------------

load("outputData/network/wTO_DM.RData")
#1. DM_hepatocytes.1
p<-bpstart(MulticoreParam(40))
fgsea_res<-pipeline(wTO_DM[[1]],tox2db,0.1,1000,chem2gene,p)
bpstop(p)
save(fgsea_res,file = 'outputData/all_compounds_fgsea/DM_hep.RData')
gc()

#2. DM_liver.1
p<-bpstart(MulticoreParam(40))
fgsea_res<-pipeline(wTO_DM[[2]],tox2db,0.1,1000,chem2gene,p)
bpstop(p)
save(fgsea_res,file = 'outputData/all_compounds_fgsea/DM1.RData')
gc()

#3. DM_liver.3
p<-bpstart(MulticoreParam(40))
fgsea_res<-pipeline(wTO_DM[[3]],tox2db,0.05,1000,chem2gene,p)
bpstop(p)
save(fgsea_res,file = 'outputData/all_compounds_fgsea/DM3.RData')
gc()

#4. DM_liver.5
p<-bpstart(MulticoreParam(40))
fgsea_res<-pipeline(wTO_DM[[4]],tox2db,0.05,1000,chem2gene,p)
bpstop(p)
save(fgsea_res,file = 'outputData/all_compounds_fgsea/DM5.RData')
gc()


# Tg gate -----------------------------------------------------------------

load("outputData/network/wTO_TGG.RData")
#1.TG-gate-human-innvitro
p<-bpstart(MulticoreParam(40))
fgsea_res<-pipeline(wTO_TGG[[1]],tox2db,0.02,700,chem2gene,p)
bpstop(p)
save(fgsea_res,file = 'outputData/all_compounds_fgsea/TGG_human_in_vitro.RData')
gc()

#2.TG-gate-rat-innvitro
p<-bpstart(MulticoreParam(40))
fgsea_res<-pipeline(wTO_TGG[[2]],tox2db,0.05,1000,chem2gene,p)
bpstop(p)
save(fgsea_res,file = 'outputData/all_compounds_fgsea/TGG_rat_in_vitro.RData')
gc()

#3.TGG_rat_in_vivo_rep.15 day
p<-bpstart(MulticoreParam(40))
fgsea_res<-pipeline(wTO_TGG[[3]],tox2db,0.02,200,chem2gene,p)
bpstop(p)
save(fgsea_res,file = 'outputData/all_compounds_fgsea/TGG_rat_in_vivo_rep_15_day.RData')
gc()

#4.TGG_rat_in_vivo_rep.29 day
p<-bpstart(MulticoreParam(40))
fgsea_res<-pipeline(wTO_TGG[[4]],tox2db,0.05,500,chem2gene,p)
bpstop(p)
save(fgsea_res,file = 'outputData/all_compounds_fgsea/TGG_rat_in_vivo_rep_29_day.RData')
gc()

#5.TGG_rat_in_vivo_rep.8 day
p<-bpstart(MulticoreParam(40))
fgsea_res<-pipeline(wTO_TGG[[5]],tox2db,0.05,700,chem2gene,p)
bpstop(p)
save(fgsea_res,file = 'outputData/all_compounds_fgsea/TGG_rat_in_vivo_rep_8_day.RData')
gc()

#6. "TGG_rat_in_vivo_single.High"
p<-bpstart(MulticoreParam(40))
fgsea_res<-pipeline(wTO_TGG[[6]],tox2db,0.05,700,chem2gene,p)
bpstop(p)
save(fgsea_res,file = 'outputData/all_compounds_fgsea/Fgsea_tg_single.high.RData')

#7."TGG_rat_in_vivo_single.Low" 
p<-bpstart(MulticoreParam(40))
fgsea_res<-pipeline(wTO_TGG[[7]],tox2db,0.05,1000,chem2gene,p)
bpstop(p)
save(fgsea_res,file = 'outputData/all_compounds_fgsea/Fgsea_tg_single.low.RData')
gc()

#8."TGG_rat_in_vivo_single.Middle"
p<-bpstart(MulticoreParam(40))
fgsea_res<-pipeline(wTO_TGG[[8]],tox2db,0.05,700,chem2gene,p)
bpstop(p)
save(fgsea_res,file = 'outputData/all_compounds_fgsea/Fgsea_tg_single.middle.RData')

# Consensus_LINCS and hepatocytes -----------------------------------------

load('outputData/network/wTO_hep_cons.RData')
#1. Consensus hepatocytes
p<-bpstart(MulticoreParam(40))
fgsea_res<-pipeline_consensus(wTO_hep_cons,tox2db,0.05,1000,chem2gene,p)
bpstop(p)
save(fgsea_res,file = 'outputData/all_compounds_fgsea/hep_cons.RData')
gc()

load('outputData/network/wTO_LINCS_cons.RData')
#2. Consensus LINCS
p<-bpstart(MulticoreParam(40))
fgsea_res<-pipeline_consensus(wTO_LINCS_cons,tox2db,0.1,700,chem2gene,p)
bpstop(p)
save(fgsea_res,file = 'outputData/all_compounds_fgsea/lincs_cons.RData')
gc()


# PPI ---------------------------------------------------------------------



#1.
p<-bpstart(MulticoreParam(40))
fgsea_res<-pipeline_ppi(tox2db,0.85,1000,chem2gene,p)
bpstop(p)
save(fgsea_res,file = 'outputData/all_compounds_fgsea/ppi.RData')
gc()





