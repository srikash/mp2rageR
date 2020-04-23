library(mp2rageR)
library(neurobase)
library(pracma)

setwd("~/Desktop/MP2RAGE_Function/data/full_FoV")
## Sa2RAGE protocol info and loading the Sa2RAGE data for B1 estimation

file_B1 = "sub01_sa2rage_wip654_wholebrain_B1map_reg-to-mp2rage.nii"
baseResolution = 128
PartialFourierInPE = 6/8
iPATpe = 2
RefLines = 24

paramList_Sa2RAGE = list(
  TR = 2.4,
  TRFLASH = 2.2e-3,
  TIs = c(58e-3, 1800e-3),
  NZslices = c(0,0),
  FlipDegrees = c(4, 10),
  averageT1 = 1.5
)
paramList_Sa2RAGE$NZslices=baseResolution*c(PartialFourierInPE-0.5,0.5)/iPATpe+c(RefLines/2,RefLines/2)*(1-1/iPATpe)

nii_B1=readnii(file_B1)


## MP2RAGE protocol info and loading the MP2RAGE dataset
file_UNI="sub01_mp2rage_1p0mm_iso_p3_UNI_Images.nii"
slicesPerSlab = 240
PartialFourierInSlice = 8/8

paramList_MP2RAGE <-
  list(
    B0 = 7.0,
    TR = 5.0,
    TRFLASH = 6.9e-3,
    TIs = c(800e-3, 2700e-3),
    NZslices = c(0, 0),
    FlipDegrees = c(4, 5)
  )
paramList_MP2RAGE$NZslices=slicesPerSlab * c(PartialFourierInSlice-0.5,0.5)


# check the properties of this MP2RAGE protocol this happens to be a
# very B1 insensitive protocol
#plotMP2RAGEproperties(MP2RAGE)

# load the MP2RAGE data
nii_UNI=readnii(file_UNI)

# performing the correction
# data_B1 = nii_B1@.Data
# data_UNI = nii_UNI@.Data

T1_B1correct(NULL,nii_B1@.Data,paramList_Sa2RAGE,nii_UNI@.Data,NULL,paramList_MP2RAGE,NULL,0.96)

# saving the data
save_untouch_nii(MP2RAGEcorr,'data/MP2RAGEcorr.nii')
save_untouch_nii(T1corrected,'data/T1corrected.nii')
%% if another technique was used to obtain the relative B1 maps
%  (1 means correctly calibrated B1)
B1=load_untouch_nii('Sa2RAGE_B1map.nii.gz');
B1.img=double(B1.img)/1000;

[ T1corrected MP2RAGEcorr] = T1B1correctpackageTFL(B1,MP2RAGEimg,[],MP2RAGE,[],0.96)

% saving the data
save_untouch_nii(MP2RAGEcorr,'MP2RAGEcorr.nii')
save_untouch_nii(T1corrected,'T1corrected.nii')
