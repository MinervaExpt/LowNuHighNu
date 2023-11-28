from collections import OrderedDict

error_bands = OrderedDict()
error_bands["Flux"]                   = ["Flux"]
error_bands["Muon"]                   = ["Muon_Energy_MINERvA",
                                         "Muon_Energy_MINOS",
                                         "Muon_Energy_Resolution",
                                         "Muon_Energy",
                                         "Muon_Theta",
                                         "BeamAngleX",
                                         "BeamAngleY"]
#error_bands["Muon_MINERvA"]           = ["Muon_Energy_MINERvA"]
#error_bands["Muon_MINOS"]             = ["Muon_Energy_MINOS"]
#error_bands["Muon_Other"]             = ["Muon_Energy_Resolution",
#                                         "BeamAngleX",
#                                         "BeamAngleY"]
error_bands["normalization"]          = ["NOMAD_normalization"]
error_bands["HighNu_Stats"]           = ["HighNu_Stat_Bin_01",
                                         "HighNu_Stat_Bin_02",
                                         "HighNu_Stat_Bin_03",
                                         "HighNu_Stat_Bin_04",
                                         "HighNu_Stat_Bin_05",
                                         "HighNu_Stat_Bin_06",
                                         "HighNu_Stat_Bin_07",
                                         "HighNu_Stat_Bin_08",
                                         "HighNu_Stat_Bin_09",
                                         "HighNu_Stat_Bin_10",
                                         "HighNu_Stat_Bin_11",
                                         "HighNu_Stat_Bin_12",
                                         "HighNu_Stat_Bin_13",
                                         "HighNu_Stat_Bin_14",
                                         "HighNu_Stat_Bin_15",
                                         "HighNu_Stat_Bin_16",
                                         "HighNu_Stat_Bin_17",
                                         "HighNu_Stat_Bin_18",
                                         "HighNu_Stat_Bin_19",
                                         "HighNu_Stat_Bin_20",
                                         "HighNu_Stat_Bin_21",
                                         "HighNu_Stat_Bin_22",
                                         "HighNu_Stat_Bin_23",
                                         "HighNu_Stat_Bin_24"]
error_bands["Model_Fits"]             = ["Low_Recoil_2p2h_Tune",
                                         "RPA_HighQ2",
                                         "RPA_LowQ2",
                                         "CCQE_3D_Fit1",
                                         "CCQE_3D_Fit2",
                                         "CCQE_3D_Fit3"]
#error_bands["Low_recoil_fits"]        = ["Low_Recoil_2p2h_Tune"]
#error_bands["RPA"]                    = ["RPA_HighQ2",
#                                         "RPA_LowQ2"]
#error_bands["CCQE_3D_Fits"]           = ["CCQE_3D_Fit1",
#                                         "CCQE_3D_Fit2", 
#                                         "CCQE_3D_Fit3"]
#error_bands["3DFits_Region1"]   = ["CCQE_3D_Fit1"]
#error_bands["3DFits_Region2"]   = ["CCQE_3D_Fit2"]
#error_bands["3DFits_Region3"]   = ["CCQE_3D_Fit3"]
error_bands["Detector"]               = ["Proton_Response",
                                         "Pion_Response",
                                         "EM_Response",
                                         "Other_Response",
                                         "MINOS_Reconstruction_Efficiency"]
#error_bands["Response"]               = ["Proton_Response",
#                                         "Pion_Response",
#                                         "EM_Response",
#                                         "Other_Response"]
#error_bands["MINOS"]                  = ["MINOS_Reconstruction_Efficiency"]
error_bands["CCQENu Other"]           = ["MEU_Proton",
                                         "Mass_Model_Proton",
                                         "Target_Mass",
                                         "Bethe_Bloch",
                                         "Birks_Response_Proton",
                                         "Crosstalk",
#                                         "Reweight_Neutron",
                                         "Michel_Efficiency",  ## CCQENu xsec only (not in highnu flux)
                                         "Proton_TrackEff",    ## CCQENu xsec only (not in highnu flux)
#                                         "Crosstalk"]
#error_bands["CCQENu Reweight"]        = ["Reweight_Neutron",
                                         "Reweight_Pion",
                                         "Reweight_Proton"]
error_bands["GENIE"]                  = ["GENIE_FrAbs_N",
                                         "GENIE_FrAbs_N",
                                         "GENIE_FrAbs_pi",
                                         "GENIE_FrCEx_N",
                                         "GENIE_FrCEx_pi",
                                         "GENIE_FrElas_N",
                                         "GENIE_FrElas_pi",
                                         "GENIE_FrInel_N",
                                         "GENIE_FrInel_pi",
                                         "GENIE_FrPiProd_N",
                                         "GENIE_FrPiProd_pi",
                                         "GENIE_MFP_N",
                                         "GENIE_MFP_pi",
                                         "GENIE_AGKYxF1pi",
                                         "GENIE_AhtBY",
                                         "GENIE_BhtBY",
                                         "GENIE_CCQEPauliSupViaKF",
                                         "GENIE_CV1uBY",
                                         "GENIE_CV2uBY",
                                         "GENIE_EtaNCEL",
                                         "GENIE_MaCCQE",
                                         "GENIE_MaCCQEshape",
                                         "GENIE_MaNCEL",
                                         "GENIE_MaRES",
                                         "GENIE_MvRES",
                                         "GENIE_NormCCQE",
                                         "GENIE_NormCCRES",
                                         "GENIE_NormDISCC",
                                         "GENIE_NormNCRES",
                                         "GENIE_RDecBR1gamma",
                                         "GENIE_Rvn1pi",
                                         "GENIE_Rvn2pi",
                                         "GENIE_Rvn3pi",
                                         "GENIE_Rvp1pi",
                                         "GENIE_Rvp2pi",
                                         "GENIE_Theta_Delta2Npi",
                                         "GENIE_VecFFCCQEshape"]
#error_bands["GENIE_FSI"]              = ["GENIE_FrAbs_N",
#                                         "GENIE_FrAbs_N",
#                                         "GENIE_FrAbs_pi",
#                                         "GENIE_FrCEx_N",
#                                         "GENIE_FrCEx_pi",
#                                         "GENIE_FrElas_N",
#                                         "GENIE_FrElas_pi",
#                                         "GENIE_FrInel_N",
#                                         "GENIE_FrInel_pi",
#                                         "GENIE_FrPiProd_N",
#                                         "GENIE_FrPiProd_pi",
#                                         "GENIE_MFP_N",
#                                         "GENIE_MFP_pi"]
#error_bands["GENIE_InteractionModel"] = ["GENIE_AGKYxF1pi",
#                                         "GENIE_AhtBY",
#                                         "GENIE_BhtBY",
#                                         "GENIE_CCQEPauliSupViaKF",
#                                         "GENIE_CV1uBY",
#                                         "GENIE_CV2uBY",
#                                         "GENIE_EtaNCEL",
#                                         "GENIE_MaCCQE",
#                                         "GENIE_MaCCQEshape",
#                                         "GENIE_MaNCEL",
#                                         "GENIE_MaRES",
#                                         "GENIE_MvRES",
#                                         "GENIE_NormCCQE",
#                                         "GENIE_NormCCRES",
#                                         "GENIE_NormDISCC",
#                                         "GENIE_NormNCRES",
#                                         "GENIE_RDecBR1gamma",
#                                         "GENIE_Rvn1pi",
#                                         "GENIE_Rvn2pi",
#                                         "GENIE_Rvn3pi",
#                                         "GENIE_Rvp1pi",
#                                         "GENIE_Rvp2pi",
#                                         "GENIE_Theta_Delta2Npi",
#                                         "GENIE_VecFFCCQEshape"]

## GENIE-Interaction-Models-Zoom1
#error_bands["GENIE_AGKYxF1pi"]          = ["GENIE_AGKYxF1pi"]
#error_bands["GENIE_AhtBY"]              = ["GENIE_AhtBY"]
#error_bands["GENIE_BhtBY"]              = ["GENIE_BhtBY"]
#error_bands["GENIE_CCQEPauliSupViaKF"]  = ["GENIE_CCQEPauliSupViaKF"]
#error_bands["GENIE_CV1uBY"]             = ["GENIE_CV1uBY"]
#error_bands["GENIE_CV2uBY"]             = ["GENIE_CV2uBY"]
#error_bands["GENIE_EtaNCEL"]            = ["GENIE_EtaNCEL"]
#error_bands["GENIE_MaCCQE"]             = ["GENIE_MaCCQE"]
 
## GENIE-Interaction-Models-Zoom1
#error_bands["GENIE_MaCCQE"]             = ["GENIE_MaCCQE"]
#error_bands["GENIE_MaCCQEshape"]        = ["GENIE_MaCCQEshape"]
#error_bands["GENIE_MaNCEL"]             = ["GENIE_MaNCEL"]
#error_bands["GENIE_MaRES"]              = ["GENIE_MaRES"]
#error_bands["GENIE_MvRES"]              = ["GENIE_MvRES"]
#error_bands["GENIE_NormCCQE"]           = ["GENIE_NormCCQE"]
#error_bands["GENIE_NormCCRES"]          = ["GENIE_NormCCRES"]
#error_bands["GENIE_NormDISCC"]          = ["GENIE_NormDISCC"]
#error_bands["GENIE_NormNCRES"]          = ["GENIE_NormNCRES"]

## GENIE-Interaction-Models-Zoom3
#error_bands["GENIE_RDecBR1gamma"]       = ["GENIE_RDecBR1gamma"]
#error_bands["GENIE_Rvn1pi"]             = ["GENIE_Rvn1pi"]
#error_bands["GENIE_Rvn2pi"]             = ["GENIE_Rvn2pi"]
#error_bands["GENIE_Rvn3pi"]             = ["GENIE_Rvn3pi"]
#error_bands["GENIE_Rvp1pi"]             = ["GENIE_Rvp1pi"]
#error_bands["GENIE_Rvp2pi"]             = ["GENIE_Rvp2pi"]
#error_bands["GENIE_Theta_Delta2Npi"]    = ["GENIE_Theta_Delta2Npi"]
#error_bands["GENIE_VecFFCCQEshape"]     = ["GENIE_VecFFCCQEshape"]

