import multicrabDatasetsCommon as common

datasets = {
    # Signal
    # "TTToHplusBWB_M90": {
    #     "dataVersion": ,
    #     "crossSection": 16.442915,
    #     "data": {
    #         "AOD": {
    #             "datasetpath": ,
    #             "number_of_jobs": 20,
    #         },
    #     }
    # },
    # "TTToHplusBWB_M100": {
    #     "dataVersion": ,
    #     "crossSection": 14.057857,
    #     "data": {
    #         "AOD": {
    #             "datasetpath": ,
    #             "number_of_jobs": 20,
    #         },
    #     }
    # },
    # "TTToHplusBWB_M120": {
    #     "dataVersion": ,
    #     "crossSection": 8.984715,
    #     "data": {
    #         "AOD": {
    #             "datasetpath": ,
    #             "number_of_jobs": 20,
    #         },
    #     }
    # },
    # "TTToHplusBWB_M140": {
    #     "dataVersion": ,
    #     "crossSection": 4.223402,
    #     "data": {
    #         "AOD": {
    #             "datasetpath": ,
    #             "number_of_jobs": 20,
    #         },
    #     }
    # },
    # "TTToHplusBWB_M160": {
    #     "dataVersion": ,
    #     "crossSection": 0.811493,
    #     "data": {
    #         "AOD": {
    #             "datasetpath": ,
    #             "number_of_jobs": 20,
    #         },
    #     }
    # },

    # QCD Winter10
    "QCD_Pt30to50_Winter10": {
        "dataVersion": "39Xredigi",
        "crossSection": 5.312e+07,
        "data": {
            "AOD": {
                "datasetpath": "/QCD_Pt_30to50_TuneZ2_7TeV_pythia6/Winter10-E7TeV_ProbDist_2010Data_BX156_START39_V8-v1/AODSIM",
                "number_of_jobs": 150,
            },
        },
    },
    # "QCD_Pt50to80_Winter10": {
    #     "dataVersion":  "39Xredigi",
    #     "crossSection": 6.359e+06,
    #     "data": {
    #         "AOD": {
    #             "datasetpath": 
    #             "number_of_jobs": 150,
    #         },
    #     },
    # },
    "QCD_Pt80to120_Winter10": {
        "dataVersion": "39Xredigi"
        "crossSection": 7.843e+05,
        "data": {
            "AOD": {
                "datasetpath": "/QCD_Pt_80to120_TuneZ2_7TeV_pythia6/Winter10-E7TeV_ProbDist_2010Data_BX156_START39_V8-v1/AODSIM",
                "number_of_jobs": 150,
            },
        },
    },
    "QCD_Pt120to170_Winter10": {
        "dataVersion": "39Xredigi",
        "crossSection": 1.151e+05,
        "data": {
            "AOD": {
                "datasetpath": "/QCD_Pt_120to170_TuneZ2_7TeV_pythia6/Winter10-E7TeV_ProbDist_2010Data_BX156_START39_V8-v1/AODSIM",
                "number_of_jobs": 150,
            },
        },
    },
    "QCD_Pt170to300_Winter10": {
        "dataVersion": "39Xredigi",
        "crossSection": 2.426e+04,
        "data": {
            "AOD": {
                "datasetpath": "/QCD_Pt_170to300_TuneZ2_7TeV_pythia6/Winter10-E7TeV_ProbDist_2010Data_BX156_START39_V8-v1/AODSIM",
                "number_of_jobs": 150,
            },
        },
    },
    "QCD_Pt300to470_Winter10": {
        "dataVersion": "39Xredigi",
        "crossSection": 1.168e+03,
        "data": {
            "AOD": {
                "datasetpath": "/QCD_Pt_300to470_TuneZ2_7TeV_pythia6/Winter10-E7TeV_ProbDist_2010Data_BX156_START39_V8-v1/AODSIM",
                "number_of_jobs": 150 # Adjusted for PATtuple file size
            },
        },
    },


    # # Electroweak (Winter10)
    # "TTJets": {
    #     "dataVersion": 
    #     "crossSection": 165,
    #     "data": {
    #         "AOD": {
    #             "datasetpath":
    #             "number_of_jobs": 100,
    #         },
    #     },
    # },
    "WJets_Winter10": {
        "dataVersion": "38Xrelval",
        "crossSection": 24640,
        "data": {
            "AOD": {
                "datasetpath": "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/Winter10-START39_V8-v1/AODSIM",
                "number_of_jobs": 1000,
                "use_server": 1,
            },
        },
    },


    # # Backgrounds for electroweak background measurement (Fall10)
    # "QCD_Pt20_MuEnriched": {
    #     "dataVersion": 
    #     "crossSection": 296600000.*0.0002855,
    #     "data": {
    #         "AOD": {
    #             "datasetpath":
    #             "number_of_jobs": 200,
    #         }
    #     },
    # },
    # "DYJetsToLL": { # Z+jets
    #     "dataVersion":
    #     "crossSection": 2321,
    #     "data": {
    #         "AOD": {
    #             "datasetpath":
    #             "number_of_jobs": 15,
    #         }
    #     },
    # },
    # "TToBLNu_s-channel": {
    #     "dataVersion": 
    #     "crossSection": 0.99,
    #     "data": {
    #         "AOD": {
    #             "datasetpath":
    #             "number_of_jobs": 10, # Adjusted for PAT on the fly
    #         }
    #     },
    # },
    # "TToBLNu_t-channel": {
    #     "dataVersion":
    #     "crossSection": 63./3.,
    #     "data": {
    #         "AOD": {
    #             "datasetpath":
    #             "number_of_jobs": 10, # Adjusted for PAT on the fly
    #         }
    #     },
    # },
    # "TToBLNu_tW-channel": {
    #     "dataVersion":
    #     "crossSection": 10.56,
    #     "data": {
    #         "AOD": {
    #             "datasetpath":
    #             "number_of_jobs": 10, # Adjusted for PAT on the fly
    #         },
    #     },
    # },
}
