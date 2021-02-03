N_TISSUES_PROJ = {"mc_tm": 12, "mc_ts24": 15, "mc_ts30": 10, "mc_other_10X": 5, "mc_hca": 35, "mc_other": 4,
                  "mc_PanglaoDB": 5}


def get_tissue_by_task_id(dataset, task_id, tasks_per_tiss):
    if dataset == "mc_tm" or dataset == "tm":
        tissues = ("Bladder", "Heart_and_Aorta", "Kidney", "Limb_Muscle", "Liver", "Lung", "Mammary_Gland", "Marrow",
              "Spleen", "Thymus", "Tongue", "Trachea")
    elif dataset == "mc_ts24" or dataset == "ts24":
        tissues = ("BAT", "Bladder", "BM", "GAT", "Heart", "Hepatocytes", "Kidney", "Lung", "MAT", "Muscle",
                   "Pancrease_exo", "SCAT", "Spleen", "Thymus", "Tongue")
    elif dataset == "mc_ts30" or dataset == "ts30":
        tissues = ("Fat", "Heart_and_Aorta", "Kidney", "Large_Intestine", "Limb_Muscle", "Liver", "Lung", "Marrow",
                   "Pancreas", "Spleen")
    elif dataset == "mc_ebi_tm" or dataset == "ebi_tm":
        tissues = ("Adipose", "Bladder", "Bone_Marrow", "Cerebellum", "Cerebral_Cortex", "Colon", "Diaphragm", "Fat",
                   "Heart_and_Aorta", "Hippocampus", "Kidney", "Limb_Muscle", "Liver", "Lung", "Mammary_Gland",
                   "Pancreas", "Skin", "Spleen", "Striatum", "Thymus", "Tongue", "Trachea")
    elif dataset == "mc_ebi" or dataset == "ebi":
        tissues = ("E-ENAD-21-breast_epithelial1", "E-ENAD-27-islet", "E-GEOD-81547-pancreas", "E-GEOD-81608-islet2",
                   "E-GEOD-83139-pancreatic_endocrine", "E-GEOD-86618-lung_epithelial", "E-GEOD-89232-dendritic_cells",
                   "E-GEOD-130148-lung", "E-MTAB-6386-B_cells", "E-MTAB-6653-lung_carcinomas",
                   "E-MTAB-6701-fetal-maternal_interface", "E-CURD-6-bone_marrow", "E-MTAB-7316-retina")
    elif dataset == "mc_mca" or dataset == "mca":
        tissues = ("Bladder", "BoneMarrow", "Brain", "Kidney", "Liver", "Lung", "MammaryGland.Involution",
                   "MammaryGland.Lactation", "MammaryGland.Pregnancy", "MammaryGland.Virgin", "MesenchymalStemCells",
                   "Muscle", "Ovary", "Pancreas", "PeripheralBlood", "Placenta", "Prostate", "SmallIntestine", "Spleen",
                   "Stomach", "Testis", "Thymus", "TrophoblastStemCells", "Uterus")
    elif dataset == "mc_other_10X" or dataset == "other_10X":
        tissues = ("adipose", "ASD_snRNAseq", "liver", "kidney2", "skin")
    elif dataset == "mc_other" or dataset == "other":
        tissues = ("colon-epi", "colon-fib", "colon-imm", "retina")
    elif dataset == "mc_hca" or dataset == "hca":
        tissues = ("Adipose", "Adrenal-Gland", "Artery", "Ascending-Colon", "Bladder", "Bone-Marrow", "Cerebellum",
                    "Cervix", "Duodenum", "Epityphlon", "Esophagus", "Fallopian-Tube", "Gall-Bladder", "Heart", "Ileum",
                    "Jejunum", "Kidney", "Liver", "Lung", "Muscle", "Omentum", "Pancreas", "Peripheral-Blood", "Pleura",
                    "Prostate", "Rectum", "Sigmoid-Colon", "Spleen", "Stomach", "Temporal-Lobe", "Thyroid", "Trachea",
                    "Transverse-Colon", "Ureter", "Uterus")
    elif dataset == "mc_PanglaoDB" or dataset == "PanglaoDB":
        tissues = ("Lung_Epithelial", "Mammary_Gland", "Pancreatic_Islets", "Prostate", "Testis", "Bone_Marrow", "Liver",
                   "Substantia_Nigra")
    elif dataset == "mc_blood" or dataset == "blood":
        tissues = ("Blood",)
    elif dataset == "mc_brain" or dataset == "brain":
        tissues = ("Brain",)
    elif dataset == "mc_brain_olfactory" or dataset == "brain_olfactory":
        tissues = ("Brain_olfactory",)
    elif dataset == "mc_heart_circulation" or dataset == "heart_circulation":
        tissues = ("Heart_Circulation",)
    else:
        return None
    return tissues[task_id // tasks_per_tiss]
