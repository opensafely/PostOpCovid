from cohortextractor import (
    codelist_from_csv,
    codelist,
    filter_codes_by_category,
    combine_codelists
)


############# Broad major surgical categories

any_proc = codelist_from_csv(
    "codelists/user-colincrooks-procedurecategory.csv",
    system="opcs4",
    column="code",
    category_column="Category"
    )        

abdominal_codes = filter_codes_by_category(any_proc,"Abdominal")

cardiac_codes = filter_codes_by_category(any_proc,"Cardiac" )
    
obstetrics_codes = filter_codes_by_category(any_proc,"Obstetrics"    )
    
orthopaedic_codes = filter_codes_by_category(any_proc,"Orthopaedic"  )    
    
thoracic_codes = filter_codes_by_category(any_proc,"Thoracic"  )    
    
vascular_codes = filter_codes_by_category(any_proc,"Vascular")


major_codes = filter_codes_by_category(codelist_from_csv(
    "codelists/user-colincrooks-procedurecategory.csv",
    system="opcs4",
    column="code",
    category_column="Major"
    ),'1')        

############## Colectomy codes


right_hemicolectomy_codes = codelist_from_csv(
  "codelists/user-colincrooks-surgery-right-hemicolectomy-opcs-4-codes.csv",
    column="Code",
    system="opcs4"
)

left_hemicolectomy_codes = codelist_from_csv(
  "codelists/user-colincrooks-surgery-left-hemicolectomy-opcs-4-codes.csv",
    column="Code",
    system="opcs4"
)

total_colectomy_codes = codelist_from_csv(
  "codelists/user-colincrooks-surgery-total-colectomy-opcs-4-codes.csv",
    column="Code",
    system="opcs4"
)

rectal_resection_codes = codelist_from_csv(
  "codelists/user-colincrooks-surgery-rectal-resection-opcs-4-codes.csv",
    column="Code",
    system="opcs4"
)

any_colorectal_resection = combine_codelists(
    left_hemicolectomy_codes,
    right_hemicolectomy_codes,
    total_colectomy_codes,
    rectal_resection_codes
)


cholecystectomy_codes = codelist_from_csv(
  "codelists/user-colincrooks-cholecystectomy.csv",
    column="Code",
    system="opcs4"
)

hipreplacement_codes = codelist_from_csv(
  "codelists/user-colincrooks-hipreplacement.csv",
    column="Code",
    system="opcs4"
)

kneereplacement_codes = codelist_from_csv(
  "codelists/user-colincrooks-kneereplacement.csv",
    column="Code",
    system="opcs4"
)

############################################################################################
# Fracture codelists reused from https://github.com/opensafely/surgery-research/ 


#### Surgical procedures #####

elective_codes = codelist_from_csv(
    "codelists/opensafely-final-orthopaedic-surgery-elective-procedure-opcs-codes-sc.csv", 
    system="opcs4",
    column="Code")

nonelective_codes = codelist_from_csv(
    "codelists/opensafely-final-orthopaedic-surgery-emergency-procedure-opcs-codes-sc.csv", 
    system="opcs4",
    column="Code")

fracture_codes = codelist_from_csv(
    "codelists/opensafely-fracture-potential-emergency-opcs-codes.csv", 
    system="opcs4",
    column="Code")

nonelective_notcomplication_codes = codelist_from_csv(
    "codelists/opensafely-surgery-emergency-opcs-4-codes-without-complication-codes.csv",
    system="opcs4",
    column="Code")

#### Trauma diagnosis codes #####

trauma_dislocation_codes = codelist_from_csv(
    "codelists/user-salmachaudhury-surgery-trauma-icd-10-diagnosis-codes-dislocations-sc.csv", 
    system="icd10",
    column="code")

trauma_spine_codes = codelist_from_csv(
    "codelists/user-salmachaudhury-surgery-trauma-icd-10-diagnosis-codes-spines-sc.csv", 
    system="icd10",
    column="code")

trauma_lower_limb_codes = codelist_from_csv(
    "codelists/user-salmachaudhury-surgery-trauma-icd-10-diagnosis-codes-lower-limb-sc.csv", 
    system="icd10",
    column="code")

trauma_pathological_fracture_codes = codelist_from_csv(
    "codelists/user-salmachaudhury-surgery-trauma-diagnosis-icd-10-codes-pathological-fracture-sc.csv", 
    system="icd10",
    column="code")

trauma_periprosthetic_fracture_codes = codelist_from_csv(
    "codelists/user-salmachaudhury-surgery-trauma-icd-10-diagnosis-codes-peri-prosthetic-fracture-sc.csv", 
    system="icd10",
    column="code")

trauma_neurovascular_injury_codes = codelist_from_csv(
    "codelists/user-salmachaudhury-surgery-trauma-diagnosis-icd-10-neurovascular-injuries.csv", 
    system="icd10",
    column="code")

trauma_open_fracture_codes = codelist_from_csv(
    "codelists/user-salmachaudhury-surgery-trauma-diagnosis-icd-10-codes-open-fractures.csv", 
    system="icd10",
    column="code")

all_trauma_codes = combine_codelists(
        trauma_dislocation_codes,
        trauma_spine_codes,
        trauma_lower_limb_codes,
        trauma_pathological_fracture_codes,
        trauma_periprosthetic_fracture_codes,
        trauma_neurovascular_injury_codes,
        trauma_open_fracture_codes
        )




######################################################################################
VTE_HES_codes = codelist_from_csv(
  "codelists/opensafely-venous-thromboembolism-current-by-type-secondary-care-and-mortality-data.csv", system='icd10', column="code"
)

VTE_Read_codes = codelist_from_csv(
  "codelists/opensafely-incident-venous-thromboembolic-disease.csv",system="ctv3",column="CTV3Code"
)

warfarin = codelist_from_csv(
    "codelists/opensafely-warfarin.csv", system="snomed", column="id"
)
LMWH = codelist_from_csv(
    "codelists/opensafely-low-molecular-weight-heparins.csv", system="snomed", column="code"
)
DOAC = codelist_from_csv(
    "codelists/opensafely-direct-acting-oral-anticoagulants-doac.csv", system="snomed", column="id"
)

any_anticoagulation = combine_codelists(
    LMWH,
    DOAC,
    warfarin
)

MI_HES_codes = codelist_from_csv(
  "codelists/opensafely-cardiovascular-secondary-care.csv", system='icd10', column="icd", category_column = "mi"
)

CCF_HES_codes = codelist_from_csv(
  "codelists/opensafely-cardiovascular-secondary-care.csv", system='icd10', column="icd", category_column = "heartfailure"
)

MI_Read_codes = codelist_from_csv(
  "codelists/opensafely-myocardial-infarction.csv",system="ctv3",column="CTV3ID"
)

CCF_Read_codes = codelist_from_csv(
  "codelists/opensafely-heart-failure.csv", system='ctv3', column="CTV3ID",
)

Stroke_HES_codes = codelist_from_csv(
  "codelists/opensafely-stroke-secondary-care.csv", system='icd10', column="icd"
)

Stroke_Read_codes = codelist_from_csv(
  "codelists/opensafely-stroke-updated.csv",system="ctv3",column="CTV3ID"
)

PVD_Read_codes = codelist_from_csv(
  "codelists/opensafely-peripheral-arterial-disease.csv",system="ctv3",column="code"
)

Dementia_Read_codes = codelist_from_csv(
  "codelists/opensafely-dementia-complete.csv",system="ctv3",column="code"
)

Chronic_Respiratory_Read_codes = codelist_from_csv(
  "codelists/opensafely-chronic-respiratory-disease.csv",system="ctv3",column="CTV3ID"
)

RA_SLE_Psoriasis_Read_codes = codelist_from_csv(
  "codelists/opensafely-ra-sle-psoriasis.csv",system="ctv3",column="CTV3ID"
)

Ulcer_or_bleed_Read_codes = codelist_from_csv(
  "codelists/opensafely-gi-bleed-or-ulcer.csv",system="ctv3",column="CTV3ID"
)

all_liver_disease_Read_codes = codelist_from_csv(
  "codelists/opensafely-chronic-liver-disease.csv",system="ctv3",column="CTV3ID"
)

cirrhosis_snomed_codes = codelist_from_csv(
  "codelists//user-colincrooks-cirrhosis.csv",system="snomed",column="code"
)


all_diabetes_Read_codes = codelist_from_csv(
  "codelists/opensafely-diabetes.csv",system="ctv3",column="CTV3ID"
)

diabetic_complication_snomed_codes = codelist_from_csv(
  "codelists/user-colincrooks-diabetic-complications.csv", system="snomed", column="code"
)

other_neuro_Read_codes = codelist_from_csv(
  "codelists/opensafely-other-neurological-conditions.csv",system="ctv3",column="CTV3ID"
)

Chronic_kidney_snomed_codes = codelist_from_csv(
  "codelists/opensafely-chronic-kidney-disease-snomed.csv",system="snomed",column="id"
)

CKD_Read_codes = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-ckd35.csv",system="snomed",column="code"
)

haem_cancer_Read_codes = codelist_from_csv(
    "codelists/opensafely-haematological-cancer.csv", system="ctv3", column="CTV3ID",
)

lung_cancer_Read_codes = codelist_from_csv(
    "codelists/opensafely-lung-cancer.csv", system="ctv3",column="CTV3ID",
)

other_cancer_Read_codes = codelist_from_csv(
    "codelists/opensafely-cancer-excluding-lung-and-haematological.csv", system="ctv3", column="CTV3ID",
)

non_haem_cancer = combine_codelists(
    lung_cancer_Read_codes,
    other_cancer_Read_codes
)

metastatic_cancer_snomed_codes = codelist_from_csv(
    "codelists/user-colincrooks-metastases.csv", system="snomed", column="code",
)

HIV_Read_codes = codelist_from_csv(
    "codelists/opensafely-hiv.csv", system="ctv3", column="CTV3ID",
)
