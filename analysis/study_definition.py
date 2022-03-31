from cohortextractor import StudyDefinition, patients, codelist, codelist_from_csv  
from codelists import *

list_dict = {"LeftHemicolectomy":left_hemicolectomy_codes, "RightHemicolectomy":right_hemicolectomy_codes,
"TotalColectomy":total_colectomy_codes,"RectalResection":rectal_resection_codes}

comorb_code_list_dict = {"MI":MI_Read_codes,"CCF":CCF_Read_codes,"Stroke":Stroke_Read_codes,"PVD":PVD_Read_codes,"Dementia":Dementia_Read_codes,
"Respiratory":Chronic_Respiratory_Read_codes,"RA_SLE_Psoriasis":RA_SLE_Psoriasis_Read_codes,"Ulcer_or_bleed":Ulcer_or_bleed_Read_codes,
"all_liver":all_liver_disease_Read_codes,"Cirrhosis":cirrhosis_snomed_codes,"all_diabetes":all_diabetes_Read_codes,
"Diabetic_Complications":diabetic_complication_snomed_codes,"Other_Neurology":other_neuro_Read_codes,"Renal":Chronic_kidney_snomed_codes,"CKD_3_5":CKD_Read_codes,
"Non_Haematology_malignancy":non_haem_cancer,"Haematology_malignancy":haem_cancer_Read_codes,"Metastases":metastatic_cancer_snomed_codes,
"HIV":HIV_Read_codes}

def with_these_vaccination_date_X(name, index_date, n, return_expectations):

    def var_signature(name, on_or_after, return_expectations):
        return {
            name: patients.with_tpp_vaccination_record(
                    returning="date",
                    target_disease_matches="SARS-2 CORONAVIRUS",
                    on_or_after=on_or_after,
                    date_format="YYYY-MM-DD",
                    find_first_match_in_period=True,
                    return_expectations=return_expectations
        ),
        }
    variables = var_signature(f"{name}_1", index_date, return_expectations)
    for i in range(2, n+1):
        variables.update(var_signature(f"{name}_{i}", f"{name}_{i-1} + 1 day", return_expectations))
    return variables


def loop_over_OPCS_codelists(code_list_dict, returning, return_expectations):

    def with_these_procedures(key,codes,returning,return_expectations):
        return {
            f"{key}_{returning}": (
                patients.admitted_to_hospital(
                    with_these_procedures=codes,
                    between=["index_date","index_date + 2 years"],
                    find_first_match_in_period=True,
                    returning=returning,
                    date_format="YYYY-MM-DD",
                    return_expectations=return_expectations,
                )
            )
        }
    variables = {}
    for key,codes in code_list_dict.items():
        variables.update(with_these_procedures(key,codes,returning,return_expectations))
    return variables

def loop_over_OPCS_codelists_discharge(code_list_dict, returning, return_expectations):

    def with_these_procedures(admission_date,key,codes,returning,return_expectations):
        return {
            f"{key}_{returning}": (
                patients.admitted_to_hospital(
                    with_these_procedures=codes,
                    between=[admission_date,"index_date + 2 years"],
                    find_first_match_in_period=True,
                    returning=returning,
                    date_format="YYYY-MM-DD",
                    return_expectations=return_expectations,
                )
            )
        }
    variables = {}
    for key,codes in code_list_dict.items():
        variables.update(with_these_procedures(f"{key}_date_admitted",key,codes,returning,return_expectations))
    return variables

def post_operative_COVID(code_list_dict, returning, return_expectations):

    def with_these_procedures(key,admission_date, admission90_date,returning,return_expectations):
        return {
            f"{key}_{returning}": (
                patients.with_test_result_in_sgss(
                    pathogen = "SARS-CoV-2",
                    returning = returning,
                    test_result = "positive",
                    date_format="YYYY-MM-DD",
                    between=[admission_date, admission90_date],
                    return_expectations = return_expectations,
                )
            )
        }
    variables = {}
    for key,codes in code_list_dict.items():
        variables.update(with_these_procedures(key,f"{key}_date_admitted",f"{key}_date_discharged + 90 days",returning,return_expectations))
    return variables

def with_emergency_readmissions(code_list_dict, returning, return_expectations):

    def with_these_procedures(key,admission_date, admission90_date,returning,return_expectations):
        return {
            f"{key}_emergency_readmit_{returning}": (
                patients.admitted_to_hospital(
                    with_admission_method = ['21', '2A', '22', '23', '24', '25', '2D'],
                    find_first_match_in_period=True,
                    returning = returning,
                    date_format="YYYY-MM-DD",
                    between=[admission_date, admission90_date],
                    return_expectations = return_expectations,
                )
            )
        }
    variables = {}
    for key,codes in code_list_dict.items():
        variables.update(with_these_procedures(key,f"{key}_date_discharged",f"{key}_date_discharged + 90 days",returning,return_expectations))
    return variables

def with_post_op_hospital_events(name,event_list,code_list_dict, returning, return_expectations):

    def with_these_procedures(name,event_list,key,admission_date, admission90_date,returning,return_expectations):
        return {
            f"{key}_{name}_HES_{returning}": (
                patients.admitted_to_hospital(
                    with_these_diagnoses = event_list,
                    find_first_match_in_period=True,
                    returning = returning,
                    date_format="YYYY-MM-DD",
                    between=[admission_date, admission90_date],
                    return_expectations = return_expectations,
                )
            )
        }
    variables = {}
    for key,codes in code_list_dict.items():
        variables.update(with_these_procedures(name,event_list,key,f"{key}_date_admitted",f"{key}_date_discharged + 90 days",returning,return_expectations))
    return variables

def with_post_op_GP_events(name,event_list,code_list_dict, returning, return_expectations):

    def with_these_procedures(name,event_list,key,admission_date, admission90_date,returning,return_expectations):
        return {
            f"{key}_{name}_GP_{returning}": (
                patients.with_these_clinical_events(
                    codelist = event_list,
                    find_first_match_in_period=True,
                    returning = returning,
                    date_format="YYYY-MM-DD",
                    between=[admission_date, admission90_date],
                    return_expectations = return_expectations,
                )
            )
        }
    variables = {}
    for key,codes in code_list_dict.items():
        variables.update(with_these_procedures(name,event_list,key,f"{key}_date_admitted",f"{key}_date_discharged + 90 days",returning,return_expectations))
    return variables

def with_post_op_GP_medications(name,event_list,code_list_dict, returning, return_expectations):

    def with_these_procedures(name,event_list,key,admission_date, admission90_date,returning,return_expectations):
        return {
            f"{key}_{name}_prescriptions_{returning}": (
                patients.with_these_medications(
                    codelist = event_list,
                    find_first_match_in_period=True,
                    returning = returning,
                    date_format="YYYY-MM-DD",
                    between=[admission_date, admission90_date],
                    return_expectations = return_expectations,
                )
            )
        }
    variables = {}
    for key,codes in code_list_dict.items():
        variables.update(with_these_procedures(name,event_list,key,f"{key}_date_admitted",f"{key}_date_discharged + 90 days",returning,return_expectations))
    return variables


def with_pre_op_GP_events(comorb_code_list_dict):

    def with_these_conditions(key,codes):
        return {
            f"pre_{key}_GP": (
                patients.with_these_clinical_events(
                    codelist = codes,
                    find_first_match_in_period=True,
                    returning = "date",
                    date_format="YYYY-MM-DD",
                    between=["1950-01-01", "index_date + 2 years"],
                    return_expectations = {"date": {"earliest": "1950-01-01", "latest": "today"},"rate": "uniform",}
                )
            )
        }
    variables = {}
    for key,codes in comorb_code_list_dict.items():
        variables.update(with_these_conditions(key,codes))
    return variables

study = StudyDefinition(
    index_date = "2020-02-01",

    default_expectations={
        "date": {"earliest": "index_date" , "latest": "index_date + 2 years"},
        "rate": "uniform",
        "incidence": 1,
    },

    #####################################
    # Population definition
    #####################################
   population=patients.satisfying(
       """
        (first_surgery_date)
        AND (first_surgery_date - dob >=(18*365)  AND first_surgery_date - dob <= (110*365))
        AND (sex = "M" OR sex = "F")
       """,
 #      registered=patients.registered_with_one_practice_between(
  #     "first_surgery_date - 365 days", "first_surgery_date"    ## Minimum should have prior registration for first surgery
   #    ),

    #   follow_up=patients.registered_as_of(
     #  "first_surgery_date + 90 days" ## Minimum should have follow up for first surgery
     #  )
   ),

     #   registered
     #   AND (follow_up OR died)
     #   AND

    **loop_over_OPCS_codelists(list_dict,returning = "date_admitted", return_expectations ={"incidence": 1,"rate" : "uniform",}),
   
    first_surgery_date=patients.minimum_of("LeftHemicolectomy_date_admitted", "RightHemicolectomy_date_admitted","TotalColectomy_date_admitted", "RectalResection_date_admitted"),

    **loop_over_OPCS_codelists_discharge(list_dict,returning = "date_discharged", return_expectations ={"incidence": 1,"rate" : "uniform",}),
   
    first_surgery_discharge_date=patients.minimum_of("LeftHemicolectomy_date_discharged", "RightHemicolectomy_date_discharged","TotalColectomy_date_discharged", "RectalResection_date_discharged"),

    region=patients.registered_practice_as_of(
        "first_surgery_date",
        returning="nuts1_region_name",
        return_expectations={
            "rate": "universal",
            "category": {
                "ratios": {
                    "North East": 0.1,
                    "North West": 0.1,
                    "Yorkshire and the Humber": 0.1,
                    "East Midlands": 0.1,
                    "West Midlands": 0.1,
                    "East of England": 0.1,
                    "London": 0.2,
                    "South East": 0.2,
                },
            },
        },
    ),

    dob=patients.date_of_birth(  
        "YYYY-MM",
        return_expectations={
            "date": {"earliest": "1950-01-01", "latest": "today"},
            "rate": "uniform",
        }
    ),

    sex=patients.sex(
        return_expectations={
            "rate": "universal",
            "category": {"ratios": {"M": 0.49, "F": 0.51}},
        }
    ),

    bmi=patients.most_recent_bmi(
        between=["index_date", "first_surgery_date"],
        minimum_age_at_measurement=18,
        include_measurement_date=True,
        date_format="YYYY-MM",
        return_expectations={
            "float": {"distribution": "normal", "mean": 28, "stddev": 8},
            "incidence": 0.80,
        }
    ),
    
    imd=patients.address_as_of(
        'first_surgery_date',
        returning="index_of_multiple_deprivation",
        round_to_nearest=100,
        return_expectations={
            "rate": "universal",
            "category": {"ratios": {"100": 0.1, "200": 0.2, "300": 0.7}},
        },
    ),


    **with_these_vaccination_date_X(
        name = "covid_vaccine_dates", 
        index_date = "index_date",
        n = 3, 
        return_expectations = {
            "rate": "uniform",
            "incidence": 0.85,
        },
    ),

    #########################################
    # Procedure details
    #########################################

    **loop_over_OPCS_codelists(list_dict,returning = "admission_method", return_expectations ={"category": {"ratios": {"11": 0.1, "21": 0.2, "22": 0.7}}, "incidence" : 1}),

    **loop_over_OPCS_codelists(list_dict,returning = "primary_diagnosis", return_expectations ={ "category": {"ratios": {"C180": 0.5, "C190": 0.5}}, "incidence" : 1,}),

    **loop_over_OPCS_codelists(list_dict,returning = "days_in_critical_care", return_expectations ={"category": {"ratios": {"5": 0.1, "6": 0.2, "7": 0.7}}, "incidence" : 0.1}),
   
    #########################################
    # Pre operative risk factors
    ##########################################

 
     MI_HES=patients.admitted_to_hospital(
        between=["index_date","index_date + 2 years"],
        returning="date_admitted",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",
            "incidence": 0.1,
        },
        with_these_diagnoses=MI_HES_codes
    ),

    CCF_HES=patients.admitted_to_hospital(
        between=["index_date","index_date + 2 years"],
        returning="date_admitted",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",
            "incidence": 0.1,
        },
        with_these_diagnoses=CCF_HES_codes
    ),

    **with_pre_op_GP_events(comorb_code_list_dict),
    ####################################
    # Follow up events
    ######################################

    **post_operative_COVID(list_dict, returning= "case_category", return_expectations = {"incidence" : 1,"category": {"ratios": {"": 0.3, "LFT_Only": 0.4, "PCR_Only": 0.2, "LFT_WithPCR": 0.1}},}),

    **post_operative_COVID(list_dict, returning= "date", return_expectations = {"incidence" : 1,"rate" : "uniform"}),

    **with_emergency_readmissions(list_dict, returning = "primary_diagnosis",  return_expectations={"rate": "uniform","category": {"ratios": {"K920": 0.5, "K921": 0.5}}}),

    **with_emergency_readmissions(list_dict, returning = "date_admitted",  return_expectations={"rate" : "uniform"}),

    **with_post_op_hospital_events(name= "VTE",event_list=VTE_HES_codes,code_list_dict=list_dict, returning="date_admitted", return_expectations={"incidence": 0.1,"rate" : "uniform",}),

    **with_post_op_GP_events(name= "VTE",event_list=VTE_Read_codes,code_list_dict=list_dict, returning="date", return_expectations={"incidence": 0.1,"rate" : "uniform",}),

    **with_post_op_GP_medications(name= "anticoagulation",event_list=any_anticoagulation,code_list_dict=list_dict, returning="date", return_expectations={"incidence": 0.1,"rate" : "uniform",}),

   # **with_post_op_GP_events(name= "CurrentCancer", event_list = non_haem_cancer,code_list_dict=list_dict, returning="date", return_expectations={"incidence": 0.1,"rate" : "uniform",}),

    died=patients.died_from_any_cause(
        between=["index_date","index_date + 3 years"],
        returning="binary_flag",
        return_expectations={
            "rate" : "uniform",
            "incidence": 0.1,
        },
    ),
	date_death_ons = patients.died_from_any_cause(
        between=["index_date","index_date + 3 years"],
		returning = "date_of_death",
		date_format = "YYYY-MM-DD",
		return_expectations = {
			"rate": "exponential_increase"
		},
	),

	date_death_cpns = patients.with_death_recorded_in_cpns(
        between=["index_date","index_date + 3 years"],
		returning = "date_of_death",
		date_format = "YYYY-MM-DD",
		return_expectations = {
			"rate": "exponential_increase"
		},
	),

)
