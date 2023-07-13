from cohortextractor import StudyDefinition, patients, codelist, codelist_from_csv,params  
from codelists import *

name_op = params["name_op"]
indexdate = params["indexdate"]

list_dict = {"Abdominal":abdominal_codes, 
"Cardiac":cardiac_codes, 
"Obstetrics":obstetrics_codes, 
"Orthopaedic":orthopaedic_codes, 
"Thoracic":thoracic_codes,
"Vascular":vascular_codes}

sub_proc_dict = {"Colectomy":any_colorectal_resection,
"HipReplacement":hipreplacement_codes,
"KneeReplacement":kneereplacement_codes,
"Cholecystectomy":cholecystectomy_codes}

comorb_code_list_dict = {"MI":MI_Read_codes,"CCF":CCF_Read_codes,"Stroke":Stroke_Read_codes,"PVD":PVD_Read_codes,"Dementia":Dementia_Read_codes,
"Respiratory":Chronic_Respiratory_Read_codes,"RA_SLE_Psoriasis":RA_SLE_Psoriasis_Read_codes,"Ulcer_or_bleed":Ulcer_or_bleed_Read_codes,
"all_liver":all_liver_disease_Read_codes,"Cirrhosis":cirrhosis_snomed_codes,"all_diabetes":all_diabetes_Read_codes,
"Diabetic_Complications":diabetic_complication_snomed_codes,"Other_Neurology":other_neuro_Read_codes,"Renal":Chronic_kidney_snomed_codes,"CKD_3_5":CKD_Read_codes,
"Non_Haematology_malignancy":non_haem_cancer,"Haematology_malignancy":haem_cancer_Read_codes,"Metastases":metastatic_cancer_snomed_codes,
"HIV":HIV_Read_codes}

n_op = 5

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


def loop_over_OPCS_codelists(code_list_dict, return_expectations):
    def with_these_procedures(key,start,codes,returning,return_expectations,i):
        return {
            f"{key}_{i}_{returning}": (
                patients.admitted_to_hospital(
                    with_these_procedures=codes,
                    between=[start,"index_date + 3 years"],
                    find_first_match_in_period=True,
                    returning=returning,
                    date_format="YYYY-MM-DD",
                    return_expectations=return_expectations,
                )
            )
        }
    
    variables = {}
    for key,codes in code_list_dict.items():
        variables.update(with_these_procedures(key,"index_date",codes,"date_admitted",return_expectations,1))
        variables.update(with_these_procedures(key,f"{key}_1_date_admitted",codes,"date_discharged",return_expectations,1))
        for i in range(2,n_op + 1):
          variables.update(with_these_procedures(key,f"{key}_{i-1}_date_discharged",codes,"date_admitted",return_expectations,i))
          variables.update(with_these_procedures(key,f"{key}_{i}_date_admitted",codes,"date_discharged",return_expectations,i))
    return variables

def loop_over_OPCS_codelists_admission_info(code_list_dict, returning, return_expectations):
    def with_these_procedures(admission_date,discharge_date,key,codes,returning,return_expectations,i):
        return {
            f"{key}_{i}_{returning}": (
                patients.admitted_to_hospital(
                    with_these_procedures=codes,
                    between=[admission_date,discharge_date],
                    find_first_match_in_period=True,
                    returning=returning,
                    date_format="YYYY-MM-DD",
                    return_expectations=return_expectations,
                )
            )
        }
    variables = {}
    for key,codes in code_list_dict.items():
      for i in range(1,n_op + 1):
        variables.update(with_these_procedures(f"{key}_{i}_date_admitted",f"{key}_{i}_date_discharged",key,codes,returning,return_expectations,i))
    return variables

def loop_over_OPCS_codelists_major_minor(code_list_dict, returning, return_expectations):
    def with_these_procedures(admission_date,discharge_date,key,major_codes,returning,return_expectations,i):
        return {
            f"{key}_{i}_Major_HES_{returning}": (
                patients.admitted_to_hospital(
                    with_these_procedures=major_codes,
                    between=[admission_date,discharge_date],
                    find_first_match_in_period=True,
                    returning=returning,
                    date_format="YYYY-MM-DD",
                    return_expectations=return_expectations,
                )
            )
        }
    variables = {}
    for key,codes in code_list_dict.items():
      for i in range(1,n_op + 1):
        variables.update(with_these_procedures(f"{key}_{i}_date_admitted",f"{key}_{i}_date_discharged",key,major_codes,returning,return_expectations,i))
    return variables

def post_operative_COVID(code_list_dict, returning, return_expectations):
    def with_these_procedures(key,admission_date, admission90_date,returning,return_expectations,i):
        return {
            f"{key}_{i}_{returning}": (
                patients.with_test_result_in_sgss(
                    pathogen = "SARS-CoV-2",
                    returning = returning,
                    restrict_to_earliest_specimen_date=False,
                    find_first_match_in_period=True,                    
                    test_result = "positive",
                    date_format="YYYY-MM-DD",
                    between=[admission_date, admission90_date],
                    return_expectations = return_expectations,
                )
            )
        }
    variables = {}
    for key,codes in code_list_dict.items():
      for i in range(1,n_op + 1):
        variables.update(with_these_procedures(key,f"{key}_{i}_date_admitted - 7 days",f"{key}_{i}_date_discharged + 90 days",returning,return_expectations,i))
    return variables


def recent_COVID(code_list_dict, returning, return_expectations):
    def with_these_procedures(key,admission_date, admission90_date,returning,return_expectations,i):
        return {
            f"{key}_{i}_recent_{returning}": (
                patients.with_test_result_in_sgss(
                    pathogen = "SARS-CoV-2",
                    returning = returning,
                    restrict_to_earliest_specimen_date=False,
                    find_last_match_in_period=True,                    
                    test_result = "positive",
                    date_format="YYYY-MM-DD",
                    between=[admission_date, admission90_date],
                    return_expectations = return_expectations,
                )
            )
        }
    variables = {}
    for key,codes in code_list_dict.items():
      for i in range(1,n_op + 1):
        variables.update(with_these_procedures(key,f"{key}_{i}_date_admitted - 42 days",f"{key}_{i}_date_admitted - 8 days",returning,return_expectations,i))
    return variables


def previous_COVID(code_list_dict, returning, return_expectations):
    def with_these_procedures(key,admission_date,returning,return_expectations,i):
        return {
            f"{key}_{i}_previous_{returning}": (
                patients.with_test_result_in_sgss(
                    pathogen = "SARS-CoV-2",
                    returning = returning,
                    restrict_to_earliest_specimen_date=False,
                    find_last_match_in_period=True,                    
                    test_result = "positive",
                    date_format="YYYY-MM-DD",
                    on_or_before=admission_date,
                    return_expectations = return_expectations,
                )
            )
        }
    variables = {}
    for key,codes in code_list_dict.items():
      for i in range(1,n_op + 1):
        variables.update(with_these_procedures(key,f"{key}_{i}_date_admitted - 43 days",returning,return_expectations,i))
    return variables


def with_emergency_readmissions(code_list_dict, returning, return_expectations):
    def with_these_procedures(key,admission_date, admission90_date,returning,return_expectations,i):
        if returning.__eq__("date_admitted"):
            return {
                 f"{key}_{i}_emergency_readmit_{returning}": (
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
        else:
            return{
                f"{key}_{i}_emergency_readmit_{returning}": (
                    patients.admitted_to_hospital(
                        with_admission_method = ['21', '2A', '22', '23', '24', '25', '2D'],
                        find_first_match_in_period=True,
                        returning = returning,
                        between=[admission_date, admission90_date],
                        return_expectations = return_expectations,
                    )
                )
            }
    variables = {}
    for key,codes in code_list_dict.items():            
      for i in range(1,n_op + 1):
        variables.update(with_these_procedures(key,f"{key}_{i}_date_discharged",f"{key}_{i}_date_discharged + 90 days",returning,return_expectations,i))
    return variables

def with_post_op_hospital_events(name,event_list,code_list_dict, returning, return_expectations):
    def with_these_procedures(name,event_list,key,admission_date, admission90_date,returning,return_expectations,i):
        return {
            f"{key}_{i}_{name}_HES_{returning}": (
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
      for i in range(1,n_op + 1):
        variables.update(with_these_procedures(name,event_list,key,f"{key}_{i}_date_admitted",f"{key}_{i}_date_discharged + 90 days",returning,return_expectations,i))
    return variables

def with_sub_procedures(subproc_list,code_list_dict, returning, return_expectations):
    def with_these_procedures(subkey,subcodes,key,admission_date, admission90_date,returning,return_expectations,i):
        return {
            f"{key}_{i}_{subkey}_HES_{returning}": (
                patients.admitted_to_hospital(
                    with_these_procedures = subcodes,
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
      for i in range(1,n_op + 1):
        for subkey,subcodes in subproc_list.items():
          variables.update(with_these_procedures(subkey,subcodes,key,f"{key}_{i}_date_admitted",f"{key}_{i}_date_admitted",returning,return_expectations,i))
    return variables

def with_post_op_GP_events(name,event_list,code_list_dict, returning, return_expectations):
    def with_these_procedures(name,event_list,key,admission_date, admission90_date,returning,return_expectations,i):
        return {
            f"{key}_{i}_{name}_GP_{returning}": (
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
      for i in range(1,n_op + 1):
        variables.update(with_these_procedures(name,event_list,key,f"{key}_{i}_date_admitted",f"{key}_{i}_date_discharged + 90 days",returning,return_expectations,i))
    return variables

def with_post_op_GP_medications(name,event_list,code_list_dict, returning, return_expectations):
    def with_these_procedures(name,event_list,key,admission_date, admission90_date,returning,return_expectations,i):
        return {
            f"{key}_{i}_{name}_prescriptions_{returning}": (
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
      for i in range(1,n_op + 1):
        variables.update(with_these_procedures(name,event_list,key,f"{key}_{i}_date_admitted",f"{key}_{i}_date_discharged + 90 days",returning,return_expectations,i))
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
                    between=["1950-01-01", "index_date + 3 years"],
                    return_expectations = {"date": {"earliest": "1950-01-01", "latest": "today"},"rate": "uniform","incidence": 0.1}
                )
            )
        }
    variables = {}
    for key,codes in comorb_code_list_dict.items():
        variables.update(with_these_conditions(key,codes))
    return variables

study = StudyDefinition(
    index_date = indexdate,

    default_expectations={
        "date": {"earliest": "index_date" , "latest": "index_date + 3 years"},
        "rate": "uniform",
        "incidence": 1,
    },

#        
#        AND registered
#        AND (follow_up OR died)
    #####################################
    # Population definition
    #####################################
    population=patients.satisfying(
        """
        has_surgery
        AND age >= 18  
        AND age <= 110
        AND (sex = "M" OR sex = "F")
        """,
       
        has_surgery=patients.admitted_to_hospital(
            with_these_procedures=any_proc, 
          #  on_or_after = "index_date")
            between = ["index_date",
                      "index_date + 7 months"]),
    ),
    
    **loop_over_OPCS_codelists(list_dict[name_op], return_expectations ={"incidence": 1,"rate" : "uniform",}),

    
    first_surgery_date=patients.minimum_of(*[s + "_1_date_admitted" for s in list_dict.keys()]),
#"LeftHemicolectomy_date_admitted", "RightHemicolectomy_date_admitted","TotalColectomy_date_admitted", "RectalResection_date_admitted"
#    **loop_over_OPCS_codelists_discharge(list_dict,returning = "date_discharged", return_expectations ={"incidence": 1,"rate" : "uniform",}),
   
    first_surgery_discharge_date=patients.minimum_of(*[s + "_1_date_discharged" for s in list_dict.keys()]),
# "LeftHemicolectomy_date_discharged", "RightHemicolectomy_date_discharged","TotalColectomy_date_discharged", "RectalResection_date_discharged"
   
   
    registered=patients.registered_as_of(
    "first_surgery_date"    ## Minimum should have prior registration for first surgery
    ),

    follow_up=patients.registered_as_of(
    "first_surgery_date + 90 days" ## Minimum should have follow up for first surgery
    ),

    dereg_date=patients.date_deregistered_from_all_supported_practices(
        on_or_after="index_date",
        date_format="YYYY-MM",
        return_expectations={
            "incidence": 0.05
        }
    ),
    age=patients.age_as_of(
        "first_surgery_date",
        return_expectations={
            "rate" : "universal",
            "int" : {"distribution" : "population_ages"}
        }
    ),

    region=patients.registered_practice_as_of(
        "first_surgery_date",
        returning="nuts1_region_name",
        return_expectations={
            "rate": "universal",
            "category": {
                "ratios": { 
                    "Other": 0.1,
                    "North East": 0.1,
                    "North West": 0.1,
                    "Yorkshire and the Humber": 0.1,
                    "East Midlands": 0.1,
                    "West Midlands": 0.1,
                    "East of England": 0.1,
                    "London": 0.1,
                    "South East": 0.1,
                    "South West": 0.1
                },
            },
        },
    ),

    dob=patients.date_of_birth(  
        "YYYY-MM",
        return_expectations={
            "date": {"earliest": "1920-01-01", "latest": "today"},
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
            "category": {"ratios": {"3000": 0.2, "7000": 0.2, "14000": 0.2, "20000": 0.2, "27000": 0.2}},
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

    **loop_over_OPCS_codelists_major_minor(list_dict[name_op],returning = "binary_flag", return_expectations ={"incidence": 0.3,"rate" : "uniform",}),

    **loop_over_OPCS_codelists_admission_info(list_dict[name_op],returning = "admission_method", return_expectations ={"category": {"ratios": {"11": 0.1, "21": 0.2, "22": 0.7}}, "incidence" : 1}),

    **loop_over_OPCS_codelists_admission_info(list_dict[name_op],returning = "primary_diagnosis", return_expectations ={ "category": {"ratios": {"C150": 0.2,"C180": 0.2, "C190": 0.2, "K570": 0.2, "K512": 0.2}}, "incidence" : 1,}),

    **loop_over_OPCS_codelists_admission_info(list_dict[name_op],returning = "days_in_critical_care", return_expectations ={"category": {"ratios": {"5": 0.1, "6": 0.2, "7": 0.7}}, "incidence" : 0.1}),
   
    **with_sub_procedures(subproc_list=sub_proc_dict,code_list_dict=list_dict[name_op], returning="binary_flag", return_expectations={"incidence": 0.1,"rate" : "uniform",}),

   
    #########################################
    # Pre operative risk factors  

    CCF_HES=patients.admitted_to_hospital(
        between=["index_date","index_date + 3 years"],
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

 #   **post_operative_COVID(list_dict, returning= "case_category", return_expectations = {"incidence" : 1,"category": {"ratios": {"": 0.3, "LFT_Only": 0.4, "PCR_Only": 0.2, "LFT_WithPCR": 0.1}},}),

    **post_operative_COVID(list_dict[name_op], returning= "date", return_expectations = {"incidence" : 1,"rate" : "uniform"}),

  #  **recent_COVID(list_dict, returning= "case_category", return_expectations = {"incidence" : 1,"category": {"ratios": {"": 0.3, "LFT_Only": 0.4, "PCR_Only": 0.2, "LFT_WithPCR": 0.1}},}),

    **recent_COVID(list_dict[name_op], returning= "date", return_expectations = {"incidence" : 0.1,"rate" : "uniform"}),

  #  **previous_COVID(list_dict, returning= "case_category", return_expectations = {"incidence" : 0.1,"category": {"ratios": {"": 0.3, "LFT_Only": 0.4, "PCR_Only": 0.2, "LFT_WithPCR": 0.1}},}),

    **previous_COVID(list_dict[name_op], returning= "date", return_expectations = {"incidence" : 0.1,"rate" : "uniform"}),

    **with_emergency_readmissions(list_dict[name_op], returning = "primary_diagnosis",  return_expectations={"rate": "uniform","category": {"ratios": {"K920": 0.5, "K921": 0.5}}}),

    **with_emergency_readmissions(list_dict[name_op], returning = "date_admitted",  return_expectations={"rate" : "uniform"}),

    **with_post_op_hospital_events(name= "VTE",event_list=VTE_HES_codes,code_list_dict=list_dict[name_op], returning="date_admitted", return_expectations={"incidence": 0.5,"rate" : "uniform",}),

    **with_post_op_GP_events(name= "VTE",event_list=VTE_Read_codes,code_list_dict=list_dict[name_op], returning="date", return_expectations={"incidence": 0.5,"rate" : "uniform",}),

    **with_post_op_GP_medications(name= "anticoagulation",event_list=any_anticoagulation,code_list_dict=list_dict[name_op], returning="date", return_expectations={"incidence": 0.5,"rate" : "uniform",}),

 #   **with_post_op_hospital_events(name = "Trauma", event_list = all_trauma_codes,code_list_dict=list_dict, returning="date_admitted", return_expectations={"incidence": 0.1,"rate" : "uniform",}),

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
  
  death_underlying_cause_ons=patients.died_from_any_cause(
      on_or_after="index_date",
      returning="underlying_cause_of_death",
      return_expectations={"category": {"ratios": {"U071":0.2, "C33":0.2, "I60":0.1, "F01":0.1 , "F02":0.05 , "I22":0.05 ,"C34":0.05, "I23":0.25}},},
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
