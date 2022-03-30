from cohortextractor import StudyDefinition, patients, codelist, codelist_from_csv  
from codelists import *


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
       registered
        AND (follow_up OR died)
        AND surgery_date
        AND (age >=18 AND age <= 110)
        AND (sex = "M" OR sex = "F")
       """,
       registered=patients.registered_with_one_practice_between(
       "surgery_date - 365 days", "surgery_date"
       ),

       follow_up=patients.registered_as_of(
       "surgery_date + 90 days"
       )
   ),


   # Define earliest colorectal surgery date only for this initial analysis

   surgery_date=patients.admitted_to_hospital(
        returning="date_admitted",
        with_these_procedures=any_colorectal_resection,
        between=["index_date","index_date + 2 years"],
        find_first_match_in_period=True,
        date_format="YYYY-MM-DD",          
        return_expectations={
            "rate" : "uniform",
            "incidence": 1,
        }
    ),

    surgery_discharge_date=patients.admitted_to_hospital(
        with_these_procedures=any_colorectal_resection,
        find_first_match_in_period=True,
        returning="date_discharged",
        between=["surgery_date","surgery_date + 2 years"],
        date_format="YYYY-MM-DD",
        return_expectations={
           "rate" : "uniform",
            "incidence": 1,
        }
    ),

    surgery_admimeth=patients.admitted_to_hospital(
        find_first_match_in_period=True,
        between=["surgery_date","surgery_date"],
        returning="admission_method",
        return_expectations={"category": {"ratios": {"11": 0.1, "21": 0.2, "22": 0.7}}, "incidence" : 1
        },
        with_these_procedures=any_colorectal_resection
    ),

    surgery_cancer=patients.admitted_to_hospital(
        find_first_match_in_period=True,
        between=["surgery_date","surgery_date"],
        returning="primary_diagnosis",
        return_expectations={            
            "category": {"ratios": {"C180": 0.5, "C190": 0.5}}, "incidence" : 0.5,
        },
        with_these_procedures=any_colorectal_resection
    ),

    first_emergency_readmission_diagnosis=patients.admitted_to_hospital(
        with_admission_method = ['21', '2A', '22', '23', '24', '25', '2D'],
        returning="primary_diagnosis",
        between=["surgery_discharge_date", "surgery_discharge_date + 1 year"],
        find_first_match_in_period=True,
        return_expectations={
            "rate": "universal",
            "category": {"ratios": {"K920": 0.5, "K921": 0.5}}
        },    
    ),
    
    index_primary_diagnosis=patients.admitted_to_hospital(
        with_these_procedures=any_colorectal_resection,
        find_first_match_in_period=True,
        returning="primary_diagnosis",
        between=["surgery_date","surgery_discharge_date"],
        return_expectations={
            "category": {"ratios": {"C180": 0.5, "C190": 0.5}},
        },     
    ),

    index_ICU_days=patients.admitted_to_hospital(
        with_these_procedures=any_colorectal_resection,
        find_first_match_in_period=True,
        returning="days_in_critical_care",
        between=["surgery_date","surgery_discharge_date"],
        return_expectations = {"category": {"ratios": {"5": 0.1, "6": 0.2, "7": 0.7}}, "incidence" : 0.1}, #{"int" : {"distribution": "poisson", "mean": 5}, "incidence" : 0.6},
    ),

    region=patients.registered_practice_as_of(
        "surgery_date",
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


    age=patients.age_as_of(
        "surgery_date",
        return_expectations={
            "rate" : "universal",
            "int" : {"distribution" : "population_ages"}
        }   
    ),

    sex=patients.sex(
        return_expectations={
            "rate": "universal",
            "category": {"ratios": {"M": 0.49, "F": 0.51}},
        }
    ),

    bmi=patients.most_recent_bmi(
        between=["index_date", "surgery_date"],
        minimum_age_at_measurement=18,
        include_measurement_date=True,
        date_format="YYYY-MM",
        return_expectations={
            "float": {"distribution": "normal", "mean": 28, "stddev": 8},
            "incidence": 0.80,
        }
    ),
    
    imd=patients.address_as_of(
        'surgery_date',
        returning="index_of_multiple_deprivation",
        round_to_nearest=100,
        return_expectations={
            "rate": "universal",
            "category": {"ratios": {"100": 0.1, "200": 0.2, "300": 0.7}},
        },
    ),


    **with_these_vaccination_date_X(
        name = "covid_vaccine_dates", 
 #       returning="date",
        index_date = "index_date",
        n = 3, 
        return_expectations = {
            "rate": "uniform",
            "incidence": 0.85,
        },
    ),

#     #########################################
#     # Procedure details
#     #########################################

    right_hemicolectomy=patients.admitted_to_hospital(
        find_first_match_in_period=True,
        between=["index_date","index_date + 2 years"],
        returning="date_admitted",
        date_format="YYYY-MM-DD",
        return_expectations={
            "incidence": 0.1,
        },
        with_these_procedures=right_hemicolectomy_codes
    ),

    right_hemicolectomy_admimeth=patients.admitted_to_hospital(
        find_first_match_in_period=True,
        between=["index_date","index_date + 2 years"],
        returning="admission_method",
        return_expectations={"category": {"ratios": {"11": 0.1, "21": 0.2, "22": 0.7}}, "incidence" : 1
        },
        with_these_procedures=right_hemicolectomy_codes
    ),

    left_hemicolectomy=patients.admitted_to_hospital(
        find_first_match_in_period=True,
        between=["index_date","index_date + 2 years"],
        returning="date_admitted",
        date_format="YYYY-MM-DD",
        return_expectations={
            "incidence": 0.1,
        },
        with_these_procedures=left_hemicolectomy_codes
    ),

    left_hemicolectomy_admimeth=patients.admitted_to_hospital(
        find_first_match_in_period=True,
        between=["index_date","index_date + 2 years"],
        returning="admission_method",
        return_expectations={"category": {"ratios": {"11": 0.1, "21": 0.2, "22": 0.7}}, "incidence" : 1
        },
        with_these_procedures=left_hemicolectomy_codes
    ),

    total_colectomy=patients.admitted_to_hospital(
        find_first_match_in_period=True,
        between=["index_date","index_date + 2 years"],
        returning="date_admitted",
        date_format="YYYY-MM-DD",
        return_expectations={
            "incidence": 0.1,
        },
        with_these_procedures=total_colectomy_codes
    ),

    total_colectomy_admimeth=patients.admitted_to_hospital(
        find_first_match_in_period=True,
        between=["index_date","index_date + 2 years"],
        returning="admission_method",
        return_expectations={"category": {"ratios": {"11": 0.1, "21": 0.2, "22": 0.7}}, "incidence" : 1
        },
        with_these_procedures=total_colectomy_codes
    ),

    rectal_resection=patients.admitted_to_hospital(
        find_first_match_in_period=True,
        between=["index_date","index_date + 2 years"],
        returning="date_admitted",
        date_format="YYYY-MM-DD",
        return_expectations={
            "incidence": 0.1,
        },
        with_these_procedures=rectal_resection_codes
    ),    

    rectal_resection_admimeth=patients.admitted_to_hospital(
        find_first_match_in_period=True,
        between=["index_date","index_date + 2 years"],
        returning="admission_method",
        return_expectations={
            "category": {"ratios": {"11": 0.1, "21": 0.2, "22": 0.7}}, 
            "incidence" : 1,
        },
        with_these_procedures=rectal_resection_codes
    ),

    #########################################
    # Pre operative risk factors
    ##########################################

 
     MI_HES=patients.admitted_to_hospital(
        between=["index_date","surgery_date"],
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
        between=["index_date","surgery_date"],
        returning="date_admitted",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",
            "incidence": 0.1,
        },
        with_these_diagnoses=CCF_HES_codes
    ),

    MI_GP=patients.with_these_clinical_events(
        MI_Read_codes,
        between=["index_date","surgery_date"],
        returning="date",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",            
            "incidence": 0.1,
        },
    ),

    CCF_GP=patients.with_these_clinical_events(
        CCF_Read_codes,
        between=["index_date","surgery_date"],
        returning="date",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",            
            "incidence": 0.1,
        },
    ),

    Stroke_HES=patients.admitted_to_hospital(
        between=["index_date","surgery_date"],
        returning="date_admitted",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",
            "incidence": 0.1,
        },
        with_these_diagnoses=Stroke_HES_codes
    ),

    Stroke_GP=patients.with_these_clinical_events(
        Stroke_Read_codes,
        between=["index_date","surgery_date"],
        returning="date",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",            
            "incidence": 0.1,
        },
    ),

    PVD_GP=patients.with_these_clinical_events(
        PVD_Read_codes,
        between=["index_date","surgery_date"],
        returning="date",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",            
            "incidence": 0.1,
        },
    ),

    Dementia_GP=patients.with_these_clinical_events(
        Dementia_Read_codes,
        between=["index_date","surgery_date"],
        returning="date",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",            
            "incidence": 0.1,
        },
    ),

    Chronic_Respiratory_GP=patients.with_these_clinical_events(
        Chronic_Respiratory_Read_codes,
        between=["index_date","surgery_date"],
        returning="date",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",            
            "incidence": 0.1,
        },
    ),


    RA_SLE_Psoriasis_GP=patients.with_these_clinical_events(
        RA_SLE_Psoriasis_Read_codes,
        between=["index_date","surgery_date"],
        returning="date",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",            
            "incidence": 0.1,
        },
    ),

    Ulcer_or_bleed_GP=patients.with_these_clinical_events(
        Ulcer_or_bleed_Read_codes,
        between=["index_date","surgery_date"],
        returning="date",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",            
            "incidence": 0.1,
        },
    ),

    all_liver_GP=patients.with_these_clinical_events(
        all_liver_disease_Read_codes,
        between=["index_date","surgery_date"],
        returning="date",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",            
            "incidence": 0.1,
        },
    ),

    cirrhosis_GP=patients.with_these_clinical_events(
        cirrhosis_snomed_codes,
        between=["index_date","surgery_date"],
        returning="date",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",            
            "incidence": 0.1,
        },
    ),
    all_diabetes_GP=patients.with_these_clinical_events(
        all_diabetes_Read_codes,
        between=["index_date","surgery_date"],
        returning="date",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",            
            "incidence": 0.1,
        },
    ),

    diabetic_complication_GP=patients.with_these_clinical_events(
        diabetic_complication_snomed_codes,
        between=["index_date","surgery_date"],
        returning="date",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",            
            "incidence": 0.1,
        },
    ),
    
    other_neuro_GP=patients.with_these_clinical_events(
        other_neuro_Read_codes,
        between=["index_date","surgery_date"],
        returning="date",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",            
            "incidence": 0.1,
        },
    ),       

    Renal_GP=patients.with_these_clinical_events(
        Chronic_kidney_snomed_codes,
        between=["index_date","surgery_date"],
        returning="date",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",            
            "incidence": 0.1,
        },
    ),        

    CKD_3_5_GP=patients.with_these_clinical_events(
        CKD_Read_codes,
        between=["index_date","surgery_date"],
        returning="date",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",            
            "incidence": 0.1,
        },
    ),        

    non_haem_cancer_GP=patients.with_these_clinical_events(
        non_haem_cancer,
        between=["index_date","surgery_date"],
        returning="date",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",            
            "incidence": 0.1,
        },
    ),   

    haem_cancer_GP=patients.with_these_clinical_events(
        haem_cancer_Read_codes,
        between=["index_date","surgery_date"],
        returning="date",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",            
            "incidence": 0.1,
        },
    ),   

    metastatic_cancer_GP=patients.with_these_clinical_events(
        metastatic_cancer_snomed_codes,
        between=["index_date","surgery_date"],
        returning="date",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",            
            "incidence": 0.1,
        },
    ),   

    HIV_GP=patients.with_these_clinical_events(
        HIV_Read_codes,
        between=["index_date","surgery_date"],
        returning="date",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
           "rate" : "uniform",            
            "incidence": 0.1,
        },
    ),   


    ####################################
    # Follow up events
    ######################################

    SARS_CoV_2_test_type = patients.with_test_result_in_sgss(
        pathogen = "SARS-CoV-2",
        returning = "case_category",
        test_result = "positive",
        between=["surgery_discharge_date", "surgery_discharge_date + 90 days"],
        return_expectations = {
            "incidence" : 1,
            "category": {"ratios": {"": 0.3, "LFT_Only": 0.4, "PCR_Only": 0.2, "LFT_WithPCR": 0.1}},
             }
    ),

    SARS_CoV_2_test_date = patients.with_test_result_in_sgss(
        pathogen = "SARS-CoV-2",
        between=["surgery_discharge_date", "surgery_discharge_date + 90 days"],
        returning = "date",
    	date_format = "YYYY-MM-DD",
        test_result = "positive",
        return_expectations = {
            "incidence" : 1,
            "rate" : "uniform"
             }
    ),
    first_emergency_readmission_date=patients.admitted_to_hospital(
        with_admission_method = ['21', '2A', '22', '23', '24', '25', '2D'],
        returning="date_admitted",
        between=["surgery_discharge_date", "surgery_discharge_date + 90 days"],
        find_first_match_in_period=True,
        date_format="YYYY-MM-DD",
        return_expectations={
            "rate" : "uniform"
        }
    ),

    VTE_HES=patients.admitted_to_hospital(
        between=['surgery_date', 'surgery_discharge_date + 90 days'],
        returning="date_admitted",
        date_format="YYYY-MM-DD",
        return_expectations={
            "incidence": 0.1,
            "rate" : "uniform",            
        },
        with_these_diagnoses=VTE_HES_codes
    ),

    VTE_GP=patients.with_these_clinical_events(
        VTE_Read_codes,
        between=['surgery_date', 'surgery_discharge_date + 90 days'],
        returning="date",
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
            "rate" : "uniform",
            "incidence": 0.1,
        },
    ),

    Anticoagulant_prescription=patients.with_these_medications(
        any_anticoagulation,
        between=['surgery_date', 'surgery_discharge_date + 90 days'],
        returning="date",
          date_format="YYYY-MM-DD",
        find_first_match_in_period=True,      
        return_expectations={
            "incidence": 0.1,
            "rate" : "uniform",
        },
    ),

    died=patients.died_from_any_cause(
        between=["surgery_date","surgery_discharge_date + 90 days"],
        returning="binary_flag",
        return_expectations={
            "rate" : "uniform",
            "incidence": 0.1,
        },
    ),
	date_death_ons = patients.died_from_any_cause(
        between=["surgery_date","surgery_discharge_date + 90 days"],
		returning = "date_of_death",
		date_format = "YYYY-MM-DD",
		return_expectations = {
			"rate": "exponential_increase"
		},
	),

	date_death_cpns = patients.with_death_recorded_in_cpns(
        between=["surgery_date","surgery_discharge_date + 90 days"],
		returning = "date_of_death",
		date_format = "YYYY-MM-DD",
		return_expectations = {
			"rate": "exponential_increase"
		},
	),

)
