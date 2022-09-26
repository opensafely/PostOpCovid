from cohortextractor import StudyDefinition, patients, codelist, codelist_from_csv, params
from codelists import *

n_op =  params["n_op"]
name_op = params["name_op"]
op_date_admit_name = name_op+"_"+n_op+"_date_admitted"
op_date_discharge_name = name_op+"_"+n_op+"_date_discharged"

def post_operative_COVID(key, i, returning, return_expectations):
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
    variables = with_these_procedures(key,"date_admitted - 7 days","date_discharged + 90 days",returning,return_expectations,i)
    return variables


def recent_COVID(key, i, returning, return_expectations):
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
    variables = with_these_procedures(key,"date_admitted - 42 days","date_discharged - 8 days",returning,return_expectations,i)
    return variables


def previous_COVID(key, i, returning, return_expectations):
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
    variables = with_these_procedures(key,"date_admitted - 43 days",returning,return_expectations,i)
    return variables


def with_emergency_readmissions(key, i, returning, return_expectations):
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
    variables = with_these_procedures(key,"date_discharged","date_discharged + 90 days",returning,return_expectations,i)
    return variables

def with_post_op_hospital_events(name,event_list,key,i, returning, return_expectations):
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
    variables = with_these_procedures(name,event_list,key,"date_admitted","date_discharged + 90 days",returning,return_expectations,i)
    return variables

def with_post_op_GP_events(name,event_list,key,i, returning, return_expectations):
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
    variables = with_these_procedures(name,event_list,key,"date_admitted","date_discharged + 90 days",returning,return_expectations,i)
    return variables

def with_post_op_GP_medications(name,event_list,key,i, returning, return_expectations):
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
    variables = with_these_procedures(name,event_list,key,"date_admitted","date_discharged + 90 days",returning,return_expectations,i)
    return variables


def admitdates(key,i):
    return {f"{key}_{i}_date_admitted":patients.with_value_from_file("output/input.csv",
          returning=op_date_admit_name,
          returning_type="date"),
          f"{key}_{i}_date_discharged":patients.with_value_from_file("output/input.csv",
          returning=op_date_discharge_name,
          returning_type="date")},


study = StudyDefinition(
    index_date = "2021-11-01",

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
        """,
        has_surgery=patients.which_exist_in_file("output/input.csv"),
    ),
    
   # **admitdates(key = name_op, i = n_op),
   date_admitted = patients.with_value_from_file("output/input.csv",
          returning=op_date_admit_name,
          returning_type="date"),
   
   date_discharged = patients.with_value_from_file("output/input.csv",
          returning=op_date_discharge_name,
          returning_type="date"),          

    age=patients.age_as_of(
         "date_admitted" ,
        return_expectations={
            "rate" : "universal",
            "int" : {"distribution" : "population_ages"}
        }
    ),

    region=patients.registered_practice_as_of(
         "date_admitted",
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
        between=["index_date",  "date_admitted"],
        minimum_age_at_measurement=18,
        include_measurement_date=True,
        date_format="YYYY-MM",
        return_expectations={
            "float": {"distribution": "normal", "mean": 28, "stddev": 8},
            "incidence": 0.80,
        }
    ),
    
    imd=patients.address_as_of(
         "date_admitted",
        returning="index_of_multiple_deprivation",
        round_to_nearest=100,
        return_expectations={
            "rate": "universal",
            "category": {"ratios": {"3000": 0.2, "7000": 0.2, "14000": 0.2, "20000": 0.2, "27000": 0.2}},
        },
    ),

    ####################################
    # Follow up events
    ######################################

 #   **post_operative_COVID(key = name_op, i = n_op, returning= "case_category", return_expectations = {"incidence" : 1,"category": {"ratios": {"": 0.3, "LFT_Only": 0.4, "PCR_Only": 0.2, "LFT_WithPCR": 0.1}},}),

    **post_operative_COVID(key = name_op, i = n_op, returning= "date", return_expectations = {"incidence" : 1,"rate" : "uniform"}),

  #  **recent_COVID(key = name_op, i = n_op, returning= "case_category", return_expectations = {"incidence" : 1,"category": {"ratios": {"": 0.3, "LFT_Only": 0.4, "PCR_Only": 0.2, "LFT_WithPCR": 0.1}},}),

    **recent_COVID(key = name_op, i = n_op, returning= "date", return_expectations = {"incidence" : 0.1,"rate" : "uniform"}),

  #  **previous_COVID(key = name_op, i = n_op, returning= "case_category", return_expectations = {"incidence" : 0.1,"category": {"ratios": {"": 0.3, "LFT_Only": 0.4, "PCR_Only": 0.2, "LFT_WithPCR": 0.1}},}),

    **previous_COVID(key = name_op, i = n_op, returning= "date", return_expectations = {"incidence" : 0.1,"rate" : "uniform"}),


    **with_emergency_readmissions(key = name_op, i = n_op, returning = "primary_diagnosis",  return_expectations={"rate": "uniform","category": {"ratios": {"K920": 0.5, "K921": 0.5}}}),

    **with_emergency_readmissions(key = name_op, i = n_op, returning = "date_admitted",  return_expectations={"rate" : "uniform"}),

    **with_post_op_hospital_events(name= "VTE",event_list=VTE_HES_codes,key = name_op, i = n_op, returning="date_admitted", return_expectations={"incidence": 0.5,"rate" : "uniform",}),

    **with_post_op_GP_events(name= "VTE",event_list=VTE_Read_codes,key = name_op, i = n_op, returning="date", return_expectations={"incidence": 0.5,"rate" : "uniform",}),

    **with_post_op_GP_medications(name= "anticoagulation",event_list=any_anticoagulation,key = name_op, i = n_op, returning="date", return_expectations={"incidence": 0.5,"rate" : "uniform",}),

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

	date_death_cpns = patients.with_death_recorded_in_cpns(
        between=["index_date","index_date + 3 years"],
		returning = "date_of_death",
		date_format = "YYYY-MM-DD",
		return_expectations = {
			"rate": "exponential_increase"
		},
	),

)
