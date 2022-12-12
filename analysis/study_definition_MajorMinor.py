from cohortextractor import StudyDefinition, patients, codelist, codelist_from_csv  
from codelists import *

n_op =  params["n_op"]
name_op = params["name_op"]
op_date_admit_name = name_op+"_"+n_op+"_date_admitted"
op_date_discharge_name = name_op+"_"+n_op+"_date_discharged"


def with_sub_procedures(major_codes, key, i,returning, return_expectations):
    def with_these_procedures(major_codes,key, i,admission_date, discharge_date,returning,return_expectations):
        return {
            f"{key}_{i}_Major_HES_{returning}": (
                patients.admitted_to_hospital(
                    with_these_procedures = major_codes,
                    find_first_match_in_period=True,
                    returning = returning,
                    date_format="YYYY-MM-DD",
                    between=[admission_date, discharge_date],
                    return_expectations = return_expectations,
                )
            )
        }
    variables = {}
        variables.update(with_these_procedures(major_codes,key,i,admission_date = f"{key}_{i}_date_admitted",discharge_date = f"{key}_{i}_date_discharged",returning,return_expectations))
    return variables


study = StudyDefinition(
    index_date = "2020-02-01",

    default_expectations={
        "date": {"earliest": "index_date" , "latest": "index_date + 3 years"},
        "rate": "uniform",
        "incidence": 1,
    },

    population=patients.satisfying(
        """
        in_surgery_cohort
        """,
       
        in_surgery_cohort=patients.which_exist_in_file("output/input.csv"),
    ),
    
   date_admitted = patients.with_value_from_file("output/input.csv",
          returning=op_date_admit_name,
          returning_type="date"),
   
   date_discharged = patients.with_value_from_file("output/input.csv",
          returning=op_date_discharge_name,
          returning_type="date"),
    
    **with_sub_procedures(major_codes=major_codes,key = name_op, i = n_op, returning="binary_flag", return_expectations={"incidence": 0.1,"rate" : "uniform",}),

)
