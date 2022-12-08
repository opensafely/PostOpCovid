from cohortextractor import StudyDefinition, patients, codelist, codelist_from_csv  
from codelists import *

def with_sub_procedures(majorminor_list,code_list_dict, returning, return_expectations):
    def with_these_procedures(subkey,majorminorcodes,key,admission_date, admission90_date,returning,return_expectations,i):
        return {
            f"{key}_{i}_{subkey}_HES_{returning}": (
                patients.admitted_to_hospital(
                    with_these_procedures = majorminorcodes,
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
        for majorminorkey,majorminorcodes in majorminor_list.items():
          variables.update(with_these_procedures(majorminorkey,majorminorcodes,key,f"{key}_{i}_date_admitted",f"{key}_{i}_date_admitted",returning,return_expectations,i))
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
    
    **with_sub_procedures(majorminor_list=majormino_dict,code_list_dict=list_dict, returning="binary_flag", return_expectations={"incidence": 0.1,"rate" : "uniform",}),

)
