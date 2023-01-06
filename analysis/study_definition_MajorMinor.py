from cohortextractor import StudyDefinition, patients, codelist, codelist_from_csv, params
from codelists import *

i =  params["n_op"]
key = params["name_op"]
op_date_admit_name = key+"_"+i+"_date_admitted"
op_date_discharge_name = key+"_"+i+"_date_discharged"


def with_these_procedures(major_codes,key, i,date_admitted, date_discharged,returning,return_expectations):
  return {
    f"{key}_{i}_Major_HES_{returning}": (
      patients.admitted_to_hospital(
      with_these_procedures = major_codes,
      find_first_match_in_period=True,
      returning = returning,
      date_format="YYYY-MM-DD",
      between=[date_admitted, date_discharged],
      return_expectations = return_expectations,
      )
    )
  }



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
    
  **with_these_procedures(major_codes,key,i, 'date_admitted','date_discharged',returning="binary_flag",return_expectations={"incidence": 0.1,"rate" : "uniform",})
)
