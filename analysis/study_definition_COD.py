from cohortextractor import StudyDefinition, patients, codelist, codelist_from_csv  
from codelists import *


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
        AND died
        """,
       
        in_surgery_cohort=patients.which_exist_in_file("output/input.csv"),
        
        died=patients.died_from_any_cause(
          on_or_after="index_date",
          returning="binary_flag",
          return_expectations={
            "rate" : "uniform",
            "incidence": 0.1,
        },
      ),
    ),
    
    

  	date_death_ons = patients.died_from_any_cause(
      on_or_after="index_date",
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
)
