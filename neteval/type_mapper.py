import pandas as pd
import numpy as np

def assign_type():
    pass

def import_type_classifications():
    pass

def clean_type_data():
    pass

def process_type_data():
    pass

def type_mapper_wrapper(no_type):
    if no_type:
        assign_type()
    else:
        import_type_classifications()
        clean_type_data()
        process_type_data()
    pass
