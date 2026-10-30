@echo off

REM Generate the table of parameter combinations
call popy run -o trials.pyml
    
REM Replicate the template with actual values instead of placeholders
call popy replicate -t trials.pyml_output\synthetic_data.csv mytut.pyml.template
    
REM Execute every replicate in sequence
call popy run -o mytut.*.pyml
