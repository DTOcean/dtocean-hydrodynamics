@echo off
(
pyuic4 ../designer/load_new.ui > ../load_new.py
pyuic4 ../designer/mdi_layout.ui > ../mdi_layout.py
pyuic4 ../designer/new_selection.ui > ../new_selection.py
pyuic4 ../designer/performance_fit.ui > ../performance_fit.py
pyuic4 ../designer/plotter.ui > ../plotter.py
pyuic4 ../designer/read_db_form.ui > ../read_db_form.py
pyuic4 ../designer/read_nemoh_form.ui > ../read_nemoh_form.py
pyuic4 ../designer/read_wamit_form.ui > ../read_wamit_form.py
pyuic4 ../designer/run_nemoh_form.ui > ../run_nemoh_form.py
)
pause