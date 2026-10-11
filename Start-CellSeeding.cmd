@echo off
rem Double-click to open the cell seeding calculator (GUI).
cd /d "%~dp0"
where pyw >nul 2>nul
if %errorlevel%==0 (
    start "" pyw -3 cell_seeding_calc.py --gui
) else (
    start "" pythonw cell_seeding_calc.py --gui
)
