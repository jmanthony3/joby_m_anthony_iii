@echo off
type nul > temp
sed -n "s%%[[:blank:]]*version = %%%%p" "pyproject.toml" | sed -n "1p" > temp
for /f %%i in ('type "temp"') do (
    set version=%%i
)
del temp
set version=%version:~1,-1%

del src\joby_m_anthony_iii\numerical_methods.log
del tests\numerical_methods.log

@REM pip install --upgrade build twine
@REM python -m build
@REM twine check ./dist/joby_m_anthony_iii-%version%*


@REM python -m twine upload --repository pypi ./dist/joby_m_anthony_iii-%version%.tar.gz
@REM timeout 5 /NOBREAK
@REM pip install --upgrade joby_m_anthony_iii==%version%

echo "Compiled, uploaded, and updated to 'joby_m_anthony_iii-%version%'. Generating documentation..."

rmdir /Q /S docs
sphinx-apidoc -f -M -F -H="joby_m_anthony_iii" -A="Joby M. Anthony III" -V="%version%" -o ./docs ./src/joby_m_anthony_iii
sed -n "/# -- Extension configuration /{=}" "docs/conf.py" | sed -n "1p" > temp
for /f %%i in ('type "temp"') do (
    set /A line_number=%%i+1
)
del temp
sed -i "%line_number%s/.*/extensions = [\n    'sphinx.ext.napoleon'\n]\nnapoleon_google_docstring = False\n/" "docs/conf.py"
sphinx-build -b html docs/ docs/_build/html
xcopy docs\_build\html docs /E
del numerical_methods.log