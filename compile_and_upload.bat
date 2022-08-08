@echo off

rem ------------------------------------------------------------
rem     METADATA
rem ------------------------------------------------------------
rem get project version number
type nul > temp
sed -n "s%%[[:blank:]]*version = %%%%p" "pyproject.toml" | sed -n "1p" > temp
for /f %%i in ('type "temp"') do (
    set version=%%i
)
del temp
rem remove double quotes from either end
set version=%version:~1,-1%

del src\joby_m_anthony_iii\numerical_methods.log
del tests\numerical_methods.log

@REM rem ------------------------------------------------------------
@REM rem     Build Python wheel
@REM rem ------------------------------------------------------------
@REM rem upload to PyPi and upgrade local
@REM pip install --upgrade build twine
@REM python -m build
@REM twine check ./dist/joby_m_anthony_iii-%version%*
@REM python -m twine upload --repository pypi ./dist/joby_m_anthony_iii-%version%.tar.gz

pip install --upgrade joby_m_anthony_iii==%version%

echo "Compiled, uploaded, and updated to 'joby_m_anthony_iii-%version%'. Generating documentation..."

rem ------------------------------------------------------------
rem     Generate/Update API Documentation
rem ------------------------------------------------------------
rmdir /Q /S docs

sphinx-apidoc -f -M -F -H="joby_m_anthony_iii" -A="Joby M. Anthony III" -V="%version%" -o ./docs ./src/joby_m_anthony_iii

sed -i "s/html_theme = .*/html_theme = 'sphinx_rtd_theme'/" "docs/conf.py"
@REM sed -i "s/html_theme = .*/html_theme = 'python_docs_theme'/" "docs/conf.py"

sed -n "/# -- Extension configuration /{=}" "docs/conf.py" | sed -n "1p" > temp
for /f %%i in ('type "temp"') do (
    set /A line_number=%%i+1
)
del temp
sed -i "%line_number%s/.*/extensions = [\n    'sphinx.ext.napoleon'\n]\nnapoleon_google_docstring = False\n/" "docs/conf.py"

sed -i "s@BUILDDIR      = _build@BUILDDIR      = ../../joby_m_anthony_iii-docs\nPDFBUILDDIR   = /tmp\nPDF           = ../manual.pdf\n@" "docs/Makefile"

(
echo.
echo latexpdf:
echo 	$^(SPHINXBUILD^) -b latex $^(ALLSPHINXOPTS^) $^(PDFBUILDDIR^)^/latex
echo 	#                                          ^^^^^^
echo 	@echo "Running LaTeX files through pdflatex..."
echo 	^make -C $^(PDFBUILDDIR^)^/latex all-pdf
echo 	#         ^^^^^^
echo 	^copy $^(PDFBUILDDIR^)^/latex^/*.pdf $^(PDF^)
echo 	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
echo 	@echo "pdflatex finished; see $(PDF)"
)>>"docs/Makefile"

echo .. include:: ../../README.md> "docs/includeme.rst"

sed -n "/:caption: Contents:/{=}" "docs/index.rst" | sed -n "1p" > temp
for /f %%i in ('type "temp"') do (
    set /A line_number=%%i+1
)
del temp
sed -i "%line_number%s/.*/   includeme\n/" "docs/index.rst"

@REM sphinx-build -b html docs/ docs/_build/html
@REM @REM sphinx-build -b html docs/ "..\\jmanthony3.github.io\\joby_m_anthony_iii\html"
@REM xcopy "docs\_build\doctrees" "..\joby_m_anthony_iii-docs\doctrees\" /E /Y && xcopy "docs\_build\html" "..\joby_m_anthony_iii-docs\html\" /E /Y
@REM cd docs && make latexpdf && cd .. && copy "docs\_build\latex\*.pdf" "manual.pdf" /Y
@REM del numerical_methods.log
cd docs && make html && cd .. && xcopy "docs\_build\doctrees" "..\joby_m_anthony_iii-docs\doctrees\" /E /Y && xcopy "docs\_build\html" "..\joby_m_anthony_iii-docs\html\" /E /Y && cd docs && make latexpdf && cd .. && copy "docs\_build\latex\*.pdf" "manual.pdf" /Y && cd "..\joby_m_anthony_iii-docs\html" && git add . && git commit -m "rebuilt docs" && git push origin gh-pages && cd "..\..\joby_m_anthony_iii"