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
sed -n "/BUILDDIR      = _build/{=}" "docs/Makefile" | sed -n "1p" > temp
for /f %%i in ('type "temp"') do (
    set /A line_number=%%i
)
del temp
sed -i "%line_number%s@.*@BUILDDIR      = ../../joby_m_anthony_iii-docs\nPDFBUILDDIR   = /tmp\nPDF           = ../manual.pdf\n@" "docs/Makefile"

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
@REM sphinx-build -b html docs/ "..\\jmanthony3.github.io\\joby_m_anthony_iii\html"
@REM del numerical_methods.log
cd docs && make html && cd .. && xcopy "docs\_build\doctrees" "..\joby_m_anthony_iii-docs\doctrees\" /E /Y && xcopy "docs\_build\html" "..\joby_m_anthony_iii-docs\html\" /E /Y && cd docs && make latexpdf && cd .. && copy "docs\_build\latex\*.pdf" "manual.pdf" /Y