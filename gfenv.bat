@ECHO OFF
if "%1" == "R" GOTO RELEASE

SET PYTHONPATH=%CD%;%CD%\..\_vs2015\Debug_x86\_gfrd\
GOTO END

:RELEASE
SET PYTHONPATH=%CD%;%CD%\..\_vs2015\Release_x86\_gfrd\

:END