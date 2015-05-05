@ECHO OFF
if "%1" == "R" GOTO RELEASE

SET PYTHONPATH=%CD%;%CD%\..\_vs2013\Debug_x86\_gfrd\
GOTO END

:RELEASE
SET PYTHONPATH=%CD%;%CD%\..\_vs2013\Release_x86\_gfrd\

:END