@echo off

del test.exe

g++ -o test.exe  lib/main.cpp

start cmd.exe /k "test.exe & echo."
