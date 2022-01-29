/*
MIT License

Copyright (c) 2022 Jakob Roar Bentzon

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <stdio.h>

#ifndef UCLIB_H
#define UCLIB_H
#if DOUBLE_PRECISION
typedef double Real;
#else
typedef float Real;
#endif
typedef double CoordReal;

#ifdef __cplusplus
extern "C"
{
#endif
#if defined(WIN32) || defined(_WINDOWS) || defined(_WINNT)
#define USERFUNCTION_EXPORT __declspec(dllexport)
#define USERFUNCTION_IMPORT __declspec(dllimport)
#else
#define USERFUNCTION_EXPORT
#define USERFUNCTION_IMPORT
#endif

    extern void USERFUNCTION_IMPORT ucarg(void *, char *, char *, int);
    extern void USERFUNCTION_IMPORT ucfunc(void *, char *, char *);
    extern void USERFUNCTION_IMPORT ucfunction(void *, char *, char *, int, ...);

    // extern "C"
    // {
         void USERFUNCTION_EXPORT helloWorld()
        {
            printf("Hello World!");
        }
    // }

    void USERFUNCTION_EXPORT uclib();
#ifdef __cplusplus
}
#endif
#endif
