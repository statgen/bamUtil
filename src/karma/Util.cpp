/*
 * Copyright (c) 2009 Regents of the University of Michigan
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include <assert.h>
#include <math.h>
#include <string.h>
#include "IntHash.h"
#include "Util.h"

String getBaseName(String pathName)
{
    if (pathName.Length()==0) return "";
    String baseName;
    int baseNameStart = pathName.FindLastChar('/');
    if (baseNameStart == -1)
    {
        baseNameStart = pathName.FindLastChar('\\');
        if (baseNameStart == -1)
        {
            baseNameStart = 0;
        }
    }
    else baseNameStart++;

    int suffixStart = pathName.FindLastChar('.');

    if (suffixStart == -1) suffixStart = pathName.Length() - 1;
    else suffixStart--;

    baseName = pathName.Mid(baseNameStart, suffixStart);
    return baseName;
}

bool IsPrime(int n)
{
    assert(n > 0);
    IntHash prime;
    IntHash notprime;

    prime.Add(1);
    prime.Add(2);
    prime.Add(3);
    prime.Add(5);
    prime.Add(7);

    if (prime.Find(n) != -1)
        return true;

    if (n % 2 == 0 || n % 3 ==0 || notprime.Find(n) != -1)
        return false;

    int sqrtn = (int) sqrt(n * 1.0);

    int x = 5;

    while (x <= sqrtn)
    {
        if (n % x == 0)
        {
            notprime.Add(n);
            return false;
        }
        x += 2;
        if (n % x == 0)
        {
            notprime.Add(n);
            return false;
        }
        x += 4;
    }

    prime.Add(n);
    return true;

}

int ClosestPrimeBelow(int n)
{
    if (IsPrime(n))
        return n;

    n -= 1; // ensure not returning n+1 (e.g, n=4);
    n |= 1;
    while (!IsPrime(n))
        n -= 2;

    return n;
}

int ClosestPrimeAbove(int n)
{
    if (IsPrime(n))
        return n;

    n |= 1;
    while (!IsPrime(n))
        n += 2;

    return n;
}

//
// XXX - refactor this - it belongs in either GenomeSequence or WordIndex
//
unsigned int MapBaseToInteger(char base)
{
    switch (base)
    {
    case 'A' :
    case 'a' :
        return 0;
    case 'T' :
    case 't' :
        return 1;
    case 'C' :
    case 'c' :
        return 2;
    case 'G' :
    case 'g' :
        return 3;
    case 'N' :
    case 'n' :
        return 4;
    default:
        return 5;
    }
}

//
// XXX - refactor this - it belongs in either GenomeSequence or WordIndex
//
char MapIntegerToBase(unsigned int i)
{
    switch (i)
    {
    case 0 :
        return 'A';
    case 1:
        return 'T';
    case 2:
        return 'C';
    case 3:
        return 'G';
    case 4:
        return 'N';
    default:
        return 'M';
    }
}

int power(int x, int y)
{
    assert(y >= 0);

    int result = 1;
    for (int i = 0; i < y; i++)
        result *= x;

    return result;
}

bool signalPoll::userQuit;

signalPoll::signalPoll()
{
    userQuit = false;
    pollingForQuit = false;
}

signalPoll::~signalPoll()
{
    disableQuit();
}


static void userQuitHandler(int i)
{
    signalPoll::userQuit = true;
}

void signalPoll::enableQuit()
{
    struct sigaction sa;
    int rc;
    memset(&sa, '\0', sizeof(sa));
    sa.sa_handler = userQuitHandler;
#if 0
    sa.sa_mask = 0;
    sa.sa_flags = 0;
    sa.sa_restorer = NULL;
#endif
    pollingForQuit = true;
    rc = sigaction(SIGINT, &sa, &oldQuitAction);
    if (rc)
    {
        perror("failed to set signal");
    }
}

void signalPoll::disableQuit()
{
    sigaction(SIGQUIT, &oldQuitAction, NULL);
    pollingForQuit = false;
    userQuit = false;
}

bool signalPoll::userSaidQuit()
{
    return userQuit;
}

