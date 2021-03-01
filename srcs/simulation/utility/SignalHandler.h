#pragma once
#ifndef SIGNALHANDLER_H_INCLUDED
#define SIGNALHANDLER_H_INCLUDED

#include <signal.h>

extern int g_shouldClose;

void signalHandler(int /** signalNumber **/)
{
    g_shouldClose = 1;
}

#endif // SIGNALHANDLER_H_INCLUDED
