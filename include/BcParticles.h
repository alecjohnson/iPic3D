/***************************************************************************
  BcParticles.h  -  Library to manage boundary conditions for particles
  -------------------
begin                : Fri Jun 4 2004
copyright            : (C) 2004 Los Alamos National Laboratory
developers           : Stefano Markidis, Giovanni Lapenta
email                : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/

#ifndef BcParticles_H
#define BcParticles_H

namespace BCparticles
{
  enum Enum
  {
    EXIT = 0,
    PERFECT_MIRROR = 1,
    REEMISSION = 2
  };
}

#endif
