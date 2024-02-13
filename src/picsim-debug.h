/**************************************************************************
 *
 * $Id: picsim-debug.h,v 1.8 2015/11/06 18:01:01 patrick Exp $
 *
 * Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>
 * All Rights Reserved.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2.  of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#ifndef PICSIM_DEBUG_H
#define PICSIM_DEBUG_H


#if !defined(DEBUG_LEVEL)

#define BEGIN_DEBUG_MASTER_OUTPUT(level)
#define END_DEBUG_MASTER_OUTPUT
#define BEGIN_DEBUG_OUTPUT(level)
#define END_DEBUG_OUTPUT
#define HEADER_DEBUG_OUTPUT1(h1)
#define HEADER_DEBUG_OUTPUT2(h1,h2)
#define VAR_DEBUG_OUTPUT1(v1)
#define VAR_DEBUG_OUTPUT2(v1,v2)
#define VAR_DEBUG_OUTPUT3(v1,v2,v3)

#else // defined(DEBUG_LEVEL)

#if defined(HAVE_MPI)

#define BEGIN_DEBUG_MASTER_OUTPUT(level)          \
MPI_Barrier(MPI_COMM_WORLD);                      \
if (rankProc == masterProc &&                     \
		debugLevel() >= level) {                      \

#define END_DEBUG_MASTER_OUTPUT                   \
	<< std::endl;                                   \
}                                                 \


#define BEGIN_DEBUG_OUTPUT(level)                 \
if (debugLevel() >= level) {                      \
	for (int ip=0; ip<nbProc; ++ip) {               \
		MPI_Barrier(MPI_COMM_WORLD);                  \
		if (ip == rankProc) {                         \


#define END_DEBUG_OUTPUT0                         \
    }                                             \
  }                                               \
}                                                 \

#define END_DEBUG_OUTPUT                          \
			<< std::endl;                               \
		}                                             \
	}                                               \
}                                                 \

#define HEADER_DEBUG_OUTPUT1(h1)                  \
std::cout << "Proc " << rankProc << ": "          \
	        << h1 << "(): " << std::endl            \


#define HEADER_DEBUG_OUTPUT2(h1,h2)               \
std::cout << "Proc " << rankProc << ": "          \
	        << h1 << "::"                           \
					<< h2 << "(): " << std::endl            \


#else // !defined(HAVE_MPI)            

#define BEGIN_DEBUG_MASTER_OUTPUT(level)          \
if (debugLevel() >= level) {                      \

#define END_DEBUG_MASTER_OUTPUT                   \
	<< std::endl;                                   \
}                                                 \


#define BEGIN_DEBUG_OUTPUT(level)                 \
if (debugLevel() >= level) {                      \

#define END_DEBUG_OUTPUT0                         \
}                                                 \

#define END_DEBUG_OUTPUT                          \
	<< std::endl;                                   \
}                                                 \

#define HEADER_DEBUG_OUTPUT1(h1)                  \
std::cout << h1 << "(): " << std::endl            \


#define HEADER_DEBUG_OUTPUT2(h1,h2)               \
std::cout << h1 << "::"                           \
	        << h2 << "(): " << std::endl            \


#endif // defined(HAVE_MPI)



#define VAR_DEBUG_OUTPUT1(v1)                     \
<< "     " #v1 "=" << v1 << std::endl             \



#define VAR_DEBUG_OUTPUT2(v1,v2)                  \
<< "     " #v1 "=" << v1                          \
<< "     " #v2 "=" << v2 << std::endl             \



#define VAR_DEBUG_OUTPUT3(v1,v2,v3)               \
<< "     " #v1 "=" << v1                          \
<< "     " #v2 "=" << v2                          \
<< "     " #v3 "=" << v3 << std::endl             \

#endif // DEBUG_LEVEL

#endif // PICSIM_DEFS_H
