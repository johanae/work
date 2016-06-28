/*
 *  distribute.h
 *  
 *
 *  Created by Johannes Langguth on 26.08.14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

void sdd_Retreive(ELLmatrixVector* Z,ELLmatrixVector* LZ,localparttype* L);
void sddMPI_localRenumbering(ELLmatrixVector* LZ, localparttype* L);
void mainloop(int TIME_STEP,ELLmatrixVector* Z);