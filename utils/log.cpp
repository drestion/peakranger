/*
 * log.cpp
 *
 *  Created on: Jun 18, 2011
 *      Author: xin
 */
#include "log.h"
using namespace std;
FILE* Output2FILE::pStream = stderr;
string Output2FILE::filename = "stderr";
bool Output2FILE::stream_is_open = false;
void Output2FILE::closeStream(){
  if(pStream!=stderr && pStream!=stdout && stream_is_open)
  fclose(pStream);
}

