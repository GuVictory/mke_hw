#ifndef DZ2_EXPORT_H
#define DZ2_EXPORT_H

#include "mke/Triangulation.h"
#include "entities/Vec.h"

#include <ostream>

void outExport(std::ostream& os, const Triangulation& t, const Vec& values);
void fileExport(const Triangulation& t, const Vec& values);

#endif //DZ2_EXPORT_H
