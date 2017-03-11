#ifndef VORO_FUNC_H_
#define VORO_FUNC_H_ 1

#include "../voro++/voro++.hh"
#include "./spaceutility.h"
using namespace voro;

double retTotalEnergy(container& con, double value_tarea) {

  double totalenergy = 0.0;

  c_loop_all cm(con);
  voronoicell_neighbor cl;
  if(cm.start()) do if(con.compute_cell(cl,cm)) {
        double ve = (cl.volume()-1.0);
        double ae = (cl.surface_area()-value_tarea);
        totalenergy += ve*ve+ae*ae;
      } while (cm.inc());

  return totalenergy;
}

#endif
