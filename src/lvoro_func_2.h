#ifndef LVORO_FUNC_2_H_
#define LVORO_FUNC_2_H_ 1

#include "./lvoro_func.h"

#include <cassert>

using std::string;

using namespace voro;

void readData(string s, vector<LVoronoiPoint>& vps, container& con, const int numparticles) {
  ifstream ifs;
  ifs.open(s.c_str());

  for (int i = 0; i != numparticles; ++i) {
    double x, y, z;
    ifs >> x; ifs >> y; ifs >> z;
    vps[i].x_ = x;
    vps[i].y_ = y;
    vps[i].z_ = z;
    con.put(i, x, y, z);
  }

  ifs.close();

}

/*
 *list = i*numparticles + j
 */
void makeTransitionPairs( container& con, vector<int>& list, const int numparticles) {
  const int listsize = 100;

  c_loop_all cmh(con);
  voronoicell_neighbor ch;

  bool breaktag = false;

  if (cmh.start()) do if (con.compute_cell(ch, cmh)) {
        int id = cmh.pid();
        vector<int> neighbors;
        ch.neighbors(neighbors);

        for (auto ni : neighbors) {
          if ( ni > id) list.push_back(id*numparticles+ni);
          if ( list.size() >= listsize ) {
            breaktag = true;
            break;
          }
        }

        if (breaktag) break;
      } while (cmh.inc());

  assert( list.size() == listsize);
}

bool checkPairness( container& con, const int pair1, const int pair2) {
  assert( pair1 < pair2);

    c_loop_all cmh(con);
  voronoicell_neighbor ch;

  bool breaktag = false;

  if (cmh.start()) do if (con.compute_cell(ch, cmh)) {
        int id = cmh.pid();
        if ( id != pair1 ) continue;

        vector<int> neighbors;
        ch.neighbors(neighbors);
        for (auto ni : neighbors) {
          if ( ni == pair2 ) return true;
        }

        return false;

      } while (cmh.inc());
}

#endif
