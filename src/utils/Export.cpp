#include "utils/Export.h"
#include <utils/GlobalIndices.h>

#include <tuple>
#include <fstream>

using namespace std;

void outExport(ostream& os, const Triangulation& t, const Vec& values)
{
    auto m = t.nodes().size();
    auto n = t.triangles().size();

    vector<Point3d> nodes(m);
    for (auto& n : t.nodes()) {
        auto [p, i] = n;
        nodes[i] = p;
    }

    os << m << " 3 3 U_x Uy Uz " << endl;
    for (Index i = 0; i < m; ++i) {
        auto& p = nodes[i];
        auto [x, y, z]
        = tie(values[_g(i, Coord::X, m)], values[_g(i, Coord::Y, m)],
              values[_g(i, Coord::Z, m)]);

        os << (i + 1) << " " << p[0] << " " << p[1] << " " << p[2] << " ";
        os << x << " " << y << " " << z << endl;
    }

    os << n << " 3 3 BC_id mat_id mat_id_Out" << endl;
    auto& tr = t.triangles();
    for (Index i = 0; i < n; ++i) {
        auto& [j, e] = tr[i];
        os << (i + 1) << " " << (e[0].index() + 1) << " " << (e[1].index() + 1)
           << " " << (e[2].index() + 1) << " " << (j + 1) << " 1 0" << endl;
    }
}

void fileExport(const Triangulation& t, const Vec& values)
{

    std::ofstream out;
    out.open("result.dat");

    auto m = t.nodes().size();
    auto n = t.triangles().size();

    vector<Point3d> nodes(m);
    for (auto& n : t.nodes()) {
        auto [p, i] = n;
        nodes[i] = p;
    }

    out << m << " 3 3 Ux Uy Uz " << endl;
    for (Index i = 0; i < m; ++i) {
        auto& p = nodes[i];
        auto [x, y, z]
        = tie(values[_g(i, Coord::X, m)], values[_g(i, Coord::Y, m)],
              values[_g(i, Coord::Z, m)]);

        out << (i + 1) << " " << p[0] << " " << p[1] << " " << p[2] << " ";
        out << x << " " << y << " " << z << endl;
    }

    out << n << " 3 3 BC_id mat_id mat_id_Out" << endl;
    auto& tr = t.triangles();
    for (Index i = 0; i < n; ++i) {
        auto& [j, e] = tr[i];
        out << (i + 1) << " " << (e[0].index() + 1) << " " << (e[1].index() + 1)
           << " " << (e[2].index() + 1) << " " << (j + 1) << " 1 0" << endl;
    }
    out.close();
}