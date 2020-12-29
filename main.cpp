#include <utils/Constants.h>
#include <utils/Export.h>
#include <generators/GlobalEqGen.h>
#include <generators/LocalEqGen.h>
#include <generators/LupFactor.h>
#include <mke/Triangulation.h>

#include <chrono>
#include <fstream>
#include <iostream>

using namespace std;

int main(int argc, char const* argv[])
{
    auto init = consts::init;
    auto t = Triangulation::cuboid({consts::xdim, consts::ydim, consts::zdim}, consts::scale);

    auto start = chrono::high_resolution_clock::now();

    auto [mat, vec] = GlobalEqGenSystem(t, LocalEqGen(gen_local));
    auto x0 = Vec(vec.size(), init);
    auto res = psolve(mat, vec, move(x0));

    auto finish = chrono::high_resolution_clock::now();
    cout << "Времы выполнения:  "
         << chrono::duration_cast<chrono::milliseconds>(finish - start)
                 .count()
         << "ms" << endl;

    t.extract_triangles();

    // Выведем в консоль
    // outExport(cout, t, res);

    // Или выведем в файл
    fileExport(t, res);
    return 0;
}
