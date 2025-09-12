#include <gmsh.h>
#include <vector>
#include <utility>

void assign_physical(const std::vector<std::pair<int,int>>& out, int volTag, int lateralTag, int portInTag, int portOutTag,
                     const std::string& volName, const std::string& lateralName, const std::string& portInName, const std::string& portOutName) {
    int vol = out[1].second;
    int lat1 = out[2].second, port1 = out[3].second;
    int lat2 = out[4].second, port2 = out[5].second;

    int pgVol = gmsh::model::addPhysicalGroup(3, {vol}, volTag);
    gmsh::model::setPhysicalName(3, pgVol, volName);
    gmsh::model::addPhysicalGroup(2, {lat1, lat2}, lateralTag);
    gmsh::model::setPhysicalName(2, lateralTag, lateralName);
    gmsh::model::addPhysicalGroup(2, {port1}, portInTag);
    gmsh::model::setPhysicalName(2, portInTag, portInName);
    gmsh::model::addPhysicalGroup(2, {port2}, portOutTag);
    gmsh::model::setPhysicalName(2, portOutTag, portOutName);
}

int main() {
    gmsh::initialize();
    gmsh::model::add("model");

    double lc = 1e-6;
    std::vector<int> pts(4);
    pts[0] = gmsh::model::occ::addPoint(0, -1000, 0, lc);
    pts[1] = gmsh::model::occ::addPoint(20000, -1000, 0, lc);
    pts[2] = gmsh::model::occ::addPoint(20000, 1000, 0, lc);
    pts[3] = gmsh::model::occ::addPoint(0, 1000, 0, lc);

    std::vector<int> lines(4);
    for (int i = 0; i < 4; ++i)
        lines[i] = gmsh::model::occ::addLine(pts[i], pts[(i+1)%4]);
    int loop = gmsh::model::occ::addCurveLoop(lines);
    int surf_kapton = gmsh::model::occ::addPlaneSurface({loop});

    // Signal
    std::vector<int> pSig(4);
    double coordsSig[4][2] = {{0,0},{20000,0},{20000,100},{0,100}};
    for (int i = 0; i < 4; ++i)
        pSig[i] = gmsh::model::occ::addPoint(coordsSig[i][0], coordsSig[i][1], 0, lc);
    std::vector<int> lSig(4);
    for (int i = 0; i < 4; ++i)
        lSig[i] = gmsh::model::occ::addLine(pSig[i], pSig[(i+1)%4]);
    int loopSig = gmsh::model::occ::addCurveLoop(lSig);
    int surfSig = gmsh::model::occ::addPlaneSurface({loopSig});

    // GND1
    std::vector<int> pG1(4), lG1(4);
    double coordsG1[4][2] = {{0,-200},{20000,-200},{20000,-100},{0,-100}};
    for (int i = 0; i < 4; ++i) pG1[i] = gmsh::model::occ::addPoint(coordsG1[i][0], coordsG1[i][1], 0, lc);
    for (int i = 0; i < 4; ++i) lG1[i] = gmsh::model::occ::addLine(pG1[i], pG1[(i+1)%4]);
    int loopG1 = gmsh::model::occ::addCurveLoop(lG1);
    int surfG1 = gmsh::model::occ::addPlaneSurface({loopG1});

    // GND2
    std::vector<int> pG2(4), lG2(4);
    double coordsG2[4][2] = {{0,200},{20000,200},{20000,300},{0,300}};
    for (int i = 0; i < 4; ++i) pG2[i] = gmsh::model::occ::addPoint(coordsG2[i][0], coordsG2[i][1], 0, lc);
    for (int i = 0; i < 4; ++i) lG2[i] = gmsh::model::occ::addLine(pG2[i], pG2[(i+1)%4]);
    int loopG2 = gmsh::model::occ::addCurveLoop(lG2);
    int surfG2 = gmsh::model::occ::addPlaneSurface({loopG2});

    gmsh::model::occ::synchronize();

    // Extrusions
    auto outKapton = gmsh::model::occ::extrude({{2, surf_kapton}}, 0, 0, -100);
    auto outAir = gmsh::model::occ::extrude({{2, surf_kapton}}, 0, 0, 100);
    auto outSig = gmsh::model::occ::extrude({{2, surfSig}}, 0, 0, 20);
    auto outG1 = gmsh::model::occ::extrude({{2, surfG1}}, 0, 0, 20);
    auto outG2 = gmsh::model::occ::extrude({{2, surfG2}}, 0, 0, 20);

    gmsh::model::occ::synchronize();

    // Assign physical groups
    assign_physical(outKapton, 3001, 0, 0, 0, "kapton", "", "", "");
    assign_physical(outAir, 3002, 0, 0, 0, "air", "", "", "");
    assign_physical(outSig, 3003, 2001, 2004, 2005, "signal", "signal_lateral", "signal_port_in", "signal_port_out");
    assign_physical(outG1, 3004, 2002, 2006, 2007, "gnd1", "gnd1_lateral", "gnd1_port_in", "gnd1_port_out");
    assign_physical(outG2, 3005, 2003, 2008, 2009, "gnd2", "gnd2_lateral", "gnd2_port_in", "gnd2_port_out");

    gmsh::model::mesh::generate(3);
    gmsh::write("model.msh");
    gmsh::finalize();
    return 0;
}

