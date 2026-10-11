#include <VoronoiMeshMaker/Boundary2D/Builders/Boundary2DBuilder.hpp>
#include <VoronoiMeshMaker/Boundary2D/Shapes/Rectangle.hpp>
#include <VoronoiMeshMaker/Sites2D/Factory/SiteFactory.hpp>
#include <VoronoiMeshMaker/Voronoi2D/Diagram/ClippedVoronoiBuilder2D.hpp>
#include <iostream>
int main(){try{auto b=vmm::b2d::make_boundary(vmm::b2d::Rectangle(vmm::b2d::Point2{-15,-5},30,10)); auto s=vmm::s2d::make_sites(b,vmm::s2d::UniformRandom2D{87552,12345},vmm::s2d::SiteValidationOptions{});auto d=vmm::vd2d::ClippedVoronoiBuilder2D::build(s,b);std::cout<<"accepted "<<d.cells.size()<<"\n";}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
