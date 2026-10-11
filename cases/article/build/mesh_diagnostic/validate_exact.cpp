#include "ArticleExactClippedVoronoiBuilder2D.hpp"
#include <VoronoiMeshMaker/Boundary2D/Builders/Boundary2DBuilder.hpp>
#include <VoronoiMeshMaker/Boundary2D/Shapes/Rectangle.hpp>
#include <VoronoiMeshMaker/Sites2D/Factory/SiteFactory.hpp>
#include <iostream>
#include <iomanip>
int main(int argc,char** argv){try{auto n=argc>1?std::stoull(argv[1]):87552;auto b=vmm::b2d::make_boundary(vmm::b2d::Rectangle(vmm::b2d::Point2{-15,-5},30,10));auto s=vmm::s2d::make_sites(b,vmm::s2d::UniformRandom2D{n,12345},vmm::s2d::SiteValidationOptions{});auto d=vmm::vd2d::ArticleExactClippedVoronoiBuilder2D::build(s,b);long double area=0;std::size_t internal=0;for(auto& c:d.cells){area+=c.signed_area();for(auto& e:c.edges)if(!e.is_boundary_edge){++internal;if(e.neighbour_site_id.value<0)throw std::runtime_error("missing neighbour");}}std::cout<<std::setprecision(17)<<"ACCEPTED cells="<<d.cells.size()<<" internal_directed_faces="<<internal<<" area="<<area<<" error="<<std::abs(area-300)<<std::endl;if(std::abs(area-300)>1e-8)throw std::runtime_error("area mismatch");}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
